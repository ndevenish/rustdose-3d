use crate::constants::KEV_TO_JOULES;
use crate::container::Container;
use crate::parser::config::BeamConfig;

/// BeamExperimental: 2D intensity grid with bilinear interpolation.
///
/// The grid is provided as a flat `Vec<Vec<f64>>` (rows × cols) representing
/// the measured beam profile. The grid is padded with a zero border and
/// normalised so that the total flux equals `photons_per_sec`.
///
/// `beam_intensity` uses bilinear interpolation to return the intensity
/// at any continuous (x, y) coordinate.
#[derive(Debug)]
pub struct BeamExperimental {
    /// Photon energy in keV.
    photon_energy: f64,
    /// Total (un-attenuated) flux in photons/s.
    total_flux: f64,
    /// Flux after container attenuation.
    attenuated_flux: f64,
    /// Pixel size X in µm.
    pix_x: f64,
    /// Pixel size Y in µm.
    pix_y: f64,
    /// Total beam width (incl. border) in µm.
    beam_x_size: f64,
    /// Total beam height (incl. border) in µm.
    beam_y_size: f64,
    /// Pulse energy in mJ.
    pulse_energy: f64,
    /// Energy FWHM for pink beam (None if monochromatic).
    energy_fwhm: Option<f64>,
    /// Whether circular (elliptical) collimation is applied.
    is_circular: bool,
    /// Original data grid (rows × cols).
    data: Vec<Vec<f64>>,
    /// Processed beam array (data + zero border, normalised × energy × flux / pixel).
    /// Populated by `generate_beam_array`.
    beam_array: Vec<Vec<f64>>,
    /// Whether `generate_beam_array` has been called.
    array_ready: bool,
}

impl BeamExperimental {
    /// Construct from config. The `data` grid must be set separately via
    /// `set_data` before `generate_beam_array` is called.
    pub fn from_config(config: &BeamConfig) -> Result<Self, String> {
        let photon_energy = config.energy.ok_or("Experimental beam requires energy")?;
        let total_flux    = config.flux.ok_or("Experimental beam requires flux")?;
        let pix_x = config.pixel_size_x.ok_or("Experimental beam requires pixel_size_x")?;
        let pix_y = config.pixel_size_y.ok_or("Experimental beam requires pixel_size_y")?;

        // Collimation: if circular, set flag
        use crate::parser::config::Collimation;
        let is_circular = matches!(config.collimation, Some(Collimation::Circular { .. }));

        Ok(BeamExperimental {
            photon_energy,
            total_flux,
            attenuated_flux: total_flux,
            pix_x,
            pix_y,
            beam_x_size: 0.0,
            beam_y_size: 0.0,
            pulse_energy: config.pulse_energy.unwrap_or(0.0),
            energy_fwhm: config.energy_fwhm,
            is_circular,
            data: Vec::new(),
            beam_array: Vec::new(),
            array_ready: false,
        })
    }

    /// Set the raw intensity grid (rows = vertical, cols = horizontal).
    pub fn set_data(&mut self, data: Vec<Vec<f64>>) {
        let ncols = data.first().map(|r| r.len()).unwrap_or(0);
        let nrows = data.len();
        self.beam_x_size = (ncols + 2) as f64 * self.pix_x;
        self.beam_y_size = (nrows + 2) as f64 * self.pix_y;
        self.data = data;
        self.array_ready = false;
    }

    fn bilinear_interpolate(v00: f64, v10: f64, v01: f64, v11: f64, x: f64, y: f64) -> f64 {
        v00 * (1.0 - x) * (1.0 - y)
            + v10 * x * (1.0 - y)
            + v01 * (1.0 - x) * y
            + v11 * x * y
    }

    fn do_bilinear(&self, coord_x: f64, coord_y: f64, off_axis_um: f64) -> f64 {
        let real_x = coord_x - off_axis_um + self.beam_x_size / 2.0;
        let real_y = coord_y + self.beam_y_size / 2.0;
        let vh = (real_x / self.pix_x - 0.5).floor() as isize;
        let vv = (real_y / self.pix_y - 0.5).floor() as isize;

        let rows = self.beam_array.len() as isize;
        let cols = self.beam_array.first().map(|r| r.len() as isize).unwrap_or(0);

        if vh >= 0 && vv >= 0 && vh < cols - 1 && vv < rows - 1 {
            let frac_x = real_x / self.pix_x - (vh as f64 + 0.5);
            let frac_y = real_y / self.pix_y - (vv as f64 + 0.5);
            let (vv, vh) = (vv as usize, vh as usize);
            Self::bilinear_interpolate(
                self.beam_array[vv][vh],
                self.beam_array[vv][vh + 1],
                self.beam_array[vv + 1][vh],
                self.beam_array[vv + 1][vh + 1],
                frac_x,
                frac_y,
            )
        } else if vh >= 0 && vv >= 0 && (vh == cols - 1 || vv == rows - 1) {
            let (vv, vh) = (vv.min(rows - 1) as usize, vh.min(cols - 1) as usize);
            self.beam_array[vv][vh]
        } else {
            0.0
        }
    }
}

impl super::Beam for BeamExperimental {
    fn generate_beam_array(&mut self) {
        if self.data.is_empty() {
            self.array_ready = true;
            return;
        }
        let nrows = self.data.len();
        let ncols = self.data[0].len();
        let mut beam_sum = 0.0;
        // Build padded array (nrows+2 × ncols+2)
        let mut arr = vec![vec![0.0f64; ncols + 2]; nrows + 2];
        for i in 1..=nrows {
            for j in 1..=ncols {
                let v = self.data[i - 1][j - 1];
                arr[i][j] = v;
                beam_sum += v;
            }
        }
        // Normalise: divide by beam_sum, scale by energy × flux / pixel area
        let scale = if beam_sum > 0.0 {
            KEV_TO_JOULES * self.photon_energy * self.attenuated_flux
                / (beam_sum * self.pix_x * self.pix_y)
        } else {
            0.0
        };
        for row in arr.iter_mut() {
            for v in row.iter_mut() {
                *v *= scale;
            }
        }
        self.beam_array = arr;
        self.array_ready = true;
    }

    fn beam_intensity(&self, coord_x: f64, coord_y: f64, off_axis_um: f64) -> f64 {
        if !self.array_ready {
            return 0.0;
        }
        let hx = self.beam_x_size / 2.0 - self.pix_x;
        let hy = self.beam_y_size / 2.0 - self.pix_y;
        if self.is_circular {
            let dx = coord_x - off_axis_um;
            if (dx / hx).powi(2) + (coord_y / hy).powi(2) > 1.0 {
                return 0.0;
            }
        } else {
            if (coord_x - off_axis_um).abs() > hx || coord_y.abs() > hy {
                return 0.0;
            }
        }
        self.do_bilinear(coord_x, coord_y, off_axis_um)
    }

    fn description(&self) -> String {
        format!(
            "Experimental beam: {:.2e} photons/s, {:.2} keV, grid {}×{} pixels",
            self.total_flux,
            self.photon_energy,
            self.data.len(),
            self.data.first().map(|r| r.len()).unwrap_or(0),
        )
    }

    fn photons_per_sec(&self) -> f64 { self.attenuated_flux }
    fn photon_energy(&self) -> f64 { self.photon_energy }
    fn pulse_energy(&self) -> f64 { self.pulse_energy }
    fn energy_fwhm(&self) -> Option<f64> { self.energy_fwhm }
    fn is_circular(&self) -> bool { self.is_circular }

    fn apply_container_attenuation(&mut self, container: &dyn Container) {
        let frac = container.attenuation_fraction();
        self.attenuated_flux = self.total_flux * (1.0 - frac);
        if container.material_name().is_some() {
            eprintln!(
                "Beam photons per second after container attenuation is {:.2e} photons per second",
                self.attenuated_flux
            );
        }
        // Regenerate beam array with new flux
        self.array_ready = false;
        self.generate_beam_array();
    }

    fn beam_minimum_dimension(&self) -> f64 {
        self.beam_x_size.min(self.beam_y_size)
    }

    fn beam_area(&self) -> f64 {
        // Approximate by pixel count × pixel area
        let ncols = self.data.first().map(|r| r.len()).unwrap_or(0);
        self.data.len() as f64 * ncols as f64 * self.pix_x * self.pix_y
    }

    fn beam_x(&self) -> Option<f64> { Some(self.beam_x_size) }
    fn beam_y(&self) -> Option<f64> { Some(self.beam_y_size) }
    fn beam_type(&self) -> &str { "Experimental" }

    fn sx(&self) -> f64 { self.beam_x_size / (2.0 * std::f64::consts::SQRT_2) }
    fn sy(&self) -> f64 { self.beam_y_size / (2.0 * std::f64::consts::SQRT_2) }
}
