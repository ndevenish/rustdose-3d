use crate::parser::config::WedgeConfig;

/// Wedge describes a single rotation wedge of an experiment.
/// All angles are stored internally in radians.
/// All positions/translations are in micrometres.
#[derive(Debug, Clone)]
pub struct Wedge {
    /// Start angle in radians.
    pub start_ang: f64,
    /// End angle in radians.
    pub end_ang: f64,
    /// Angular resolution in radians per image.
    pub ang_res: f64,
    /// Total exposure time in seconds.
    pub exposure_time: f64,
    /// Start position offset X in µm.
    pub start_x: f64,
    /// Start position offset Y in µm.
    pub start_y: f64,
    /// Start position offset Z in µm.
    pub start_z: f64,
    /// Translation on X axis per radian in µm.
    pub trans_x: f64,
    /// Translation on Y axis per radian in µm.
    pub trans_y: f64,
    /// Translation on Z axis per radian in µm.
    pub trans_z: f64,
    /// Off-axis rotation distance in µm.
    pub off_axis_um: f64,
    /// Maximum resolution in Å.
    pub max_resolution: f64,
}

impl Wedge {
    /// Create a Wedge from parsed configuration.
    /// Input angles in WedgeConfig are in degrees; we convert to radians.
    pub fn from_config(config: &WedgeConfig) -> Self {
        let start_ang = config.start_ang.to_radians();
        let end_ang = config.end_ang.to_radians();
        let exposure_time = config.exposure_time.unwrap_or(0.0);

        // Angular resolution: default 2 degrees, but check it's not too coarse
        let ang_res_deg = config.angular_resolution.unwrap_or(2.0);
        let total_range = config.end_ang - config.start_ang;
        let ang_res = if ang_res_deg > 0.0 && (total_range / ang_res_deg) > 10.0 {
            ang_res_deg.to_radians()
        } else {
            2.0_f64.to_radians()
        };

        // Translations: input is per degree, stored as per radian
        // Converting "per degree" to "per radian" means multiplying by (180/π)
        let trans_x = config.translate_x.unwrap_or(0.0).to_degrees(); // toDegrees = * 180/π
        let trans_y = config.translate_y.unwrap_or(0.0).to_degrees();
        let trans_z = config.translate_z.unwrap_or(0.0).to_degrees();

        Wedge {
            start_ang,
            end_ang,
            ang_res,
            exposure_time,
            start_x: config.start_offset_x.unwrap_or(0.0),
            start_y: config.start_offset_y.unwrap_or(0.0),
            start_z: config.start_offset_z.unwrap_or(0.0),
            trans_x,
            trans_y,
            trans_z,
            off_axis_um: config.rot_ax_beam_offset.unwrap_or(0.0),
            max_resolution: config.max_resolution.unwrap_or(2.0),
        }
    }

    /// Returns the start position vector [x, y, z] in µm.
    pub fn start_vector(&self) -> [f64; 3] {
        [self.start_x, self.start_y, self.start_z]
    }

    /// Returns the translation vector [x, y, z] at a given angle (radians) in µm.
    pub fn translation_vector(&self, delta_phi: f64) -> [f64; 3] {
        let delta = delta_phi - self.start_ang;
        [
            self.trans_x * delta,
            self.trans_y * delta,
            self.trans_z * delta,
        ]
    }

    /// Returns total seconds of exposure.
    pub fn total_sec(&self) -> f64 {
        self.exposure_time
    }

    /// Returns a human-readable description of the wedge.
    pub fn description(&self) -> String {
        let mut s = format!(
            "Collecting data for a total of {:.1}s from phi = {:.1} to {:.1} deg.\n",
            self.exposure_time,
            self.start_ang.to_degrees(),
            self.end_ang.to_degrees()
        );

        if self.start_x != 0.0
            || self.start_y != 0.0
            || self.start_z != 0.0
            || self.trans_x != 0.0
            || self.trans_y != 0.0
            || self.trans_z != 0.0
        {
            s += &format!(
                "Start is offset by [{}, {}, {}] um [x,y,z].\n\
                 Helical scanning is at [{}, {}, {}] um/deg in [x,y,z]\n",
                self.start_x,
                self.start_y,
                self.start_z,
                self.trans_x.to_radians(), // convert back to per-degree for display
                self.trans_y.to_radians(),
                self.trans_z.to_radians()
            );
        }

        s
    }
}
