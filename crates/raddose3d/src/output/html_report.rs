use std::io::Write;

use crate::beam::Beam;
use crate::crystal::Crystal;
use crate::output::{DebugWriter, ExposureSummary};
use crate::wedge::Wedge;

const ABS_EN_THRESHOLD: f64 = 0.95;
/// Maximum voxels per axis before downsampling (keeps HTML file size reasonable).
const MAX_VOXELS_PER_AXIS: usize = 91;

/// Per-wedge summary snapshot captured during the experiment.
#[derive(Debug)]
struct WedgeSummary {
    avg_dwd: f64,
    last_dwd: f64,
    elastic_yield: f64,
    diffraction_efficiency: f64,
    avg_dose_whole: f64,
    avg_dose_exposed: f64,
    max_dose: f64,
    dose_threshold: f64,
    avg_dose_threshold: f64,
    dose_contrast: f64,
    used_volume: f64,
    absorbed_energy: f64,
    dose_inefficiency: f64,
    dose_inefficiency_pe: f64,
    /// Normalised histogram fractions (11 bins).
    dose_hist: [f64; 11],
}

/// Outputs a self-contained HTML report with Plotly.js visualizations.
///
/// Includes:
/// - 3D isosurface plot of the dose distribution
/// - Summary statistics table (per wedge)
/// - Dose histogram bar chart
#[derive(Debug)]
pub struct OutputDoseStateHTML {
    writer: DebugWriter,
    crystal_info: String,
    ddm_name: String,
    beam_desc: String,
    wedge_summaries: Vec<WedgeSummary>,
}

impl OutputDoseStateHTML {
    pub fn new(writer: Box<dyn Write + Send + Sync>) -> Self {
        OutputDoseStateHTML {
            writer: DebugWriter(writer),
            crystal_info: String::new(),
            ddm_name: String::new(),
            beam_desc: String::new(),
            wedge_summaries: Vec::new(),
        }
    }

    fn write_html(&mut self, crystal: &dyn Crystal) {
        let size = crystal.cryst_size_voxels();
        let size_um = crystal.cryst_size_um();

        // Compute downsampling steps if grid is large.
        let step = [
            size[0].div_ceil(MAX_VOXELS_PER_AXIS),
            size[1].div_ceil(MAX_VOXELS_PER_AXIS),
            size[2].div_ceil(MAX_VOXELS_PER_AXIS),
        ];
        let nx = size[0].div_ceil(step[0]);
        let ny = size[1].div_ceil(step[1]);
        let nz = size[2].div_ceil(step[2]);

        // Build flat arrays for Plotly isosurface: x, y, z, value on the full regular grid.
        // Non-crystal voxels get dose 0.
        let total = nx * ny * nz;
        let mut xs = Vec::with_capacity(total);
        let mut ys = Vec::with_capacity(total);
        let mut zs = Vec::with_capacity(total);
        let mut vals = Vec::with_capacity(total);
        let mut max_dose: f64 = 0.0;

        for ii in 0..nx {
            let i = (ii * step[0]).min(size[0] - 1);
            let x = if nx > 1 {
                size_um[0] * ii as f64 / (nx - 1) as f64
            } else {
                0.0
            };
            for jj in 0..ny {
                let j = (jj * step[1]).min(size[1] - 1);
                let y = if ny > 1 {
                    size_um[1] * jj as f64 / (ny - 1) as f64
                } else {
                    0.0
                };
                for kk in 0..nz {
                    let k = (kk * step[2]).min(size[2] - 1);
                    let z = if nz > 1 {
                        size_um[2] * kk as f64 / (nz - 1) as f64
                    } else {
                        0.0
                    };
                    xs.push(x);
                    ys.push(y);
                    zs.push(z);
                    let dose = if crystal.is_crystal_at(i, j, k) {
                        crystal.get_dose(i, j, k)
                    } else {
                        0.0
                    };
                    if dose > max_dose {
                        max_dose = dose;
                    }
                    vals.push(dose);
                }
            }
        }

        // Compute sensible isosurface levels based on actual dose range.
        let iso_levels = compute_iso_levels(max_dose);

        // Use the last wedge summary for the histogram (cumulative).
        let last_summary = self.wedge_summaries.last();

        // --- Write HTML ---
        let _ = self.writer.write_all(b"<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n\
            <meta charset=\"UTF-8\">\n\
            <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n\
            <title>RADDOSE-3D Dose Report</title>\n\
            <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\" charset=\"utf-8\"></script>\n\
            <style>\n\
            * { box-sizing: border-box; margin: 0; padding: 0; }\n\
            body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;\n\
                   background: #f5f5f5; color: #333; padding: 20px; }\n\
            h1 { margin-bottom: 4px; }\n\
            .subtitle { color: #666; margin-bottom: 20px; font-size: 0.95em; }\n\
            .card { background: #fff; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.12);\n\
                    padding: 20px; margin-bottom: 20px; }\n\
            .card h2 { margin-bottom: 12px; font-size: 1.2em; }\n\
            .layout { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }\n\
            .right-col { display: flex; flex-direction: column; gap: 20px; }\n\
            @media (max-width: 1000px) { .layout { grid-template-columns: 1fr; } }\n\
            table { border-collapse: collapse; width: 100%; }\n\
            th, td { text-align: left; padding: 8px 12px; border-bottom: 1px solid #eee; }\n\
            th { background: #f8f8f8; font-weight: 600; white-space: nowrap; }\n\
            td.num { text-align: right; font-variant-numeric: tabular-nums; }\n\
            #dose3d { width: 100%; height: 500px; }\n\
            #histogram { width: 100%; height: 300px; }\n\
            .info { font-size: 0.9em; color: #666; margin-top: 8px; }\n\
            </style>\n</head>\n<body>\n");

        let _ = writeln!(self.writer, "<h1>RADDOSE-3D Dose Report</h1>");
        let _ = write!(self.writer, "<p class=\"subtitle\">");
        // Crystal info
        let crystal_info_escaped = html_escape(&self.crystal_info);
        let _ = write!(self.writer, "{}", crystal_info_escaped);
        if !self.ddm_name.is_empty() {
            let _ = write!(self.writer, " &mdash; {}", html_escape(&self.ddm_name));
        }
        if !self.beam_desc.is_empty() {
            let _ = write!(self.writer, "<br>{}", html_escape(&self.beam_desc));
        }
        let _ = writeln!(self.writer, "</p>");

        // --- Two-column layout: table (left), 3D + histogram (right) ---
        let _ = writeln!(self.writer, "<div class=\"layout\">");

        // Left column: summary table
        let _ = writeln!(self.writer, "<div class=\"card\">");
        let _ = writeln!(self.writer, "<h2>Summary Statistics</h2>");
        let _ = writeln!(self.writer, "<table>");
        let _ = writeln!(
            self.writer,
            "<thead><tr><th>Metric</th><th>Value</th></tr></thead><tbody>"
        );
        if let Some(s) = last_summary {
            write_table_row(
                &mut self.writer,
                "Avg Diffraction Weighted Dose",
                &format!("{:.6} MGy", s.avg_dwd),
            );
            write_table_row(
                &mut self.writer,
                "Last DWD",
                &format!("{:.6} MGy", s.last_dwd),
            );
            write_table_row_html(
                &mut self.writer,
                "Elastic Yield",
                &format!("{} photons", html_scientific(s.elastic_yield, 2)),
            );
            write_table_row_html(
                &mut self.writer,
                "Diffraction Efficiency",
                &format!(
                    "{} photons/MGy",
                    html_scientific(s.diffraction_efficiency, 2)
                ),
            );
            write_table_row(
                &mut self.writer,
                "Avg Dose (Whole Crystal)",
                &format!("{:.6} MGy", s.avg_dose_whole),
            );
            write_table_row(
                &mut self.writer,
                "Avg Dose (Exposed Region)",
                &format!("{:.6} MGy", s.avg_dose_exposed),
            );
            write_table_row(
                &mut self.writer,
                "Max Dose",
                &format!("{:.6} MGy", s.max_dose),
            );
            write_table_row(
                &mut self.writer,
                &format!("Dose Threshold ({:.0}%)", ABS_EN_THRESHOLD * 100.0),
                &format!("{:.6} MGy", s.dose_threshold),
            );
            write_table_row(
                &mut self.writer,
                "Avg Dose (Threshold)",
                &format!("{:.6} MGy", s.avg_dose_threshold),
            );
            write_table_row(
                &mut self.writer,
                "Dose Contrast",
                &format!("{:.2}", s.dose_contrast),
            );
            write_table_row(
                &mut self.writer,
                "Used Volume",
                &format!("{:.1}%", s.used_volume),
            );
            write_table_row_html(
                &mut self.writer,
                "Absorbed Energy",
                &format!("{} J", html_scientific(s.absorbed_energy, 2)),
            );
            write_table_row(
                &mut self.writer,
                "Dose Inefficiency",
                &format!("{:.1} 1/g", s.dose_inefficiency),
            );
            write_table_row(
                &mut self.writer,
                "Dose Inefficiency PE",
                &format!("{:.1} 1/g", s.dose_inefficiency_pe),
            );
        }
        let _ = writeln!(self.writer, "</tbody></table></div>");

        // Right column: 3D plot + histogram stacked vertically
        let _ = writeln!(self.writer, "<div class=\"right-col\">");

        // Isosurface plot
        let _ = writeln!(self.writer, "<div class=\"card\">");
        let _ = writeln!(self.writer, "<h2>3D Dose Distribution</h2>");
        let _ = writeln!(self.writer, "<div id=\"dose3d\"></div>");
        let _ = writeln!(
            self.writer,
            "<p class=\"info\">Grid: {}x{}x{} voxels ({}x{}x{} original). \
             Crystal: {:.1} x {:.1} x {:.1} &micro;m. \
             Max dose: {:.4} MGy.</p>",
            nx, ny, nz, size[0], size[1], size[2], size_um[0], size_um[1], size_um[2], max_dose,
        );
        let _ = writeln!(self.writer, "</div>");

        // Histogram
        let _ = writeln!(self.writer, "<div class=\"card\">");
        let _ = writeln!(self.writer, "<h2>Dose Histogram</h2>");
        let _ = writeln!(self.writer, "<div id=\"histogram\"></div>");
        let _ = writeln!(self.writer, "</div>");

        let _ = writeln!(self.writer, "</div>"); // close .right-col
        let _ = writeln!(self.writer, "</div>"); // close .layout

        // --- JavaScript ---
        let _ = writeln!(self.writer, "<script>");

        // Isosurface data
        let _ = write!(self.writer, "var xData=");
        write_json_f64_array(&mut self.writer, &xs);
        let _ = writeln!(self.writer, ";");
        let _ = write!(self.writer, "var yData=");
        write_json_f64_array(&mut self.writer, &ys);
        let _ = writeln!(self.writer, ";");
        let _ = write!(self.writer, "var zData=");
        write_json_f64_array(&mut self.writer, &zs);
        let _ = writeln!(self.writer, ";");
        let _ = write!(self.writer, "var vData=");
        write_json_f64_array(&mut self.writer, &vals);
        let _ = writeln!(self.writer, ";");

        // Create isosurface traces — one per level for independent opacity/color control.
        let _ = writeln!(self.writer, "var traces=[];");
        for (i, level) in iso_levels.iter().enumerate() {
            let _ = writeln!(
                self.writer,
                "traces.push({{type:'isosurface',x:xData,y:yData,z:zData,value:vData,\
                 isomin:{val},isomax:{val},\
                 surface:{{show:true,count:1}},caps:{{x:{{show:false}},y:{{show:false}},z:{{show:false}}}},\
                 colorscale:[[0,'{color}'],[1,'{color}']],showscale:{showscale},\
                 opacity:{opacity},name:'{val:.2} MGy'}});",
                val = level.value,
                color = level.color,
                opacity = level.opacity,
                showscale = if i == 0 { "true" } else { "false" },
            );
        }

        let _ = writeln!(
            self.writer,
            "Plotly.newPlot('dose3d',traces,{{\
             scene:{{xaxis:{{title:'X (&micro;m)'}},yaxis:{{title:'Y (&micro;m)'}},zaxis:{{title:'Z (&micro;m)'}},\
             aspectmode:'data'}},\
             margin:{{l:0,r:0,t:30,b:0}},\
             title:'Dose Isosurfaces (MGy)'\
             }},{{responsive:true}});"
        );

        // Histogram
        if let Some(s) = self.wedge_summaries.last() {
            let hist = &s.dose_hist;
            let bins = hist.len(); // 11
            let _ = write!(self.writer, "var histLabels=[");
            // Bin labels
            let _ = write!(self.writer, "'<0.1'");
            for i in 1..bins - 1 {
                let lo = ExposureSummary::dose_hist_break(i);
                let hi = ExposureSummary::dose_hist_break(i + 1);
                let _ = write!(self.writer, ",'{:.1}-{:.1}'", lo, hi);
            }
            let _ = write!(
                self.writer,
                ",'>={:.1}'",
                ExposureSummary::dose_hist_break(bins - 1)
            );
            let _ = writeln!(self.writer, "];");

            let _ = write!(self.writer, "var histValues=[");
            for (i, v) in hist.iter().enumerate() {
                if i > 0 {
                    let _ = write!(self.writer, ",");
                }
                let _ = write!(self.writer, "{:.4}", v * 100.0);
            }
            let _ = writeln!(self.writer, "];");

            let _ = writeln!(
                self.writer,
                "Plotly.newPlot('histogram',[{{x:histLabels,y:histValues,type:'bar',\
                 marker:{{color:'steelblue'}}}}],\
                 {{xaxis:{{title:'Dose Range (MGy)'}},yaxis:{{title:'Voxels (%)'}},\
                 margin:{{l:50,r:20,t:30,b:60}},title:'Dose Distribution'}},\
                 {{responsive:true}});"
            );
        }

        let _ = writeln!(self.writer, "</script>");
        let _ = writeln!(self.writer, "</body>\n</html>");
    }
}

impl super::Output for OutputDoseStateHTML {
    fn publish_crystal(&mut self, crystal: &dyn Crystal) {
        self.crystal_info = crystal.crystal_info();
        self.ddm_name = crystal.ddm().name().to_string();
    }

    fn publish_beam(&mut self, beam: &dyn Beam) {
        self.beam_desc = beam.description();
    }

    fn publish_wedge(
        &mut self,
        _wedge: &Wedge,
        summary: &ExposureSummary,
        _crystal: Option<&dyn Crystal>,
    ) {
        let dwd = summary.avg_diffracted_dose();
        let diff_eff = if dwd > 0.0 {
            summary.wedge_elastic() / dwd
        } else {
            0.0
        };
        self.wedge_summaries.push(WedgeSummary {
            avg_dwd: dwd,
            last_dwd: summary.last_dwd(),
            elastic_yield: summary.wedge_elastic(),
            diffraction_efficiency: diff_eff,
            avg_dose_whole: summary.avg_dose_whole_crystal(),
            avg_dose_exposed: summary.avg_dose_exposed_region(),
            max_dose: summary.max_dose(),
            dose_threshold: summary.abs_dose_threshold(ABS_EN_THRESHOLD),
            avg_dose_threshold: summary.avg_dose_threshold(ABS_EN_THRESHOLD),
            dose_contrast: summary.dose_contrast(ABS_EN_THRESHOLD),
            used_volume: summary.used_volume_fraction(),
            absorbed_energy: summary.abs_energy_total(),
            dose_inefficiency: summary.dose_inefficiency(),
            dose_inefficiency_pe: summary.dose_inefficiency_pe(),
            dose_hist: summary.dose_hist_normalised(),
        });
    }

    fn close(&mut self, crystal: Option<&dyn Crystal>) {
        match crystal {
            Some(c) => self.write_html(c),
            None => {
                let _ = writeln!(
                    self.writer,
                    "<!-- OutputDoseStateHTML: No crystal object has been seen. -->"
                );
            }
        }
        let _ = self.writer.flush();
    }
}

/// Isosurface level configuration.
struct IsoLevel {
    value: f64,
    color: &'static str,
    opacity: f64,
}

/// Compute sensible isosurface levels based on the actual max dose.
fn compute_iso_levels(max_dose: f64) -> Vec<IsoLevel> {
    if max_dose <= 0.0 {
        return vec![];
    }
    // Use 3 levels: ~10%, ~50%, ~85% of max dose.
    // But also honour the classic 0.1 / 20 / 30 MGy levels if they fit.
    let classic = [
        IsoLevel {
            value: 0.1,
            color: "lightblue",
            opacity: 0.2,
        },
        IsoLevel {
            value: 20.0,
            color: "royalblue",
            opacity: 0.5,
        },
        IsoLevel {
            value: 30.0,
            color: "red",
            opacity: 0.9,
        },
    ];

    // If max dose > 30, the classic levels all fit.
    if max_dose > 30.0 {
        return classic.into_iter().collect();
    }

    // Otherwise pick levels at fractions of max dose.
    vec![
        IsoLevel {
            value: (max_dose * 0.10).max(1e-6),
            color: "lightblue",
            opacity: 0.2,
        },
        IsoLevel {
            value: max_dose * 0.50,
            color: "royalblue",
            opacity: 0.5,
        },
        IsoLevel {
            value: max_dose * 0.85,
            color: "red",
            opacity: 0.9,
        },
    ]
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn write_table_row(w: &mut DebugWriter, label: &str, value: &str) {
    let _ = writeln!(
        w,
        "<tr><td>{}</td><td class=\"num\">{}</td></tr>",
        html_escape(label),
        html_escape(value),
    );
}

/// Write a table row where `value_html` already contains HTML (e.g. superscripts).
fn write_table_row_html(w: &mut DebugWriter, label: &str, value_html: &str) {
    let _ = writeln!(
        w,
        "<tr><td>{}</td><td class=\"num\">{}</td></tr>",
        html_escape(label),
        value_html,
    );
}

/// Format a number in HTML scientific notation: `1.23 × 10<sup>4</sup>`.
fn html_scientific(value: f64, precision: usize) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    let exp = value.abs().log10().floor() as i32;
    let mantissa = value / 10f64.powi(exp);
    if exp == 0 {
        format!("{:.prec$}", mantissa, prec = precision)
    } else {
        format!(
            "{:.prec$}&thinsp;&times;&thinsp;10<sup>{}</sup>",
            mantissa,
            exp,
            prec = precision,
        )
    }
}

/// Write a Vec<f64> as a compact JSON array, rounding to 4 decimal places.
fn write_json_f64_array(w: &mut DebugWriter, data: &[f64]) {
    let _ = write!(w, "[");
    for (i, v) in data.iter().enumerate() {
        if i > 0 {
            let _ = write!(w, ",");
        }
        // Use compact formatting: 0 for zero, otherwise 4 decimal places.
        if *v == 0.0 {
            let _ = write!(w, "0");
        } else {
            let _ = write!(w, "{:.4}", v);
        }
    }
    let _ = write!(w, "]");
}
