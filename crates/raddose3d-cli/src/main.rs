use std::io::{self, Read};

use clap::Parser;
use raddose3d::beam;
use raddose3d::crystal;
use raddose3d::experiment::Experiment;
use raddose3d::output::{
    OutputDWDs, OutputDoseStateHTML, OutputFinalDoseStateCSV, OutputFinalDoseStateR, OutputRDECSV,
    OutputSummaryCSV, OutputSummaryText, OutputVoxelDose, OutputVoxelFluences,
};
use raddose3d::wedge::Wedge;
use raddose3d::writer::{file_writer, stdout_writer, TeeWriter};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[command(
    name = "raddose3d",
    about = "Radiation dose modelling for crystallography",
    long_about = None,
    // Disable clap's built-in --version so we handle -V ourselves (matching Java).
    disable_version_flag = true,
    // Disable clap's built-in --help so we can add -? alias ourselves.
    disable_help_flag = true,
)]
struct Cli {
    /// Input file path ('-' reads from stdin)
    #[arg(short = 'i', long = "in", value_name = "FILE")]
    input: Option<String>,

    /// Prefix for output files
    #[arg(
        short = 'p',
        long = "prefix",
        value_name = "PREFIX",
        default_value = "output-"
    )]
    prefix: String,

    /// Path to RADDOSE executable (for CoefCalcRaddose)
    #[arg(short = 'r', long = "raddose", value_name = "PATH")]
    raddose_exe: Option<String>,

    /// Test run — print message and exit without running simulation
    #[arg(short = 't', long = "test")]
    test: bool,

    /// Show version information
    #[arg(short = 'V', long = "version")]
    version: bool,

    /// Show help
    #[arg(short = '?', long = "help", action = clap::ArgAction::Help)]
    help: Option<bool>,
}

fn main() {
    let cli = Cli::parse();

    if cli.version {
        println!("RADDOSE-3D version {}", VERSION);
        println!("Please cite:");
        println!("  Zeldin, Gerstel, Garman. (2013). J. Appl. Cryst. 46, 1225-1230.");
        println!("  http://dx.doi.org/10.1107/S0021889813011461");
        return;
    }

    if cli.test {
        println!("Test run. No actual calculations will take place.");
        return;
    }

    if let Some(path) = &cli.raddose_exe {
        println!("raddose executable set to {}", path);
    }

    let input_path = match cli.input {
        Some(p) => p,
        None => {
            eprintln!("No input file specified. Use -i <file> or -i - for stdin.");
            eprintln!("Run with -? for usage information.");
            std::process::exit(1);
        }
    };

    // Read input text from file or stdin.
    let input = if input_path == "-" {
        println!("Read from console");
        let mut buf = String::new();
        io::stdin().read_to_string(&mut buf).unwrap_or_else(|e| {
            eprintln!("Error reading from stdin: {}", e);
            std::process::exit(1);
        });
        buf
    } else {
        println!("Load File: {}", input_path);
        std::fs::read_to_string(&input_path).unwrap_or_else(|e| {
            eprintln!("Error reading input file '{}': {}", input_path, e);
            std::process::exit(1);
        })
    };

    let mut config = match raddose3d_parser::parse(&input) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Parse error: {}", e);
            std::process::exit(1);
        }
    };

    // Resolve relative file paths in the config relative to the input file's directory.
    if input_path != "-" {
        if let Some(base) = std::path::Path::new(&input_path).parent() {
            config.resolve_paths(base);
        }
    }

    let prefix = &cli.prefix;
    println!("No output specifications given. Using defaults.");

    let mut experiment = Experiment::new();
    add_default_observers(&mut experiment, prefix);

    use raddose3d::parser::config::ConfigItem;
    for item in &config.items {
        match item {
            ConfigItem::Crystal(c) => match crystal::create_crystal(c) {
                Ok(crystal) => experiment.set_crystal(crystal),
                Err(e) => {
                    eprintln!("Crystal error: {}", e);
                    std::process::exit(1);
                }
            },
            ConfigItem::Beam(b) => match beam::create_beam(b) {
                Ok(beam) => experiment.set_beam(beam),
                Err(e) => {
                    eprintln!("Beam error: {}", e);
                    std::process::exit(1);
                }
            },
            ConfigItem::Wedge(w) => {
                let wedge = Wedge::from_config(w);
                experiment.expose_wedge(&wedge);
            }
        }
    }

    experiment.close();
}

/// Attach the default set of output observers, mirroring Java's `setDefaultObservers()`.
///
/// Default outputs:
/// - `{prefix}Summary.csv`
/// - `{prefix}Summary.txt` + stdout (tee)
/// - `{prefix}DoseState.csv`
/// - `{prefix}DoseState.R`
/// - `{prefix}RDE.csv`
/// - `{prefix}DWDs.csv`
/// - `{prefix}VoxelDose.csv`
/// - `{prefix}VoxelFluences.csv`
/// - Progress indicator → stderr
fn add_default_observers(experiment: &mut Experiment, prefix: &str) {
    // SummaryCSV
    match file_writer(format!("{}Summary.csv", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputSummaryCSV::new(w))),
        Err(e) => eprintln!("Could not initialize OutputSummaryCSV: {}", e),
    }

    // SummaryText — tee to stdout and file
    match file_writer(format!("{}Summary.txt", prefix)) {
        Ok(f) => {
            let tee = Box::new(TeeWriter::new(vec![stdout_writer(), f]));
            experiment.add_observer(Box::new(OutputSummaryText::new(tee)));
        }
        Err(e) => {
            eprintln!("Could not open {}Summary.txt: {}", prefix, e);
            experiment.add_observer(Box::new(OutputSummaryText::new(stdout_writer())));
        }
    }

    // FinalDoseStateCSV
    match file_writer(format!("{}DoseState.csv", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputFinalDoseStateCSV::new(w))),
        Err(e) => eprintln!("Could not initialize OutputFinalDoseStateCSV: {}", e),
    }

    // FinalDoseStateR
    match file_writer(format!("{}DoseState.R", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputFinalDoseStateR::new(w))),
        Err(e) => eprintln!("Could not initialize OutputFinalDoseStateR: {}", e),
    }

    // RDECSV
    match file_writer(format!("{}RDE.csv", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputRDECSV::new(w))),
        Err(e) => eprintln!("Could not initialize OutputRDECSV: {}", e),
    }

    // DWDs
    match file_writer(format!("{}DWDs.csv", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputDWDs::new(w))),
        Err(e) => eprintln!("Could not initialize OutputDWDs: {}", e),
    }

    // VoxelDose
    match file_writer(format!("{}VoxelDose.csv", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputVoxelDose::new(w))),
        Err(e) => eprintln!("Could not initialize OutputVoxelDose: {}", e),
    }

    // VoxelFluences
    match file_writer(format!("{}VoxelFluences.csv", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputVoxelFluences::new(w))),
        Err(e) => eprintln!("Could not initialize OutputVoxelFluences: {}", e),
    }

    // HTML report with Plotly visualizations
    match file_writer(format!("{}Report.html", prefix)) {
        Ok(w) => experiment.add_observer(Box::new(OutputDoseStateHTML::new(w))),
        Err(e) => eprintln!("Could not initialize OutputDoseStateHTML: {}", e),
    }

    // Progress indicator is now printed inline by expose_rd3d().
}
