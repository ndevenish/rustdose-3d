use clap::Parser;

#[derive(Parser)]
#[command(name = "raddose3d", about = "Radiation dose modelling for crystallography")]
struct Cli {
    /// Input file path
    #[arg(short, long)]
    input: String,
}

fn main() {
    let cli = Cli::parse();

    let input = std::fs::read_to_string(&cli.input)
        .unwrap_or_else(|e| {
            eprintln!("Error reading input file '{}': {}", cli.input, e);
            std::process::exit(1);
        });

    let config = match raddose3d_parser::parse(&input) {
        Ok(config) => config,
        Err(e) => {
            eprintln!("Parse error: {}", e);
            std::process::exit(1);
        }
    };

    match raddose3d::experiment::Experiment::run_from_config(&config) {
        Ok(_experiment) => {}
        Err(e) => {
            eprintln!("Simulation error: {}", e);
            std::process::exit(1);
        }
    }
}
