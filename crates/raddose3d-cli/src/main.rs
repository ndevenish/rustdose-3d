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

    match raddose3d_parser::parse(&input) {
        Ok(config) => {
            println!("Parsed {} crystal(s), {} beam(s), {} wedge(s)",
                config.crystals.len(), config.beams.len(), config.wedges.len());
        }
        Err(e) => {
            eprintln!("Parse error: {}", e);
            std::process::exit(1);
        }
    }
}
