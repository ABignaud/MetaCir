#[macro_use]
extern crate clap;

use clap::App;
use std::process;

fn main() -> Result<(), &'static str> {
    // Parse args with clap
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    // Check if fasta assembly is given
    if !matches.is_present("fasta_file") {
        return Err("You need to give the fasta file.");
    }

    // Parse subcommands
    let result = match matches.subcommand_name() {
        Some("tr") => {
            let ref matches = matches.subcommand_matches("tr").unwrap();
            circularity::run_tr(&matches)
        }
        Some("sg") => {
            let ref matches = matches.subcommand_matches("sg").unwrap();
            circularity::run_sg(matches)
        }
        Some("hic") => {
            let ref matches = matches.subcommand_matches("hic").unwrap();
            circularity::run_hic(matches)
        }
        Some("all") => {
            let ref matches = matches.subcommand_matches("all").unwrap();
            circularity::run_all(matches)
        }
        None => {
            eprintln!("A subcommand is required. Do `circosight --help` to have more information.");
            process::exit(1);
        }
        _ => unreachable!(),
    };

    // Raise an error if one transmitted.
    if let Err(e) = result {
        eprintln! {"Circosight error: {}", e};
        process::exit(1);
    }
    Ok(())
}
