#[macro_use]
extern crate clap;

use clap::App;
use std::process;

fn main() -> Result<(), &'static str> {
    // Parse args with clap
    let yaml = load_yaml!("../cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    // Check if fasta assembly is given
    if !matches.is_present("fasta_file") {
        return Err("You need to give the fasta file.");
    }

    // Parse global args.
    let fasta_file = matches.value_of("fasta_file").unwrap();
    let min_size: usize = matches
        .value_of("min_size")
        .unwrap()
        .parse()
        .expect("min_size should be an integer.");
    let threads: usize = matches
        .value_of("threads")
        .unwrap()
        .parse()
        .expect("threads should be an integer.");
    let out_file: Option<String> = match matches.value_of("out_file") {
        Some(value) => Some(value.to_string()),
        None => None,
    };

    // Parse subcommands
    let result = match matches.subcommand_name() {
        Some("tr") => {
            let ref matches = matches.subcommand_matches("tr").unwrap();
            circularity::run_tr(
                &matches, fasta_file, min_size, threads, out_file,
            )
        }
        Some("sg") => {
            let ref matches = matches.subcommand_matches("sg").unwrap();
            circularity::run_sg(
                &matches, fasta_file, min_size, threads, out_file,
            )
        }
        Some("hic") => {
            let ref matches = matches.subcommand_matches("hic").unwrap();
            circularity::run_hic(
                &matches, fasta_file, min_size, threads, out_file,
            )
        }
        Some("all") => {
            let ref matches = matches.subcommand_matches("all").unwrap();
            circularity::run_all(
                &matches, fasta_file, min_size, threads, out_file,
            )
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
