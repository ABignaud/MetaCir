use clap::ArgMatches;
use std::error::Error;
use std::process;
// use log::{trace, info, error};

mod tr;

pub fn run_tr<'a>(args: &ArgMatches<'a>) -> Result<(), Box<dyn Error>> {
    // Parse useful arguments.
    let fasta_file = args.value_of("fasta_file").unwrap();
    let min_size: usize = args
        .value_of("min_size")
        .unwrap()
        .parse()
        .expect("min_size should be an integer.");
    let seed_size: usize = args
        .value_of("seed_size")
        .unwrap()
        .parse()
        .expect("seed_size should be an integer.");
    let out_file: Option<String> = match args.value_of("out_file") {
        Some(value) => Some(value.to_string()),
        None => None,
    };

    // Raise a warning if more than one threads given as the other one are
    // useless.
    let threads: u16 = args
        .value_of("threads")
        .unwrap()
        .parse()
        .expect("threads should be an integer.");
    if threads != 1 {
        eprintln!("WARNING: TR search is not parallelized only 1 thread would be use.")
    };

    // Run the main tr function.
    if let Err(e) = tr::main(fasta_file, min_size, seed_size, out_file) {
        eprintln!("Writing error: {}", e);
        process::exit(1)
    }
    Ok(())
}

pub fn run_sg<'a>(args: &ArgMatches<'a>) -> Result<(), Box<dyn Error>> {
    if let Some(threads) = args.value_of("threads") {
        println!("Number of threads: {}", threads);
    }
    let bam_files: Vec<_> = args.values_of("bam_files").unwrap().collect();
    println!("Bam files: {:?}", bam_files);
    println!("Running Shotgun reads pipeline.");
    Ok(())
}

pub fn run_hic<'a>(_args: &ArgMatches<'a>) -> Result<(), Box<dyn Error>> {
    println!("Not implemented yet, ask the developper to do his job.");
    Ok(())
}

pub fn run_all<'a>(_args: &ArgMatches<'a>) -> Result<(), Box<dyn Error>> {
    println!("Not implemented yet, ask the developper to do his job.");
    Ok(())
}

// pub fn run(args: Config) -> Result<(), Box<dyn Error>> {
//     // SHOTGUN RUN
//     // let align = sg::extract_pairs(&args.name, &args.bam_list);
//     // let ratios = sg::build_ratio(align);
//     // if let Err(e) = sg::write_ratio(ratios) {
//     // eprintln!("Writing error: {}", e);
//     // process::exit(1);
//     // }
// }
