use clap::ArgMatches;
use std::error::Error;
use std::process;
// use log::{trace, info, error};

mod hic;
mod sg;
mod tr;

/// Parse the argument for terminal repeat detection main function.
pub fn run_tr<'a>(
    args: &ArgMatches<'a>,
    fasta_file: &str,
    min_size: usize,
    threads: usize,
    out_file: Option<String>,
) -> Result<(), Box<dyn Error>> {
    // Parse tr arguments.
    let seed_size: usize = args
        .value_of("seed_size")
        .unwrap()
        .parse()
        .expect("seed_size should be an integer.");

    // Raise a warning if more than one threads given as the other one are
    // useless.
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

/// Parse the arguments for circular contigs detection based on the shotgun
/// reads main function.
pub fn run_sg<'a>(
    args: &ArgMatches<'a>,
    fasta_file: &str,
    min_size: usize,
    threads: usize,
    out_file: Option<String>,
) -> Result<(), Box<dyn Error>> {
    // Parse sg arguments.
    let bam_files: Vec<&str> = args.values_of("bam_files").unwrap().collect();

    // Run the main tr function.
    if let Err(e) = sg::main(bam_files, fasta_file, min_size, out_file, threads)
    {
        eprintln!("Writing error: {}", e);
        process::exit(1)
    }
    Ok(())
}

/// Parse the arguments for circular contigs detection based on the hic reads
/// main function.
pub fn run_hic<'a>(
    args: &ArgMatches<'a>,
    fasta_file: &str,
    min_size: usize,
    threads: usize,
    out_file: Option<String>,
) -> Result<(), Box<dyn Error>> {
    // Parse hic arguments.
    let contig_data_file = args.value_of("contig_data_file").unwrap();
    let enzyme = args.value_of("enzyme").unwrap();
    let pairs_files: Vec<&str> =
        args.values_of("pairs_files").unwrap().collect();
    let plot_dir: &str = match args.value_of("plot") {
        Some(value) => value,
        None => "None",
    };
    let tmp_dir = args.value_of("tmp_dir").unwrap();

    // Run the main hic function
    if let Err(e) = hic::main(
        pairs_files,
        contig_data_file,
        enzyme,
        fasta_file,
        min_size,
        out_file,
        plot_dir,
        threads,
        tmp_dir,
    ) {
        eprintln!("Writing error: {}", e);
        process::exit(1)
    }
    Ok(())
}

/// Parse the arguments for the all three previous pipeline.
pub fn run_all<'a>(
    _args: &ArgMatches<'a>,
    _fasta_file: &str,
    _min_size: usize,
    _threads: usize,
    _out_file: Option<String>,
) -> Result<(), Box<dyn Error>> {
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
