use std::env;
use std::error::Error;
use std::process;

mod sg;

pub struct Config {
    name: String,
    bam_list: String,
}

impl Config {
    // Add the arguments given by the user to the struct.
    pub fn new(mut args: env::Args) -> Result<Config, &'static str> {
        args.next();

        let name = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a contig name."),
        };

        let bam_list = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a bam file."),
        };

        Ok(Config { name, bam_list })
    }
}

pub fn run(args: Config) -> Result<(), Box<dyn Error>> {
    let align = sg::extract_pairs(&args.name, &args.bam_list);
    let ratios = sg::build_ratio(align);
    if let Err(e) = sg::write_ratio(ratios) {
        eprintln!("Writing error: {}", e);
        process::exit(1);
    }
    Ok(())
}

// TESTS:
#[cfg(test)]
#[ignore]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

