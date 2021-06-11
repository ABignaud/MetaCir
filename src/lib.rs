use std::error::Error;
use std::process;

mod sg;

pub struct Args {
    name: String,
    bam_list: String,
}

impl Args {
    // Add the arguments given by the user to the struct.
    pub fn new(args: &[String]) -> Result<Args, &str> {
        // Raise an error if not enough arguments
        if args.len() < 3 {
            return Err("Not enough arguments.");
        }

        // Assign the values to the struct.
        // TODO: Avoid the use of clone.
        let name = args[1].clone();
        let bam_list = args[2].clone();

        Ok(Args { name, bam_list })
    }
}

pub fn run(args: Args) -> Result<(), Box<dyn Error>> {
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

