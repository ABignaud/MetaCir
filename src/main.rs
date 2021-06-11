use std::env;
use std::process;

use circularity::Args;

fn main() {
    let args: Vec<String> = env::args().collect();

    let args = Args::new(&args).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    if let Err(e) = circularity::run(args) {
        eprintln!("Application error: {}", e);
        process::exit(1);
    }
}