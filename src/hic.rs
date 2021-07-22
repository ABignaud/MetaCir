use std::error::Error;
use std::process::Command;

const FILE: &'static str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/", "scripts/hic.py");

pub fn main(
    pairs_files: Vec<&str>,
    contig_data_file: &str,
    cov_threshold: &str,
    enzyme: &str,
    fasta_file: &str,
    min_size: usize,
    out_file: Option<String>,
    plot_dir: &str,
    threads: usize,
    tmp_dir: &str,
) -> Result<(), Box<dyn Error>> {
    // Defined path of the python script.
    // println!{file!()}
    // let path = Path::new(file!())
    //     .parent()
    //     .unwrap()
    //     .canonicalize()
    //     .expect("issue with the python script path.")
    //     .join("hic.py");
    // println!("{:?}", path);

    // Define default option for out_file
    let out_file = out_file.unwrap_or("./circular_hic.tsv".to_string());
    let out_file: &str = &out_file[..];
    let threads: &str = &format!["{}", threads];
    let min_size: &str = &format!["{}", min_size];
    let pairs_files: &str = &pairs_files.join(",");

    // Launch the main python function.
    let _cmd = Command::new(FILE)
        .args(&[
            "--enzyme",
            enzyme,
            "--cov-threshold",
            cov_threshold,
            "--min-size",
            min_size,
            "--plot",
            plot_dir,
            "--out-file",
            out_file,
            "--tmp-dir",
            tmp_dir,
            "--threads",
            threads,
            fasta_file,
            contig_data_file,
            pairs_files,
        ])
        .spawn()
        .expect("failed to execute hic detection.");
    Ok(())
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     /// Test
//     #[test]
//     fn test_test() {
//         test();
//     }
// }
