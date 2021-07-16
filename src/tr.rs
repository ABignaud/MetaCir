//! Module to detect Terminal repeats from the fasta sequences of the contigs.
//! Two different situation are possible:
//!     DTR : Direct Terminal Repeats (Both repeats are in the sens)
//!     ITR : Inverted Terminal Repeats (One repeat is the reverse complement of
//!     the other.)

use bio::io::fasta;
use fnv::FnvHashMap;
use lazy_static::lazy_static;
use std::borrow::Borrow;
use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use std::str::from_utf8;

// Here are some functon to build the reverse complement of a sequence. A
// similar implementation is available in the bio rust module.

// Define translation dictionnary for reverse complement (IUPAC alphabet
// supported)
lazy_static! {
    static ref COMPLEMENT: Vec<u8> = {
        let mut comp: Vec<u8> = (0..=255).collect();
        for (&dna, &rev_dna) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
            comp[dna as usize] = rev_dna;
            comp[dna as usize + 32] = rev_dna + 32;  // lowercase variants
        }
        comp
    };
}

/// Return complement of a given DNA base (IUPAC alphabet supported).
fn complement(dna: u8) -> u8 {
    COMPLEMENT[dna as usize]
}

/// Return reverse complement of given sequence (IUPAC alphabet supported).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.into_iter()
        .rev()
        .map(|dna| complement(*dna.borrow()))
        .collect()
}

// Here are some useful struc for the terminal repeat description.

/// Define structure for Terminal repeat data
#[derive(Debug)]
pub struct TerminalRepeat<'a> {
    trtype: TrType,
    seq: &'a str,
    size: usize,
    count: usize,
    n_freq: f32,
    mode_freq: f32,
    flag: Option<String>,
}

/// Implement PartialEq of TerminalRepeat. The equality is just based on the
/// type, the sequence and the count as the other values directly depends on the
/// sequence.
impl PartialEq for TerminalRepeat<'_> {
    fn eq(&self, other: &Self) -> bool {
        (self.seq == other.seq)
            & (self.count == other.count)
            & (self.trtype == other.trtype)
    }
}

/// Define TR type with both possible values Direct Terminal Repeat (DTR) and
/// Inverted Terminal Repeat (ITR).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TrType {
    Dtr,
    Itr,
}

/// Compute the ratio of the most present character in a text file. For us it's
/// useful to find the mod ratio of the repeated sequence.
pub fn max_char_counts(text: &str) -> f32 {
    // Create a FnvHashMap.
    let mut counts: FnvHashMap<char, f32> = FnvHashMap::default();
    let mut max_count: f32 = 0.;
    // Count the numbers of occurences of each character.
    for character in text.chars() {
        *counts.entry(character).or_insert(0.) += 1.;
        // Update max count if necessary.
        max_count = {
            if counts[&character] > max_count {
                counts[&character]
            } else {
                max_count
            }
        }
    }
    // Divide by size to have the ratio.
    max_count / text.len() as f32
}

/// Compute the metadata of a terminal repeat sequence
pub fn compute_tr_metadata<'a>(
    seq: &'a str,
    tr_seq: &'a str,
    tr_type: TrType,
) -> TerminalRepeat<'a> {
    let size = tr_seq.len();
    let count = match tr_type {
        TrType::Dtr => seq.matches(tr_seq).count(),
        TrType::Itr => {
            seq.matches(tr_seq).count()
                + seq
                    .matches(
                        from_utf8(&reverse_complement(tr_seq.as_bytes()))
                            .unwrap(),
                    )
                    .count()
        }
    };
    let n_freq = tr_seq.matches("N").count() as f32 / size as f32;
    let mode_freq = max_char_counts(tr_seq);
    TerminalRepeat {
        trtype: tr_type,
        seq: tr_seq,
        size: size,
        count: count,
        n_freq: n_freq,
        mode_freq: mode_freq,
        flag: None,
    }
}

/// Extract Terminal repeat sequence
pub fn fetch_tr<'a>(
    seq: &'a str,
    seed_size: usize,
) -> Option<TerminalRepeat<'a>> {
    // Inititialize the values.
    let mut sub_seq = &seq[0..seed_size];
    let seq_size = seq.len();

    // Check if the sequence is repeated as DTR.
    let pos = seq.rfind(sub_seq);
    let dtr: Option<TerminalRepeat> = match pos {
        Some(pos) => {
            if pos > seq_size / 2 {
                sub_seq = &seq[pos..seq_size];
                if seq.find(sub_seq).unwrap() == 0 {
                    // Compute metadata of the Terminal Repeat.
                    let tr_type = TrType::Dtr;
                    Some(compute_tr_metadata(seq, sub_seq, tr_type))
                // Case that the end of the sequences does not map at beggining
                // of the sequence.
                } else {
                    None
                }
            // Case of a repeat in the first half of the sequence.
            } else {
                None
            }
        }
        // This case is impossible as the sequence should find at least itself.
        None => None,
    };

    // Check if there is an ITR.
    // Reverse complement the sequence
    let rev_seq: String = from_utf8(&reverse_complement(seq.as_bytes()))
        .unwrap()
        .to_string();
    let itr: Option<TerminalRepeat> = {
        if sub_seq == &rev_seq[0..seed_size] {
            let mut size: usize = seed_size + 1;
            // Elongate the sequence until it's not the same.
            while (&seq[size..size + 1] == &rev_seq[size..size + 1])
                & (size <= 200)
            {
                size += 1
            }
            sub_seq = &seq[0..size];
            let tr_type = TrType::Itr;
            // Compute metadata
            Some(compute_tr_metadata(seq, sub_seq, tr_type))
        } else {
            None
        }
    };

    // Check which is the longest terminal repeat to keep at the end.
    let dtr_size: usize = match &dtr {
        Some(tr) => tr.size,
        None => 0,
    };
    let itr_size: usize = match &itr {
        Some(tr) => tr.size,
        None => 0,
    };
    if dtr_size >= itr_size {
        dtr
    } else {
        itr
    }
}

/// Main function to search for DTR/ITR on a whole fasta assembly.
pub fn main(
    fasta_file: &str,
    min_size: usize,
    seed_size: usize,
    out_file: Option<String>,
) -> Result<(), Box<dyn Error>> {
    // Start writer depending on the parameters to the stdout or the output
    // file.
    let out: Box<dyn Write> = match out_file {
        Some(out_file) => Box::new(File::create(out_file)?),
        None => Box::new(io::stdout()),
    };
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(out);

    // Write Header
    wtr.write_record(&[
        "Contig",
        "Terminal Repeat",
        "Size",
        "Count",
        "Frequence of N",
        "Mode frequence",
        "Sequence",
    ])?;

    // start fasta reader
    let reader = fasta::Reader::new(File::open(fasta_file)?);

    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");
        let seq: &str = from_utf8(record.seq()).unwrap();
        if seq.len() < min_size {
            wtr.write_record(&[
                record.id(),
                "None",
                "0",
                "0",
                "0.",
                "0.",
                "None",
            ])?;
        } else {
            let tr: Option<TerminalRepeat> = fetch_tr(seq, seed_size);
            match tr {
                Some(tr) => wtr.write_record(&[
                    record.id(),
                    &format!("{:?}", tr.trtype),
                    &format!("{}", tr.size),
                    &format!("{}", tr.count),
                    &format!("{}", tr.n_freq),
                    &format!("{}", tr.mode_freq),
                    tr.seq,
                ])?,
                None => wtr.write_record(&[
                    record.id(),
                    "None",
                    "0",
                    "0",
                    "0.",
                    "0.",
                    "None",
                ])?,
            };
        }
        // flush writer
        // wtr.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    /// Test reverse complement function
    #[test]
    fn test_reverse_complement() {
        let seq = b"NBaTCG";
        let rev_seq: Vec<u8> = reverse_complement(seq);
        assert_eq!(rev_seq, b"CGAtVN");
    }

    /// Test max_char_counts function
    #[test]
    fn test_max_char_counts() {
        let seq = "ATCGATACGTAGCTAGCATC";
        let result = max_char_counts(seq);
        assert_eq!(result, 0.3);
    }

    /// Test fetch_tr function
    #[test]
    fn test_fecth_tr() {
        fn test_fecth_tr_code(
            seq: &str,
            expected_result: TerminalRepeat,
        ) -> () {
            let seed_size: usize = 10;
            let result = fetch_tr(seq, seed_size).unwrap();
            assert_eq!(result, expected_result);
            assert_eq!(result.size, expected_result.size);
            assert_eq!(result.n_freq, expected_result.n_freq);
            assert!(
                (result.mode_freq - expected_result.mode_freq).abs() < 0.01
            );
            ()
        }

        // DTR test
        let seq = "ATCGATACGTAGCTAGCATCGACAGTCGATATCGATACGTAGCTAGCATC";
        let expected_result = TerminalRepeat {
            trtype: TrType::Dtr,
            seq: "ATCGATACGTAGCTAGCATC",
            size: 20,
            count: 2,
            n_freq: 0.,
            mode_freq: 0.30,
            flag: None,
        };
        test_fecth_tr_code(seq, expected_result);

        // ITR test
        let seq = "ATCGATACGTAGCTAGCATCGACAGTCGATATTGCTAGCTACGTATCGAT";
        let expected_result = TerminalRepeat {
            trtype: TrType::Itr,
            seq: "ATCGATACGTAGCTAGCA",
            size: 18,
            count: 2,
            n_freq: 0.,
            mode_freq: 0.33,
            flag: None,
        };
        test_fecth_tr_code(seq, expected_result);
    }
}
