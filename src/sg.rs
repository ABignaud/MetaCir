use std::cmp;
use std::convert::TryFrom;
use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use bio::io::fasta;
use peroxide::pnorm;
use peroxide::fuga::*;

/// Structure for a read with only start position and the sense of the reads.
#[derive(Debug)]
pub struct Read {
    position: u32,
    forward: bool,
}

/// Structure with the alignement information of the contig and the alignment
/// position.
pub struct Align {
    ref_len: u32,
    read_size: u32,
    reads: Vec<Read>,
}

/// Function to extract the pairs mapping on one contig from the whole bam
/// alignment.
pub fn extract_pairs(name: &str, bam_files: &Vec<&str>) -> Option<Align> {
    // Create an empty vector and max read size of 0.
    let mut read_size: u32 = 0;
    let mut ref_len: u32 = 0;
    let mut reads: Vec<Read> = Vec::new();

    // Iterates on the list of bam files:
    for bam_file in bam_files {
        // Build bam reader
        let mut reader = bam::IndexedReader::from_path(bam_file).unwrap();
        // Extract contig id and contig_length.
        if let Some(ref_id) = reader.header().reference_id(name) {
            ref_len = reader.header().reference_len(ref_id).unwrap();
            // Iterate on the alignement region (the whole contig).
            for record in
                reader.fetch(&bam::Region::new(ref_id, 0, ref_len)).unwrap()
            {
                let record = record.unwrap();
                // Look at the read size to compute max read size
                read_size = cmp::max(read_size, record.query_len());
                // Detect strand of the read based on the flag.
                let flag = record.flag();
                if flag.mate_is_mapped() {
                    if flag.is_reverse_strand() && !flag.mate_is_reverse_strand() {
                        reads.push(Read {
                            position: u32::try_from(record.start()).unwrap(),
                            forward: false,
                        });
                    } else if !flag.is_reverse_strand()
                        && flag.mate_is_reverse_strand()
                    {
                        reads.push(Read {
                            position: u32::try_from(record.start()).unwrap(),
                            forward: true,
                        });
                    }
                }
            }
        }
    }
    if reads.len() > 0 {
        let align = Align {
            ref_len,
            read_size,
            reads,
        };
        Some(align)
    } else {
        None
    }
}

/// Structure to handle ratio data.
#[derive(Debug, Default)]
pub struct PlusRatio {
    position: Vec<u32>,
    count_forward: Vec<u32>,
    count_total: Vec<u32>,
    ratio: Vec<Option<f32>>,
}

impl PlusRatio {

    /// Compute the start end to keep as it could have no coverage at the
    /// beginning/end. Remove extremities if there have a total count with a
    /// distance at the mean greater than 1.5 time the standard deviation. The
    /// last value at the end is removed as we expect no value as we took left
    /// position.
    fn start_end(&self) -> (usize, usize) {
        let mut start: usize = 0;
        let mut end: usize =  self.count_total.len() - 2;
        let mean: f32 = self.count_total.iter().sum::<u32>() as f32 / (end + 2) as f32;
        let mut std: f32 = 0.;
        for i in self.count_total.iter() {
            std += (*i as f32 - mean).powf(2.);
        };
        let threshold = 1.5 * (std / (end + 2) as f32).sqrt();
        while (self.count_total[start] as f32 - mean).abs() > threshold {
            start += 1
        };
        while (self.count_total[end] as f32 - mean).abs() > threshold {
            end -= 1
        };
        if start >= 3 {
            start = 2
        };
        if end <= self.count_total.len() - 5{
            end = self.count_total.len()
        };
        (start, end)
    }

    /// Implement the mean of the ratio. The two extremities are removed.
    fn mean(&self, start: usize, end: usize) -> (f32, f32) {
        let mut sum: f32 = 0.;
        let mut n: f32 = 0.;
        for i in &self.ratio[start..end] {
            if let Some(i) = i {
                sum += i;
                n += 1.;
            };
        }
        if n == 0. {
            (sum, n)
        } else {
            (sum / n, n)
        }
    }

    /// Implement the standard deviation of the ratio as the mean.
    fn std(&self, start: usize, end: usize, mean: f32, n: f32) -> f32 {
        let mut sum_std: f32 = 0.;
        for i in &self.ratio[start..end] {
            if let Some(i) = i {
                sum_std += (i - mean).powf(2.);
            }
        }
        (sum_std / n).sqrt() 
    }

    fn total_contact(&self) -> u32 {
        self.count_total.iter().sum()
    }
}

/// Function to build the ratio from the read mapping score
pub fn build_ratio(align: Align) -> PlusRatio {
    // Choose a trim length of the contig of the read size.
    let trim_len: u32 = align.read_size;
    // Choose a bin_length of the read size and define the number of bins
    let binned_len: u32 = (align.ref_len / trim_len) + 2;
    // Compute the size of the residual binned.
    let res_bin: u32 = binned_len % trim_len;
    // Build of final vectors os size of length.
    let mut position: Vec<u32> = vec![0; binned_len as usize];
    let mut count_forward: Vec<u32> = vec![0; binned_len as usize];
    let mut count_total: Vec<u32> = vec![0; binned_len as usize];
    let mut ratio: Vec<Option<f32>> = vec![Some(0.0); binned_len as usize];

    for read in align.reads {
        // Play a trick on index to have the residual bin in the middle of the 
        // contig.
        let index: usize = {
            if read.position < (trim_len * binned_len / 2) {
                (read.position / trim_len) as usize 
            } else {
                ((read.position + res_bin) / trim_len) as usize
            }
        };
        count_total[index] += 1;
        // Count only for forward reads
        if read.forward {
            count_forward[index] += 1;
        }
    }
    // Compute ratio
    for (i, (forward, total)) in
        count_forward.iter().zip(&count_total).enumerate()
    {
        if (i as u32) < (binned_len as u32 / 2) { 
            position[i] = i as u32 * trim_len;
        } else {
            position[i] = i as u32 * trim_len + res_bin;
        }
        if *total > 30 {
            ratio[i] = Some(*forward as f32 / *total as f32);
        } else {
            ratio[i] = None
        }
    }
    let ratios = PlusRatio {
        position,
        count_forward,
        count_total,
        ratio,
    };
    ratios
}

/// Structure for saving results.
#[derive(Debug, PartialEq)]
pub struct Score {
    circular: String,
    contact: u32,
    score: Option<f64>,
    flag: String,
}

/// Compute the score. 
fn compute_score(ratio: f32, mean: f32, std: f32) -> f64 {
    let score = if ratio > mean {
        2. * pnorm!(2. * mean - ratio, mean , std)
    } else {
        2. * pnorm!(ratio, mean, std)
    };
    score
}

/// Fonction to compute the score. 
pub fn build_score(ratio: PlusRatio) -> Score {
    
    // Look for start and end position.
    let (start, end) = ratio.start_end();

    // compute mean and count of non zeros values.
    let (mean, n)  = ratio.mean(start + 1, end);

    // Look for the values at the extremities:
    let left = match ratio.ratio[start] {
        Some(value) => value,
        None => -1.,
    };
    let right = match ratio.ratio[end] {
        Some(value) => value,
        None => -1.,
    };

    // Case of not enough coverage on the contig or partial coverage.
    if n < ((end - start) as f32) / 2. {
        Score {
            circular: "ND".to_string(),
            contact: ratio.total_contact(),
            score: None,
            flag: "Not enough coverage on at least half the contig.".to_string(),
        }
    }

    // Case where there is no coverage at the extremities.
    else if (right == -1.) | (left == -1.) {
        Score {
            circular: "ND".to_string(),
            contact: ratio.total_contact(),
            score: None,
            flag: "No coverage in at least one extremity".to_string(),
        }
    }

    // Compute the score.
    else {
        // compute standard deviation.
        let std = ratio.std(start + 1, end, mean, n);
        let score = (compute_score(left, mean, std) + compute_score(right, mean, std)) / 2.;
        let circular = if score > 0.5 {
            "true".to_string()
        } else {
            "false".to_string()
        };
        // println!("{:?}", ratio);
        Score {
            circular: circular,
            contact: ratio.total_contact(),
            score: Some(score),
            flag: "Ok".to_string(),
        }
    }
}

/// Main function to detect circular contigs based on shotgun reads.
pub fn main(
    bam_files: Vec<&str>,
    fasta_file: &str,
    min_size: usize,
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
        "Circular",
        "Pairs_mapped",
        "Score",
        "Flag"
    ])?;

    // start fasta reader
    let reader = fasta::Reader::new(File::open(fasta_file)?);
    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");
        let seq = record.seq();
        if seq.len() < min_size {
            wtr.write_record(&[
                record.id(),
                "ND",
                "ND",
                "ND",
                "Too small.",
            ])?;
        } else {
            let align = extract_pairs(record.id(), &bam_files);
            match align {
                Some(align) => {
                    let ratio = build_ratio(align);
                    let score = build_score(ratio);
                    let score_value = match score.score {
                        Some(score) => format!("{:?}", score),
                        None => "ND".to_string(),
                    };
                    wtr.write_record(&[
                        record.id(),
                        &format!("{}", score.circular),
                        &format!("{}", score.contact),
                        &format!("{}", score_value),
                        &format!("{}", score.flag),
                    ])?;
                },
                None => {
                    wtr.write_record(&[
                        record.id(),
                        "ND",
                        "0",
                        "ND",
                        "No reads mapped on the contig.",
                    ])?;
                },
            };
        }
        // flush writer
        wtr.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    /// Test build score
    #[test]
    fn test_build_score() {
        let ratios = PlusRatio { 
            position: vec![0, 75, 150, 225, 300, 375, 450, 525, 600, 675],
            count_forward: vec![267, 142, 128, 121, 101, 68, 90, 96, 75, 50],
            count_total: vec![410, 324, 255, 262, 209, 218, 196, 199, 172, 123],
            ratio: vec![Some(0.6512195), Some(0.4382716), Some(0.5019608), Some(0.46183205), None, Some(0.3119266), None, Some(0.48241207), Some(0.4360465), None],
        };
        let expected_result = Score {
            circular: "true".to_string(),
            contact: 2368,
            score: Some(0.9747624989970725),
            flag: "Ok".to_string(),
        };
        let (start, end) = ratios.start_end();
        let (mean, n) = ratios.mean(start + 1, end);
        assert_eq!(start, 1);
        assert_eq!(end, 8);
        assert_eq!(mean, 0.43953288);
        assert_eq!(ratios.std(start + 1, end, mean, n), 0.07502747);
        assert_eq!(ratios.total_contact(), 2368);
        let result: Score = build_score(ratios);
        assert_eq!(result, expected_result);
    }
}
