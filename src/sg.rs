use itertools::izip;
use std::cmp;
use std::convert::TryFrom;
use std::error::Error;
use std::io;

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
pub fn extract_pairs(name: &String, bam_files_list: &String) -> Align {
    // Create an empty vector and max read size of 0.
    let mut read_size: u32 = 0;
    let mut ref_len: u32 = 0;
    let mut reads: Vec<Read> = Vec::new();

    // Extract bam files
    let bam_files_vec = extract_bam_files(bam_files_list);

    // Iterates on the list of bam files:
    for bam_file in bam_files_vec {
        // Build bam reader
        let mut reader = bam::IndexedReader::from_path(bam_file).unwrap();
        // Extract contig id and contig_length.
        let ref_id = reader.header().reference_id(name).unwrap();
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
                        position: u32::try_from(record.calculate_end())
                            .unwrap()
                            - 1,
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
    let align = Align {
        ref_len,
        read_size,
        reads,
    };
    align
}

/// Structure to handle ratio data.
pub struct PlusRatio {
    position: Vec<u32>,
    count_forward: Vec<u32>,
    count_total: Vec<u32>,
    ratio: Vec<Option<f32>>,
}

/// Function to build the ratio from the read mapping score
pub fn build_ratio(align: Align) -> PlusRatio {
    // Choose a trim length of the contig of half the read size.
    let trim_len: u32 = align.read_size / 2;
    // Choose a bin_length of the read size and define the number of bins
    let binned_len: u32 = align.ref_len / trim_len + 1;
    // Build of final vectors os size of length.
    let mut position: Vec<u32> = vec![0; binned_len as usize];
    let mut count_forward: Vec<u32> = vec![0; binned_len as usize];
    let mut count_total: Vec<u32> = vec![0; binned_len as usize];
    let mut ratio: Vec<Option<f32>> = vec![Some(0.0); binned_len as usize];

    for read in align.reads {
        if read.forward {
            count_forward[(read.position / trim_len) as usize] += 1;
        }
        count_total[(read.position / trim_len) as usize] += 1;
    }
    for (i, (forward, total)) in
        count_forward.iter().zip(&count_total).enumerate()
    {
        position[i] = i as u32 * trim_len;
        if *total != 0 {
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

pub fn write_ratio(ratios: PlusRatio) -> Result<(), Box<dyn Error>> {
    // start writer
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());

    // Iterates on ratio data:
    for (pos, fw, tot, rat) in izip!(
        &ratios.position,
        &ratios.count_forward,
        &ratios.count_total,
        &ratios.ratio
    ) {
        match rat {
            Some(rat) => wtr.write_record(&[
                pos.to_string(),
                fw.to_string(),
                tot.to_string(),
                rat.to_string(),
            ])?,
            None => wtr.write_record(&[
                pos.to_string(),
                fw.to_string(),
                tot.to_string(),
                "NA".to_string(),
            ])?,
        }
    }
    // flush writer
    wtr.flush()?;
    Ok(())
}

// Function to parse a comma separated bam files in a vector.
fn extract_bam_files(bam_list: &String) -> Vec<&str> {
    let bam_files_vec: Vec<&str> = bam_list.split(",").collect();
    bam_files_vec
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bam_files_extraction() {
        let content: String = "file1.bam,file-2.bam,file'3'.bam".to_string();
        let result: Vec<&str> = vec!["file1.bam", "file-2.bam", "file'3'.bam"];
        assert_eq!(extract_bam_files(&content), result);
    }
}
