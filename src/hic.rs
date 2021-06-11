use std::error::Error;
use serde::Deserialize;

// use std::fs;
// use sprs::CsMat;

// Define a structure for a contig.
pub struct Contig {
    name: String,
    matrix_file: String,
    frag_file: String,
}

#[derive(Debug, Deserialize)]
pub struct Record {
    row: u64,
    col: u64,
    val: u64,
}




impl Contig {
    // Add the arguments given by the user to the struct.
    pub fn new(args: &[String]) -> Result<Contig, &str> {
        // Raise an error if not enough arguments
        if args.len() < 4 {
            return Err("Not enough arguments.");
        }

        // Assign the values to the struct.
        // TODO: Avoid the use of clone.
        let name = args[1].clone();
        let matrix_file = args[2].clone();
        let frag_file = args[3].clone();

        Ok(Contig { name, matrix_file, frag_file })
    }

    // Read matrix based on the file given.
    pub fn get_matrix(&self) -> Result<(), csv::Error> {
        let csv = "row\tcol\tval\n7\t1\t12\n8\t0\t13\n";

        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(csv.as_bytes());
        // println!("{:?}", reader);
        println!{"testA"};
        for result in reader.deserialize::<Record>() {
            // The iterator yields Result<StringRecord, Error>, so we check the
            // error here.
            println!{"test"};
            println!("{:?}", result?);
        }
        Ok(())
    }
    // rows = matrix_data.iloc[:,0] 
    // cols = matrix_data.iloc[:,1]
    // values = matrix_data.iloc[:,2]
    // N = max(max(rows), max(cols))
    // matrix = sparse.coo_matrix((values, (rows, cols)), shape=(N+1,N+1))
    // self.matrix = (sparse.triu(matrix) + sparse.triu(matrix.T, k=1)
}