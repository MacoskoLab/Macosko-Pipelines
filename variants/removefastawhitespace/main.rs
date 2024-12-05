use std::path::Path;
use std::fs::File;
use std::io::{Read, BufRead, BufReader, Write};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        panic!("Usage: removefastawhitespace fasta.fa");
    }
    
    // define input/output paths
    assert!(Path::new("fasta_whitespace").is_dir(), "fasta_whitespace not found");
    assert!(Path::new("fasta_nowhitespace").is_dir(), "fasta_nowhitespace not found");
    let input_filename = ["fasta_whitespace", &args[1]].join("/");
    let output_filename = ["fasta_nowhitespace", &args[1]].join("/");
    
    // make sure input exists (and output doesn't exist)
    assert!(Path::new(&input_filename).exists(), "ERROR: input fasta file does not exist ({})", &args[1]);
    assert!(!Path::new(&output_filename).exists(), "ERROR: output fasta file already exists ({})", &args[1]);
    
    // make sure the input file is just text and newlines
    let reader = BufReader::new(File::open(&input_filename).expect("ERROR: could not load input fasta"));
    for (pos, byte_result) in reader.bytes().enumerate() {
        let byte = byte_result.expect("ERROR: could not read byte");
        match byte {
          0x0A => continue, // newline
          0x20..=0x7E => continue, // text
          _ => panic!("ERROR: disallowed byte 0x{:02X} at position {}", byte, pos),
        }
    }
    
    // initialize reader and writer
    let infile = File::open(&input_filename).expect("ERROR: could not load input fasta");
    let mut outfile = File::create(&output_filename).expect("ERROR: could not create output");
    let mut buffer = String::new(); // for fasta sequences
    
    eprintln!("Removing whitespace...");
    
    // loop
    for line_result in BufReader::new(infile).lines() {
        let line = line_result.expect("ERROR: could not read line");
        if line.starts_with('>') {
            if !buffer.is_empty() {
                writeln!(outfile, "{}", buffer).expect("ERROR: could not write line");
                buffer.clear();
            }
            writeln!(outfile, "{}", line.trim()).expect("ERROR: could not write line");
        } else {
            buffer.push_str(line.trim());
        }
    }
    writeln!(outfile, "{}", buffer).expect("ERROR: could not write line");
    
    // make sure the output file is just text and newlines
    let reader = BufReader::new(File::open(&output_filename).expect("ERROR: could not load output fasta"));
    for (pos, byte_result) in reader.bytes().enumerate() {
        let byte = byte_result.expect("ERROR: could not read byte");
        match byte {
          0x0A => continue, // newline
          0x20..=0x7E => continue, // text
          _ => panic!("ERROR: disallowed byte 0x{:02X} at position {}", byte, pos),
        }
    }
    
    // do a check
    let mut numheader = 0; let mut numlines = 0;
    let reader = BufReader::new(File::open(&output_filename).expect("ERROR: could not load output fasta"));
    for byte_result in reader.bytes() {
        let byte = byte_result.expect("ERROR: could not read byte");
        match byte {
          b'>' => numheader+=1,
          b'\n' => numlines+=1,
          _ => continue,
        }
    }
    assert!(numheader*2 == numlines, "ERROR: number of lines is not number of headers x 2");
    
    eprintln!("Done! Run this command to create the fasta index file:\nsamtools faidx fasta_nowhitespace/{}", &args[1]);
}
