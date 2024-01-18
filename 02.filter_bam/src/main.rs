use clap::Parser;
use regex::Regex;
use std::process::{Command, Stdio};
use std::io::{BufRead, BufReader};

#[derive(Debug)]
#[derive(Parser)]

/// A rust implementation of filter_bam.py
struct Cli {
    /// Input BAM file, should not be sorted or name-sorted, do NOT sort it by coordinate
    bam: String,
    /// MAPQ cutoff, read pairs with both MAPQ >= this value will be kept
    mapq: u8,
    /// When this parameter is added, either end of a read pair having a MAPQ >= `mapq` will be kept
    #[arg(long)]
    single_end_mapq_filtering: bool,
    /// Edit distance cutoff, read pairs with single-end NM >= this value will be removed
    #[arg(long)]
    nm: Option<u32>,
    /// Remove PCR duplicates. Note that the duplicates should have been marked in the BAM file (flag 1024)
    #[arg(long)]
    remove_dup: bool,
    /// Remove singletons in bam file
    #[arg(long)]
    remove_singletons: bool,
    /// Threads for reading bam file
    #[arg(long, default_value_t = 8)]
    threads: usize,
}


fn parse_bam(args: &Cli) {

    let threads = (args.threads).to_string();
    let mut parameters = vec!["view", &args.bam, "-h", "-@", &threads];
    // skip read pair when DUP (1024 -> '0b10000000000')
    if args.remove_dup {
        parameters.push("-F 1024")
    }
    
    let command = Command::new("samtools")
        .args(&parameters)
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to execute samtools view");

    let stdout = command.stdout.expect("Failed to open stdout");
    let reader = BufReader::new(stdout);
    
    let re = Regex::new(r"NM:i:(\d+)").unwrap();
    let mut lines = reader.lines();
    let mut next_line = lines.next();
    while let Some(line1) = next_line {
        let line1_str = line1.as_ref().expect("Failed to read the line");
        let first_char = line1_str.as_bytes()[0];
        if first_char == b'@' {
            println!("{line1_str}");
        } else {
            if let Some(line2) = lines.next() {
                let line2_str = line2.as_ref().expect("Failed to read the line");
                let cols1: Vec<&str> = line1_str.split('\t').collect();
                let cols2: Vec<&str> = line2_str.split('\t').collect();
                // singletons are found, panic or remove them
                if cols1[0] != cols2[0] {
                    if args.remove_singletons {
                        next_line = Some(line2);
                        continue;
                    } else {
                        panic!("BAM may be coord-sorted or has singletons. Sort it by read name or try --remove-singletons")    
                    }
                }
                // filter NM
                if let Some(nm) = args.nm {

                    if let Some(match1) = re.captures(line1_str) {
                        if let Some(match2) = re.captures(line2_str) {
                            let nm1: u32 = match1[1].parse().unwrap();
                            let nm2: u32 = match2[1].parse().unwrap();
                            if nm1 >= nm || nm2 >= nm {
                                next_line = lines.next();
                                continue;
                            }
                        }
                    
                    }
                }
                // remove read pair if MAPQ < cutoff and output filtered records
                let mapq1: u8 = cols1[4].parse().unwrap();
                let mapq2: u8 = cols2[4].parse().unwrap();
                if args.single_end_mapq_filtering {
                    if mapq1 >= args.mapq || mapq2 >= args.mapq {
                        println!("{line1_str}\n{line2_str}");
                    }
                } else {
                    if mapq1 >= args.mapq && mapq2 >= args.mapq {
                        println!("{line1_str}\n{line2_str}");
                    }
                }
            } else {
                break;
            }
        }
        next_line = lines.next();
    }
}


fn main() {
    // parse arguments
    let args = Cli::parse();

    parse_bam(&args);
}
