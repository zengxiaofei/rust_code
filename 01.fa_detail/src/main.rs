#![allow(unused)]

use clap::Parser;
use std::io::{self, BufRead};
use std::fs::File;

#[derive(Parser)]

/// A rust implementation of fa_detail2.py (for 64-bit system)
struct Cli {
    /// Input fasta file
    fasta: std::path::PathBuf,
}

struct Bases {
    a: usize,
    t: usize,
    c: usize,
    g: usize,
    n: usize,
}

impl Bases {

    fn new() -> Self {
        Self {
            a: 0,
            t: 0,
            c: 0,
            g: 0,
            n: 0,
        }
    }

}


fn stat_seq(seq_id: &String, seq: &String, bases: &mut Bases) -> (usize, usize, Vec<usize>) {

    let seq_upper = seq.to_uppercase();
    let a_num = seq_upper.chars().filter(|&c| c == 'A').count();
    let t_num = seq_upper.chars().filter(|&c| c == 'T').count();
    let c_num = seq_upper.chars().filter(|&c| c == 'C').count();
    let g_num = seq_upper.chars().filter(|&c| c == 'G').count();
    let n_num = seq_upper.chars().filter(|&c| c == 'N').count();
    
    let scaffold_len: usize = seq_upper.chars().count();
    let gapless_scaffold_len = scaffold_len - n_num;

    if a_num + t_num + c_num + g_num + n_num != scaffold_len {
        println!("Seq {seq_id} has base(s) not in A, T, C, G, N");
    }

    bases.a += a_num;
    bases.t += t_num;
    bases.c += c_num;
    bases.g += g_num;
    bases.n += n_num;

    let ctgs: Vec<&str> = seq_upper.split('N').filter(|ctg| ctg.chars().count() > 0).collect();
    let ctg_lens: Vec<usize> = ctgs.iter().map(|ctg| ctg.chars().count()).collect();

    (scaffold_len, gapless_scaffold_len, ctg_lens)
}


fn parse_fasta(fasta_path: &std::path::PathBuf) -> (Vec<usize>, Vec<usize>, Vec<usize>, Bases) {

    let fasta = File::open(fasta_path).expect("Unable to read the file");
    let reader = io::BufReader::new(fasta);
    
    let mut seq_id  = String::new();
    let mut seq = String::new();
    
    let mut bases = Bases::new();

    // for returns of stat_seq function
    let mut scaffold_len = 0;
    let mut gapless_scaffold_len = 0;
    let mut ctg_lens: Vec<usize> = Vec::new();
    
    // for returns of parse_fasta function
    let mut scaffold_lens: Vec<usize> = Vec::new();
    let mut gapless_scaffold_lens: Vec<usize> = Vec::new();
    let mut all_ctg_lens: Vec<usize> = Vec::new();
    
    for line in reader.lines() {
        
        let line = line.expect("Unable to read the line");
        
        if line.is_empty() {
            continue;
        }

        let first_char = line.as_bytes()[0];
        
        if first_char == b'>' {
            seq_id = line[1..].to_string();
            if !seq.is_empty() {
                (scaffold_len, gapless_scaffold_len, ctg_lens) = stat_seq(&seq_id, &seq, &mut bases);
                scaffold_lens.push(scaffold_len);
                gapless_scaffold_lens.push(gapless_scaffold_len);
                all_ctg_lens.extend(ctg_lens);
            }
            seq.clear();
        } else {
            seq.push_str(&line);
        }
    }
    (scaffold_len, gapless_scaffold_len, ctg_lens) = stat_seq(&seq_id, &seq, &mut bases);
    scaffold_lens.push(scaffold_len);
    gapless_scaffold_lens.push(gapless_scaffold_len);
    all_ctg_lens.extend(ctg_lens);

    (scaffold_lens, gapless_scaffold_lens, all_ctg_lens, bases)
}


fn calculate_nx(mut lengths: Vec<usize>, tag: &str) {

    lengths.sort_by(|a, b| b.cmp(a));
    
    let total_len: usize = lengths.iter().sum();
    let num = lengths.len();
    let mut len_sum = 0;
    let (mut n10, mut n20, mut n30, mut n40, mut n50, mut n60, mut n70, mut n80, mut n90) = (0, 0, 0, 0, 0, 0, 0, 0, 0);

    println!("#### {tag} statistics ####");
    println!("Nx\tNumber\tLength");

    for (i, &len) in lengths.iter().enumerate() {
        len_sum += len;
        if len_sum as f64 / total_len as f64 >= 0.1 && n10 == 0 {
            n10 = len;
            println!("N10\t{}\t{n10}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.2 && n20 == 0 {
            n20 = len;
            println!("N20\t{}\t{n20}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.3 && n30 == 0 {
            n30 = len;
            println!("N30\t{}\t{n30}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.4 && n40 == 0 {
            n40 = len;
            println!("N40\t{}\t{n40}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.5 && n50 == 0 {
            n50 = len;
            println!("N50\t{}\t{n50}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.6 && n60 == 0 {
            n60 = len;
            println!("N60\t{}\t{n60}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.7 && n70 == 0 {
            n70 = len;
            println!("N70\t{}\t{n70}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.8 && n80 == 0 {
            n80 = len;
            println!("N80\t{}\t{n80}", i+1);
        }
        if len_sum as f64 / total_len as f64 >= 0.9 && n90 == 0 {
            n90 = len;
            println!("N90\t{}\t{n90}", i+1);
        }
    }

    println!("longest {0}: {1}, shortest {0}: {2}", tag, lengths[0], lengths[num-1]);
    println!("total number: {}, total length: {}\n", num, total_len)
}


fn main() {
    // parse arguments
    let args = Cli::parse();

    // parse FASTA file
    let mut scaffold_lens: Vec<usize> = Vec::new();
    let mut gapless_scaffold_lens: Vec<usize> = Vec::new();
    let mut all_ctg_lens: Vec<usize> = Vec::new();
    let mut bases = Bases::new();
    (scaffold_lens, gapless_scaffold_lens, all_ctg_lens, bases) = parse_fasta(&args.fasta);

    // output base statistics
    println!("#### base statistics ####");
    println!("A: {}, T: {}, C: {}, G: {}, N: {}", bases.a, bases.t, bases.c, bases.g, bases.n);
    let gc_ratio = (bases.c + bases.g) as f64 / (bases.a + bases.t + bases.c + bases.g) as f64;
    println!("GC%: {}\n", gc_ratio);
    
    // calculate Nx values
    calculate_nx(all_ctg_lens, "contig");
    calculate_nx(scaffold_lens, "scaffold");
    calculate_nx(gapless_scaffold_lens, "gapless scaffold");
}
