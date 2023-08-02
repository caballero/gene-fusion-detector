#!/usr/bin/python3

# Script to parse a PAF aligment, it selects query sequences longer than the
# match region, in aims to select long-reads with gene fusions.

import sys
import csv
import gzip
import logging
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# PAF fields
fields = [
    "q_name",
    "q_len",
    "q_start",
    "q_end",
    "strand",
    "t_name",
    "t_len",
    "t_start",
    "t_end",
    "matches",
    "block_len",
    "map_qual",
    "tp",
    "cm",
    "s1",
    "s2",
    "dv",
    "rl"
]
    
def filter_aligns(in_file, out_file, min_len):
    tot_counter = 0
    sel_counter = 0
    sel_reads   = {}
    logging.debug("Parsing aligned PAF file")
    logging.debug(f"Minimal unmatched sequence: {min_len}")
    with open(out_file, "w") as out_fh:
        logging.debug(f"Output file: {out_file}")
        writer=csv.DictWriter(out_fh, 
                              delimiter  = "\t", 
                              quoting    = csv.QUOTE_NONE,
                              fieldnames = fields)
        with open(in_file, "r") as in_fh:
            logging.debug(f"Input file: {in_file}")        
            reader=csv.DictReader(in_fh, 
                                  delimiter  = "\t", 
                                  fieldnames = fields)
            for row in reader:
                tot_counter += 1
                match_diff = int(row["q_len"]) - int(row["t_len"])
                if match_diff >= min_len:
                    sel_counter += 1
                    name  = row["q_name"]
                    start = int(row["q_start"])
                    end   = int(row["q_end"])
                    if name in sel_reads:
                        if start < sel_reads[name]["start"]:
                            sel_reads[name]["start"] = start
                        if end > sel_reads[name]["end"]:
                            sel_reads[name]["end"] = end
                    else:
                        sel_reads[name] = {}
                        sel_reads[name]["start"] = start
                        sel_reads[name]["end"]   = end
                    writer.writerow(row)

    logging.debug(f"Read {tot_counter} matches")
    logging.debug(f"Selected {sel_counter} matches")
    
    if sel_counter > 0:
        return sel_reads
    else:
        logging.error("No filtered reads found, finishing")
        sys.exit(0)

def extract_seqs(fq_file, fa_file, reads, mask):
    logging.debug(f"Extracting sequences from {fq_file}")
    logging.debug(f"Final Fasta is {fa_file}")
    logging.debug(f"Masking mode: {mask}")
    count = 0
    with open(fa_file, "w") as out_fh:
        with gzip.open(fq_file, "rt") as in_fh:
            for record in SeqIO.parse(in_fh, "fastq"):
                if record.id in reads:
                    seq = list(record.seq)
                    id  = record.id
                    if mask:
                        start = reads[id]["start"]
                        end   = reads[id]["end"]
                        for i in range(start, end):
                            seq[i] = 'N'
                    new_record = SeqRecord(
                                    Seq("".join(seq)),
                                    id=id,
                                    description=''
                                )
                    SeqIO.write(new_record, out_fh, "fasta")
                    count += 1
    logging.debug(f"Wrote {count} records")


def main():
    # Minimal sequence size after match removal
    min_len = 50

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", 
                        help = "increase output verbosity",
                        action = "store_true")
    parser.add_argument("-i", "--input", 
                        help = "input file, PAF table",
                        type = str,
                        required = True)
    parser.add_argument("-o", "--output", 
                        help = "output file, filtered hits",
                        type = str,
                        required = True)
    parser.add_argument("-q", "--fastq", 
                        help = "input fastq file (gzip compressed)",
                        type = str,
                        required = True)
    parser.add_argument("-f", "--fasta", 
                        help = "output fasta file",
                        type = str,
                        required = True)
    parser.add_argument("-m", "--mask", 
                        help = "masks aligned regions in fasta output",
                        action = "store_true")
    parser.add_argument("-l", "--minlen", 
                        help = f"minimal length for unmatched region, default: {min_len}",
                        type = int,
                        default = min_len)
    
    
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    sel_reads = filter_aligns(in_file  = args.input,
                              out_file = args.output,
                              min_len  = args.minlen)
    
    extract_seqs(fq_file = args.fastq,
                 fa_file = args.fasta,
                 reads   = sel_reads,
                 mask    = args.mask)


if __name__ == "__main__":
    main()
