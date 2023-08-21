#!/usr/bin/python3

# Script to reads from a FastQ file.

import gzip
import logging
import argparse
from Bio import SeqIO

def get_reads(in_file):
    logging.debug(f"reading read list from {in_file}")
    count = 0
    reads = ()
    with open(in_file, "r") as fh:
        for line in fh:
            id = line.rstrip()
            reads.append(id)
            count += 1
    logging.debug(f"found {count} reads ids")
    return reads

def extract_seqs(in_file, out_file, reads):
    logging.debug(f"Extracting sequences from {in_file}")
    logging.debug(f"Final Fastq is {out_file}")
    count = 0
    with gzip.open(out_file, "wt") as out_fh:
        with gzip.open(in_file, "rt") as in_fh:
            for record in SeqIO.parse(in_fh, "fastq"):
                if record.id in reads:
                    SeqIO.write(record, out_fh, "fastq")
                    count += 1
    logging.debug(f"Wrote {count} records")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", 
                        help = "increase output verbosity",
                        action = "store_true")
    parser.add_argument("-i", "--input", 
                        help = "input file, compressed Fastq",
                        type = str,
                        required = True)
    parser.add_argument("-o", "--output", 
                        help = "output file, compressed Fastq",
                        type = str,
                        required = True)
    parser.add_argument("-l", "--list", 
                        help = "read list file (one per line)",
                        type = str,
                        required = True)    
    
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    reads = get_reads(in_file  = args.list)
    
    extract_seqs(in_file  = args.input,
                 out_file = args.output,
                 reads    = reads)

    logging.debug("Completed")

if __name__ == "__main__":
    main()
