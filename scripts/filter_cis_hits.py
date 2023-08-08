#!/usr/bin/python

import re
import sys
import csv
import gzip
import logging
import argparse

# PAF fields
paf_fields = [
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

gtf_fields = [
    "chrom",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes"
]

gtf_trans_match = re.compile('transcript_id "(.+?)"')

def load_annotations(gtf_file):
    ann = {}
    counter = 0
    logging.debug(f"reading gene annotations from {gtf_file}")
    with gzip.open(gtf_file, "rt") as fh:
        reader=csv.DictReader(filter(lambda row: row[0]!='#', fh), 
                              delimiter  = "\t", 
                              fieldnames = gtf_fields)
        for row in reader:
            if row["feature"] == "transcript":
                res = re.search(gtf_trans_match, row['attributes'])
                if res:
                    counter += 1
                    trans_id = res.group(1)
                    ann[trans_id] = {}
                    ann[trans_id]["chrom"]  = row['chrom']
                    ann[trans_id]["start"]  = row['start']
                    ann[trans_id]["end"]    = row['end']
                    ann[trans_id]["strand"] = row['strand']
                    
                else:
                    logging.error(f"cannot find transcipt_id in {row}")

    logging.debug(f"found {counter} transcripts")
    if counter >= 1:
        return ann
    else:
        logging.error("No annotations found, aborting")
        sys.exit(1)     


                
def main():
    # Minimal sequence size after match removal
    max_len = 5000000 

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", 
                        help = "increase output verbosity",
                        action = "store_true")
    parser.add_argument("-i", "--input1", 
                        help = "first input file, PAF table from map 1",
                        type = str,
                        required = True)
    parser.add_argument("-j", "--input2", 
                        help = "second input file, PAF table from map 2",
                        type = str,
                        required = True)
    parser.add_argument("-o", "--output", 
                        help = "output file, filtered hits",
                        type = str,
                        required = True)
    parser.add_argument("-g", "--gtf", 
                        help = "gene annotation file, GTF format (gzipped)",
                        type = str,
                        required = True)
    parser.add_argument("-m", "--maxlen", 
                        help = f"maximal length for matches distance, default: {max_len}",
                        type = int,
                        default = max_len)
    
    
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    annot = load_annotations(args.gtf)


if __name__ == "__main__":
    main()
