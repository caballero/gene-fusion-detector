#!/usr/bin/python

# Script to compare 2 PAF mappings to identify possible gene fusions in long-read data

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

gtf_match = re.compile('gene_id "(.+?)"; transcript_id "(.+?)"')

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
                res = re.search(gtf_match, row['attributes'])
                if res:
                    counter += 1
                    gene_id = res.group(1)
                    trans_id = res.group(2)
                    ann[trans_id] = {}
                    ann[trans_id]['gene_id'] = gene_id
                    ann[trans_id]["chrom"]   = row['chrom']
                    ann[trans_id]["start"]   = int(row['start'])
                    ann[trans_id]["end"]     = int(row['end'])
                    ann[trans_id]["strand"]  = row['strand']
                    
                else:
                    logging.error(f"cannot find transcipt_id in {row}")

    logging.debug(f"found {counter} transcripts")
    if counter >= 1:
        return ann
    else:
        logging.error("No annotations found, aborting")
        sys.exit(1)     


def load_matches(input_file, annot):
    logging.debug(f"loading mapping from {input_file}")
    matches = {}
    counter = 0
    with open(input_file, "r") as fh:
        reader=csv.DictReader(fh,
                              delimiter="\t", 
                              fieldnames=paf_fields)
        for row in reader:
            read  = row['q_name']
            name  = row['t_name'].split("|")
            trans = name[0]
            if trans in annot:
                matches[read] = {}
                matches[read]['trans']   = trans
                matches[read]['chrom']   = annot[trans]['chrom']
                matches[read]['start']   = annot[trans]['start']
                matches[read]['end']     = annot[trans]['end']
                matches[read]['strand']  = annot[trans]['strand']
                matches[read]['gene_id'] = annot[trans]['gene_id']
                counter += 1
            else:
                logging.error(f"No annotation for {trans}")
    logging.debug (f"loaded {counter} annotation for reads")
    if counter > 0:
        return matches
    else:
        logging.error(f"No annotations found from {input_file}")
        sys.exit(1)


def get_distance(s1, e1, s2, e2):
    if s1 < s2 and e1 <= s2:   # gen1 -> gen2
        return s2 - e1
    elif s2 < s1 and e2 <= s1: # gen2 -> gen1
        return s1 - e2 
    else:                      # overlapped regions
        return -1


def filter_matches(out_file, input_file, matches, annot, maxlen):
    logging.debug(f"comparing mapping from {input_file}, outfile is {out_file}")
    counter = 0
    out_header = ["#read_id",
                  "g1",
                  "t1",
                  "t1_chrom",
                  "t1_start",
                  "t1_end",
                  "t1_strand",
                  "dist",
                  "g2",
                  "t2",
                  "t2_chrom",
                  "t2_start",
                  "t2_end",
                  "t2_strand"]
    with open(out_file, "w") as out:
        writer = csv.writer(out,
                            delimiter="\t",
                            quoting=csv.QUOTE_NONE)
        writer.writerow(out_header)
        with open(input_file, "r") as fh:
            reader=csv.DictReader(fh,
                                delimiter="\t", 
                                fieldnames=paf_fields)
            for row in reader:
                read  = row['q_name']
                name  = row['t_name'].split("|")
                m2_trans = name[0]
                if m2_trans in annot:
                    if read in matches:
                        m1_trans  = matches[read]['trans'] 
                        m1_chr    = matches[read]['chrom']
                        m1_start  = matches[read]['start']
                        m1_end    = matches[read]['end']
                        m1_strand = matches[read]['strand']
                        m1_gene   = matches[read]['gene_id']
                        m2_chr    = annot[m2_trans]['chrom']
                        m2_start  = annot[m2_trans]['start']
                        m2_end    = annot[m2_trans]['end']
                        m2_strand = annot[m2_trans]['strand']
                        m2_gene   = annot[m2_trans]['gene_id']
                        if not m1_chr == m2_chr:
                            continue

                        dist = get_distance(m1_start,
                                            m1_end,
                                            m2_start,
                                            m2_end)
                        if dist <= maxlen:
                            counter += 1
                            writer.writerow([read,
                                             m1_gene,
                                             m1_trans,
                                             m1_chr,
                                             str(m1_start),
                                             str(m1_end),
                                             m1_strand,
                                             str(dist),
                                             m2_gene,
                                             m2_trans,
                                             m2_chr,
                                             str(m2_start),
                                             str(m2_end),
                                             m2_strand])        
                    else:
                        logging.error(f"read {read} is not a previous match")
                else:
                    logging.error(f"No annotation for {m2_trans}")

    if counter > 0:
        logging.debug (f"found {counter} reads as putative fusions")
    else:
        logging.error(f"No annotations found from {input_file}")
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

    annotation    = load_annotations(args.gtf)
    
    first_matches = load_matches(args.input1, annotation)
    
    filter_matches(args.output,
                   args.input2,
                   first_matches,
                   annotation,
                   args.maxlen)
    
    logging.debug("Completed")


if __name__ == "__main__":
    main()


