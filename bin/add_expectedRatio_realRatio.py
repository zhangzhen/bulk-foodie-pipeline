#! /usr/bin/env python

import sys
import pysam
import pickle
import argparse

compl_dict = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}

def compl(string):
    return_str = ""
    for i in string[::-1]:
        return_str += compl_dict[i]
    return return_str

def main():
    parser = argparse.ArgumentParser(description="Add expected and real ratio to input file.")
    parser.add_argument("input_file", help="Input file path")
    parser.add_argument("--output_file", required=True, help="Output file path")
    parser.add_argument("--expected_ratio_file", required=True, help="Expected ratio file path (txt)")
    parser.add_argument("--genome_fasta", required=True, help="Genome FASTA file path")
    args = parser.parse_args()

    expected_ratio_dict = {}
    with open(args.expected_ratio_file) as f:
        for l in f:
            l = l.strip().split('\t')
            context, exp_con, rep_uncon, exp_ratio = l[0], int(l[1]), int(l[2]), float(l[3])
            expected_ratio_dict[context] = (exp_con, rep_uncon, exp_ratio)

    fa_file = pysam.FastaFile(args.genome_fasta)

    with open(args.input_file, "r") as fin, open(args.output_file, "w") as fout:
        for l in fin:
            ch, st, ed, con, uncon, ratio, CorG = l.strip().split("\t")
            st, ed, con, uncon, ratio = int(st), int(ed), int(con), int(uncon), float(ratio)
            ref_seq_by_chr = fa_file.fetch(reference=ch, start=st-4, end=st+5).upper()
            ref_seq_by_chr_compl = compl(ref_seq_by_chr)
            C_G = ref_seq_by_chr[4]
            exp_con, exp_uncon, exp_ratio = -1, -1, "NA"
            if C_G == "C" and "N" not in ref_seq_by_chr:
                exp_con, exp_uncon, exp_ratio = expected_ratio_dict[ref_seq_by_chr[:7]]
            elif C_G == "G" and "N" not in ref_seq_by_chr:
                exp_con, exp_uncon, exp_ratio = expected_ratio_dict[ref_seq_by_chr_compl[:7]]
            elif C_G != "C" and C_G != "G":
                raise ValueError
            else:
                continue
            real_ratio = ratio / exp_ratio
            fout.write("\t".join([
                ch, str(st), str(ed), str(con), str(uncon), ref_seq_by_chr, str(ratio), str(exp_ratio), str(real_ratio)
            ]) + "\n")

if __name__ == "__main__":
    main()

    
