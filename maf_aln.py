#!/usr/bin/env python3

import argparse
import csv
import sys

from Bio import AlignIO
from collections import Counter

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Reads in a nucleotide alignment (FASTA format) and outputs 
    the minor allele frequencies (MAF) of all variable sites """,
    epilog=f"Example of use: {sys.argv[0]} -a SSaTrack_Island_fullMT.fa\
    -o SSaTrack_Island_fullMT_MAF.csv",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-a", "--alignment", type=str, required=True, help="Alignment (FASTA)"
)
parser.add_argument(
    "-o", "--out_name", type=str, required=True, help="Prefix for output file"
)

if len(sys.argv) == 1:  # help message when script called without enough argument
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
##############################################################################


def main():
    mt_alignment = get_mt_alignment()
    nd_pos = get_nd(mt_alignment)
    aln_maf = get_freq(nd_pos)
    csv_writer(aln_maf)
    return None


def get_mt_alignment():
    # reads in the alignment and returns it to a variable
    alignment_file = args.alignment
    mt_alignment = AlignIO.read(alignment_file, "fasta")
    return mt_alignment


def get_nd(mt_alignment):
    # reads in alignment to get all possible nucleotides at each position by removing any N
    nd_pos = {}
    aln_length = mt_alignment.get_alignment_length()
    for i in range(0, aln_length):  # reads in all positions of the alignment
        nt_list = []
        for aln in mt_alignment:  # reads all nucleotides possible at this position
            nt_list.append(aln[i])
        new_nt_list = [y for y in nt_list if y != "N"]  # remove N
        nd_pos[i] = new_nt_list
    return nd_pos


def get_freq(nd_pos):
    aln_maf = {}
    # counts nucleotides at each position and returns a dict of the one with polymorphism
    for i in nd_pos:
        if len(Counter(nd_pos[i])) == 1:
            pass
        else:
            count = Counter(nd_pos[i])
            nd, maf = zip(*count.items())
            # aln_maf[i + 1] = (nd[0], maf[0], nd[1], maf[1]) to get both alleles
            if maf[0] > maf[1]:
                aln_maf[i + 1] = (nd[1], maf[1])  # +1 to set the correct position
            else:
                aln_maf[i + 1] = (nd[0], maf[0])
    return aln_maf


def csv_writer(aln_maf):
    f_out_name = args.out_name
    with open(f_out_name, "w", newline="") as f_out:
        w = csv.writer(f_out)
        for key, val in aln_maf.items():
            w.writerow([key, val])


main()
