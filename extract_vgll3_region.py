#!/usr/bin/env python3

import argparse
import sys

from Bio import SeqIO

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Extract genomic region from Ssa FASTA complete sequence through given boundaries""",
    epilog=f"Example of use: {sys.argv[0]} -r Ssa25_salmo_salar.fa -g vgll3 -s 28650000 -e 28670000",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-r", "--reference_fa", type=str, required=True, help="reference (FASTA)"
)
parser.add_argument(
    "-g", "--gene", type=str, help="name of region to extract", required=True
)
parser.add_argument(
    "-s", "--start", type=int, help="start of region to extract", required=True
)
parser.add_argument(
    "-e", "--end", type=int, help="end of region to extract", required=True
)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
##############################################################################


def get_extract():
    """gets the region from given boundaries and write new fasta file from name given"""
    sschr = SeqIO.read(args.reference_fa, "fasta")
    name = args.gene
    start = args.start - 1
    stop = args.end - 1
    # 28.650.000-28.670.000 from Ayllon et al. and Barson et al. 2015
    # AKGD00000000.4 NCBI Davidson et al. 2010
    extract = sschr[start:stop]
    SeqIO.write(extract, f"{name}_{args.start}_{args.end}.fa", "fasta")
    return None


get_extract()
