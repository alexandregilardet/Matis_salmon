#!/usr/bin/env python3

import argparse
import csv
import sys

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Creates FASTA file for every probe in probes_124SNPs.csv and its reverse
    complement as well as a csv compiling all probes for use of grep in shell programs""",
    epilog=f"Example of use: {sys.argv[0]} -p probes_124SNPs.csv -of probes_124SNPs_with_revcomp_fasta -on probes_124SNPs_with_revcomp.csv",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-p", "--in_probes_csv", type=str, required=True, help="input probes in csv"
)
parser.add_argument(
    "-of",
    "--output_folder",
    type=str,
    required=True,
    help="folder to write generated files",
)
parser.add_argument(
    "-on", "--out_name", type=str, required=True, help="csv output file"
)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
##############################################################################


def main():
    probes, revcomp_probes = get_probes()
    write_fasta_files(probes, revcomp_probes)
    csv_writer(probes, revcomp_probes)
    return None


def get_probes():
    """get probes from probes_124SNPs.csv file"""
    probes = {}
    revcomp_probes = {}
    with open(args.in_probes_csv, "r") as f_probes_in:
        for line in f_probes_in.readlines():
            if line[0] == "#":
                pass  # pass the information lines
            else:
                probe_id = line.split(";")[0]
                probe_seq = line.strip().split(";")[1]  # remove /n after seq
                rev_comp_seq = reverse_comp(probe_seq)
                probes[probe_id] = probe_seq
                new_probe_id = f"{probe_id}_reverse_complement"
                revcomp_probes[new_probe_id] = rev_comp_seq
    # return two dict with id and either sequence or reverse complement sequence
    return probes, revcomp_probes


def reverse_comp(sequence):
    """write reverse complement sequence"""
    my_dna = Seq(sequence, generic_dna)
    rev_comp_my_dna = str(my_dna.reverse_complement())
    rev_comp_my_dna = rev_comp_my_dna.replace("[", "]")
    rev_comp_my_dna = rev_comp_my_dna.replace("]", "[", 1)
    return rev_comp_my_dna


def write_fasta_files(probes, revcomp_probes):
    """write FASTA file for both the probe and its reverse complement"""
    out_folder = f"{args.output_folder}/"
    for probe in probes:
        fname = out_folder + probe + ".fasta"
        sequence = probes[probe]
        file_content = f">{probe}\n{sequence}"
        with open(fname, "w") as f_out:
            f_out.write(file_content)
    for revcomp_probe in revcomp_probes:
        fname = out_folder + revcomp_probe + ".fasta"
        revcomp_sequence = revcomp_probes[revcomp_probe]
        file_content = f">{revcomp_probe}\n{revcomp_sequence}"
        with open(fname, "w") as f_out:
            f_out.write(file_content)
    return None


def csv_writer(probes, revcomp_probes):
    """compile all probes in a csv for later grep"""
    f_out_name = args.out_name
    with open(f_out_name, "w", newline="") as f_out:
        w = csv.writer(f_out)
        for key, val in probes.items():
            w.writerow([key, val])
        for key2, val2 in revcomp_probes.items():
            w.writerow([key2, val2])


main()
