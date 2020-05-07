#!/usr/bin/env python3

import argparse
import sys
import glob

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Script to convert the axiom probe sequences to ambigutity codes 
    in order to be able to map to reference sequence in Geneious""",
    epilog=f"Example of use: {sys.argv[0]} -f all_probes_fasta -o all_probes_with_ambiguity_fasta",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-f", "--fasta_folder", type=str, required=True, help="input probes (fasta)"
)
parser.add_argument(
    "-o",
    "--output_folder",
    type=str,
    required=True,
    help="folder to write generated files",
)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
##############################################################################

AMBIG_DICT = {
    "A/C": "M",
    "A/G": "R",
    "A/T": "W",
    "C/G": "S",
    "C/T": "Y",
    "G/T": "K",
    "-/A": "A",
    "-/AT": "AT",
}

# if len(sys.argv[0:]) != 1:
#    sys.exit(f'Script takes exactly 1 argument, probe fasta file with SNP in [/] format')
# probe_fasta = sys.argv[1]
# was first written to perform it on only one FASTA probe at a time


def main():
    path = f"{args.fasta_folder}/*.fasta"
    folder = glob.glob(path)
    for probe in folder:
        with open(probe, "r") as f_in:
            for line in f_in.readlines():
                if line[0] == ">":
                    header = line.strip()
                    name = header.replace(">", "")
                else:
                    seq = add_ambig(line.strip())
        # print(f"{header}\n{seq}")
        write_fasta_files(name, header, seq)


def add_ambig(seq_line):
    snp = seq_line.split("[")[1].split("]")
    amigb_code = AMBIG_DICT[snp[0]]
    sequence = seq_line.split("[")[0] + amigb_code + seq_line.split("]")[1]
    return sequence


def write_fasta_files(name, header, seq):
    out_folder = f"{args.output_folder}/"
    fname = out_folder + name + ".fasta"
    file_content = f"{header}\n{seq}\n"
    with open(fname, "w") as f_out:
        f_out.write(file_content)
    return None


main()
