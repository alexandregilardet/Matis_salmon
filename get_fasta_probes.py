#!/usr/bin/env python3

import argparse
import sys

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Creates FASTA file for every probe in SSaTrack_Annotation_for_sharing.csv file""",
    epilog=f"Example of use: {sys.argv[0]} -p SsaTrack_Annotation_for_sharing.csv",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument("-p", "--in_affy_csv", type=str, required=True, help="Probes")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
##############################################################################


def main():
    probes_seq = get_probes()
    write_fasta_files(probes_seq)
    return None


def get_probes():
    """get probes from SSaTrack_Annotation_for_sharing.csv file"""
    probes_seq = {}
    with open(args.in_affy_csv, "r") as f_affy_in:
        for line in f_affy_in.readlines():
            if line[0] == "#":
                pass  # pass the information lines
            elif line.replace('"', "").split(";")[0] == "Probe Set ID":
                pass  # pass the header line
            else:
                probe_set_id = line.replace('"', "").split(";")[0]
                probe_seq = line.replace('"', "").split(";")[6]
                probes_seq[probe_set_id] = probe_seq
    return probes_seq


def write_fasta_files(probes_seq):
    out_folder = "all_probes_fasta/"
    for probe in probes_seq:
        fname = out_folder + probe + ".fasta"
        sequence = probes_seq[probe]
        file_content = f">{probe}\n{sequence}\n"
        with open(fname, "w") as f_out:
            f_out.write(file_content)
    return None


main()
