#!/usr/bin/env python3

import argparse
import sys
import pandas as pd

from collections import Counter

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Reads in csv of genotypes for each individual at certain probes and returns 
    frequencies of each allele""",
    epilog=f"Example of use: {sys.argv[0]} -g wild_farmed_Island_gts.csv",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-g", "--genotypes", type=str, required=True, help="csv with genotypes"
)

if len(sys.argv) == 1:  # help message when script called without enough argument
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
##############################################################################


def main():
    data = read_csv_as_panda()
    header = list(data)  # list of probes
    rows = data.shape[0]  # number of rows
    listdata = data.transpose().values.tolist()  # gives each column as a list
    # list of lists each containing all genotypes for a probe
    parse_gts(listdata, header, rows)
    return None


def read_csv_as_panda():
    data = pd.read_csv(args.genotypes, sep=";", header=0, index_col="individual")
    data = data.drop(["river", "short_name"], axis=1)  # drop columns not to be used
    return data


def parse_gts(listdata, header, rows):
    h = 0  # counter for header
    for probe in listdata:
        s = ","  # separator for join
        s = s.join(probe)  # string joined all genotypes separated by a , for a probe
        count = dict(Counter(s))  # obj dict
        #print(f'diploid state {Counter(probe)}')
        del count[","]
        #print(count)
        if count.get("N") == None:
            N = 0
        else:
            N = count.get("N")
            del count["N"]
        ind = (rows * 2) - N  # total number of alleles minus the N calls
        nd = list(count.keys())
        nb = list(count.values())
        if len(nd) == 2:
            #nd1 = nd[0]
            #nd2 = nd[1]
            #nb1 = nb[0]
            #nb2 = nb[1]
            #print(f"{header[h]} {nd1} {round(nb1/ind, 3)} {nd2} {round(nb2/ind, 3)} / {ind}")
            #round to only 3 decimals
            MA = min(count, key = count.get)
            nb = count.get(MA)
            print(f'{header[h]} {MA} {round(nb/ind, 3)} / {ind}')
        elif len(nd) == 1:
            nd1 = nd[0]
            nb1 = nb[0]
            print(f"{header[h]} {nd1} {round(nb1/ind, 3)} / {ind}")
            print('there is a fixed site')
        #print(header)
        print(" ")
        h += 1
    return None


main()
