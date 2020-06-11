#!/usr/bin/env python3

# Credits: Dr Sæmundur Sveinsson

import argparse
import sys

from collections import Counter

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Reads in the SNP calls from the salmon Affymetrix chip, the reference
    file (SsaTrack_Annotation_for_sharing.csv) and information from the specific probes. 
    Outputs a file with the genotype of each individual at the specified positions.""",
    epilog=f"Example of use: {sys.argv[0]}\
    -a ../../data/Stofnfiskur samples/SsaTrack_Annotation_for_sharing.csv\
    -g  ../../data/Stofnfiskur samples/SsaTrack_Island_1-8.calls_mod.txt\
    -i probes_sex.csv\
    -o sex_gts.csv",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-a", "--annotation", type=str, required=True, help="Annotation file (csv)"
)
parser.add_argument(
    "-g",
    "--genotypes",
    type=str,
    required=True,
    help="Genotype calls from CIAGENE (tsv file)",
)
parser.add_argument(
    "-i", "--in_probes", type=str, help="List of probes of interest", required=True
)
parser.add_argument(
    "-o", "--out_fname", type=str, required=True, help="Prefix for output file"
)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
##############################################################################


def main():
    probe_info = get_probes()
    probe_genotypes = get_probe_genotypes(probe_info)
    write_gt_output(probe_info, probe_genotypes)


def get_probes():
    # reads in the in_probes file and returns a dict with probes to use and relative info
    probes = {}
    header = []
    data_lines = []
    with open(args.in_probes, "r") as f_in:
        for line in f_in.readlines():
            if line.split(" ")[0] == "###":
                pass
            elif line.split(" ")[0] == "#":
                header = line.strip().split(";")
            else:
                data_lines.append(line.strip())
    # print(header)
    # print(data_lines)
    probes = parse_probe_info(data_lines, header)
    return probes


def parse_probe_info(data_lines, header):
    probes = {}
    for data_line in data_lines:
        probe_dict = {}
        probe, chrom, source = data_line.split(";")
        # probe_dict[probe] = probe
        probe_dict["chrom"] = chrom
        probe_dict["source"] = source
        probes[probe] = probe_dict
        # embeded dict with probe as main key, then chrom and source as main values
        refernce_alleles = get_probe_allele(probe)
        probe_dict["0"] = refernce_alleles[0] + refernce_alleles[0]
        probe_dict["1"] = refernce_alleles[0] + refernce_alleles[1]
        probe_dict["2"] = refernce_alleles[1] + refernce_alleles[1]
        probe_dict["-1"] = "NN"
        # add all possible alleles combinations
    # print(probes)
    # print(len(probes))
    return probes


def get_probe_genotypes(probe_info):
    # Reads in the raw SNP genotype file and parses information into dictionaries
    probes = probe_info.keys()
    individuals = get_individuals()
    # print(len(probes))
    genotypes = parse_genotypes(probes, individuals)
    return genotypes


def get_individuals():
    # Returns a list of individuals in the same order as in genotype file
    individuals = []
    with open(args.genotypes, "r") as f_in:
        for line in f_in.readlines():
            if line.split("\t")[0] == "probeset_id":
                individuals = line.strip().replace(".CEL", "").split("\t")[1:]
                # print(individuals, len(individuals))
                return individuals


def parse_genotypes(probes, individuals):
    # Reads in the raw SNP genotype file and parses information into dictionaries
    probe_individual_gts = {}
    # print(probes)
    with open(args.genotypes, "r") as f_in:
        for line in f_in.readlines():
            if line[0] == "#":
                pass
            elif line.split("\t")[0] == "probeset_id":
                pass
            elif line.split()[0].strip() in list(probes):  # probes of interest
                probe_id = line.split()[0]
                probe_individual_gts[probe_id] = parse_gt_line(
                    probe_id, individuals, line
                )
    # print(probe_individual_gts)
    # print("Retrieved", len(probe_individual_gts.keys()), "probes")
    return probe_individual_gts


def parse_gt_line(probe_id, individuals, data_line):
    # Parses the data line an returns a dict with the genotypes and information per probe
    line_data_dict = {}
    gts = data_line.split()[1:]
    # Run check
    if len(individuals) != len(gts):
        sys.exit(f"Somethin wrong with data in {probe_id}, check data")
    # print(individuals)
    individ_gt = {}
    # print("gt_length", len(gts))
    for i in range(0, len(gts)):
        individual = individuals[i]
        gt = gts[i]
        individ_gt[individual] = gt
    line_data_dict[probe_id] = individ_gt
    return line_data_dict


def parse_ref_allele(data_line):
    # reads a data line from the SsaTrack Annotation file and outputs a list with two values
    probe = data_line.split(",")[0]
    ref_allele_A = data_line.split(",")[7]
    ref_allele_B = data_line.split(",")[8]
    # Check if this makes senes
    if ref_allele_A == ref_allele_B:
        sys.exit(f"Allele A and B are the same, check annotation file, {probe}")
    check_ref_alleles(ref_allele_A, probe)
    check_ref_alleles(ref_allele_B, probe)
    return [ref_allele_A, ref_allele_B]


def get_probe_allele(probe_id):
    # Reads the SsaTrack_Annotation_for_sharing.csv file and returns the corresponding alleles
    with open(args.annotation, "r") as f_in:
        for line in f_in.readlines():
            if line.split(",")[0] == probe_id:
                ref_alleles = parse_ref_allele(line)
    # print(ref_alleles)
    return ref_alleles


def check_ref_alleles(ref_allele, probe):
    if len(ref_allele) == 1:
        pass
    elif ref_allele in "ATCG":
        pass
    else:
        sys.exit(f"Something wrong with reference allele ({ref_allele}) of {probe}")
    return None


def write_gt_output(probe_info, probe_genotypes):
    # write the output file
    probes = probe_info.keys()
    # print(probes)
    individuals = get_individ_from_dict(probe_genotypes, probes)
    river_dict = get_river_dict(individuals)
    new_individ_name_dict = make_individ_name_dict(individuals, river_dict)
    gt_output_dict = get_gt_output_dict(
        probe_info, probe_genotypes, individuals, probes
    )
    # probes_in_file = list(gt_output_dict[list(individuals)[0].keys()])
    # print(probes_in_file)
    with open(args.out_fname, "w") as f_out:
        probes_header = ";".join(gt_output_dict["header"])
        f_out.write(f"river;individual;short_name;{probes_header}\n")
        for individual in individuals:
            # new_individ_name = change_individ_nama(individuals)
            river = river_dict[individual]
            gt = gt_output_dict[individual]
            f_out.write(
                f"{river};{individual};{new_individ_name_dict[individual]};{gt}\n"
            )
    return None


def get_individ_from_dict(probe_genotypes, probes):
    # get individual
    # print(probe_genotypes)
    probe = list(probes)[0]  # need to set to a probe for later access in dict?
    individuals = probe_genotypes[probe][probe].keys()
    # print(individuals)
    return individuals


def get_river_dict(individuals):
    # get river from individual name
    river_dict = {}
    for individual in individuals:
        river_short = individual.split("_")[0]
        river_dict[individual] = river_short
    # print(river_dict)
    return river_dict


def make_individ_name_dict(individuals, river_dict):
    odd_names = {"Hyb": "BoT"}
    individ_name_dict = {}
    river_counts = get_river_counts(individuals, river_dict)
    river_numbers = get_river_numbers(individuals, river_dict, river_counts)
    x = 0
    for individual in individuals:
        individual_id = ""
        river = individual.split("_")[0]
        river_individ_numbers = river_numbers[x]
        # river_count = river_counts[river]
        if river in odd_names.keys():
            individual_id = odd_names[river]
        else:
            individual_id = individual.split("_")[1]
        individ_name_dict[
            individual
        ] = f"{river}{individual_id[0:4]}-{river_individ_numbers}"
        x += 1
    return individ_name_dict


def get_river_numbers(individuals, river_dict, river_counts):
    # river_numbers = {}
    rivers = list(set(river_dict.values()))
    rivers.sort()
    rivers_numbers_list = []
    for river in rivers:
        river_count = river_counts[river]
        numbers = []
        if river == "Hyb":  # only 1 Hyb individual
            numbers = ["01"]
        else:
            for i in range(0, river_count):
                number = f"{i+1:02}"
                numbers.append(number)
        # print(river, river_numbers)
        rivers_numbers_list.append(numbers)
    print(rivers_numbers_list)
    rivers_flat_list = []
    for sublist in rivers_numbers_list:
        for item in sublist:
            rivers_flat_list.append(item)
    print(rivers_flat_list, len(rivers_flat_list))
    return rivers_flat_list


def get_river_counts(individuals, river_dict):
    rivers = []
    odd_names = {"Hyb": "HybBoT", "KroF172679": "Kro", "KroF172703": "Kro"}
    for individual in individuals:
        river = river_dict[individual]
        if river in odd_names.keys():
            rivers.append(odd_names[river])
        else:
            rivers.append(river)
    river_counts = Counter(rivers)
    # print(rivers, len(rivers))
    return river_counts


def get_gt_output_dict(probe_info, probe_genotypes, individuals, probes):
    gt_output_dict = {}
    missing_probes = []
    # print(len(probe_genotypes.keys()))
    for individual in individuals:
        probes_in_file = []
        individ_gt = []
        # print(individual)
        for probe in probes:
            probe_individ_gt = ""
            # print(probe)
            # print(probe_genotypes)
            # print(probe_genotypes[probe][probe][individual])
            if probe_genotypes.get(probe) != None:
                snp_call = probe_genotypes.get(probe).get(probe).get(individual)
                probe_individ_gt = probe_info[probe][snp_call]
                # print(probe_individ_gt)
                # individ_gt[individual][probe] = probe_individ_gt
                individ_gt.append(probe_individ_gt)
                probes_in_file.append(probe)
                pass
            elif probe_genotypes.get(probe) == None:
                missing_probes.append(probe)
        # print(individual, ";".join(individ_gt))
        gt_output_dict[individual] = ";".join(individ_gt)
        # print(individ_gt)
    print(
        f"#Warning {len(set(missing_probes))} probes are missing from {args.genotypes}. They are the following:"
    )
    for missing_probe in set(missing_probes):
        print(f"#{missing_probe}")
    gt_output_dict["header"] = probes_in_file
    # print(probes_in_file)
    return gt_output_dict


main()
