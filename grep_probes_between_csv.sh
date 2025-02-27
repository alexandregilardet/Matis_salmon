#!/bin/sh

### get all the common probes between shorter sequences in SsaTrack_Annotation_for_sharing.csv ###
### grepped on longer sequences in probes_124SNPs_with_revcomp.csv                             ###

#grep -v "#" ../../data/Stofnfiskur\ samples/SsaTrack_Annotation_for_sharing.csv | grep -v "Probe" |
#for probe in $(cut -f 7 -d ";"); do
#    echo "Grepping $probe";
#    grep -F $probe probes_124SNPs_with_revcomp.csv |
#    tee -a common_probes.txt;
#done

#-v grep all except this
#-f cut column
#-F interprets patterns literally and not as regular expression
#-a forces tee to append instead of overwrite


### get common probes from sea age probes/sex specific probes ###
### grepped in probes_124SNPs_with_revcomp.csv                ###

SEA_AGE_PROBES="CCAATGGCTCTTATCTCCATGGAAACCAGCGTGAT[A/G]TGGCCTCCCAGACTGTCCCCGCCCACAAGAGAAA
CTTGGTGTTTTTATGTTTACCTTGTCCACAGGTGA[A/C]CCCGTTGAGCATGTGCGCTTTAGTGAACATGTGT
GGATAGAGCAGCGGCTGCGAAAAATAGGTATGTAT[A/G]CTGCAGTCTGTAGCTTACCTTTGTATTTTATCAA
AATTTCTTTTGCCTTTTTGAAGTTTGAAATAAAAT[G/T]TCTAGTTTGCAGAATGTACCACTAAATAAATGGG
GTTGCTTGTGTTTGTGTGTTTACTATCATAAAATA[G/T]CATTTTTGGGCTCCTAGCCTGGCACTTGTGCATC
TTACAAATAATATAATACAGTGTGGTGCTCTGTAT[A/G]GTTAGTTTTATTTTCCACTATAATATTTAACTTC
AATAGTTTCTCCTCTGTTGTCATCCAGAATTAATC[A/C]GATTGTATTCCTCCAGTACAGAACAGCTGTGTGG
TCTTTACAGACCGAAGTTTATCCAAAAAGTTGTAG[A/G]ACGTTTCCGCGCGTGTGTTCTGGTGAAAGTATGA
TATGCTGTCTGACAAAATCACAATTTTAGTAGTTC[C/T]TCAAAGTAAATAATGTAAGTCATGCATTTACATCT
AGCAGCCAGAGCCCAGGGATACACAGTGATGTGAA[C/G]ACAGAGGTGGAGCAGGCCAGCCCCCCCACCCTGAC
TTATAGAATGGCAAGTCGCGTTAAACAACAAAGGA[G/T]CTGATTGGCTGACACTGCCATTCCAAAATTGCTG
TCATTCCCTTTTATTCCTTCAGGGATAGACATAAA[A/C]ACAAGAATATGTTAAGAGAAGCTGAAAACTTCAC
GGGCAGATCAGAACACCTTTGATTAAAATACCAAC[A/C]CAGATAATCAACATGTGGGGGTGGATTCTAATGA
GTGACTGACAGATGTTTCTGAGGTTCTCTGTCACA[C/T]GTTCAGAAGGGGGTTGATCAACTGCACCTGCTTC
GTCCCTGCCAACACAAACCCATCTCTTCTCATCTC[A/C]AGCCTGGGCCTAAAGAATTTGCCATTACTTAAAA
TGCCCAATTAAACCGGCCATGAAATGACACAGACT[A/G]TTCACCAGGCCCACAGAGGCTGAAAGATTGAGAT
CTGTACCCTTCAGCCTCACTGTATTTCAATGACAC[A/C]GTCTTCAAGGTTAGCAGATTGGAGGGTCCCATTA
TGTTTGGTTTGGGCGGCTTTCCCATGTCCAGTAAC[C/T]TTAGGTTTATTTTGCTTGCTACATCAGTTGTTTT
ATGCTCTATTGTCCCTTATGTACTGAAGAATAGTT[A/C]ATGTGAGTGTGGGATCATGATTTACAAGGCCATT
TACGGGGAGGGCTGCGAGTTCAGCGGATGGTTTGG[A/G]GCACACCTGTCCTCAGCTCCCTTCATGCATGGTC
"

SEX_PROBES="TGCTTGCCTGAGTTACCTAAGGCGGAAGACAATTC[C/T]ATGTAGTCATTGCTCTATACAGTACTGTGCGTCTC
TTATTTATGAGAGGTATTGGCATGTTAAATGCACC[A/G]AGCTGTCTGTCTAAACTACTGGCACACAGCTCAGA
AACAGCTTGGTGCATTCAACATGTCAATACCTCTC[A/G]TAAAAATAAGTAGTGTTGAAGTCGATCTCTCCTCC
TTCAATCTCTCTTCCACTTTGAGCCATGAGAGATT[A/G]ACATGAATATCATTATTGTTAGCTCCCTGTGTATT
CCATGCTTTAATTTCTTTGATAAATTTACGTCATT[C/T]AGCAGACGCTCTTATCCAGAGCGACTTACAAATAG
TAAGTCGCTCTGGATAAGAGCGTCTGCTAAATGAC[G/T]TAAATGTAAATGTCAGAACCTGAAGGACCTGGTGC
TCACATTTTGTTACATTACAGCCTTATTCTAAAAT[G/T]GAGTAAATAAATAAACATCTTCATCAATCTACACA
CTCTGGTCTGATGAAACCAACATTGAACTCTTTGG[A/C]CTGAATACCAAGAATCACATCAGGAGGAAATCTGG
CGTAGGAATTATGCCAGGTTTGCTTCAGCGGTGAC[A/G]CTTGGCATTCAGGCCAAAGAGTTCCATCTTGGTTT
AGCATGCTGTTGCCACCACCATGCTTCACTGTTGG[G/T]ATGATGCCAGATTTCTTCCATATGTGACGCTTGGC
TGGTGAGTACTCCCATTCTACTCAGCAGATCCTCT[A/C]AAACTCTATCAGGTTGGATGGGGAGCATCGCCGCA
CACCTTTGACAGCGATTGCATCCTCTTGGGTATGA[C/T]GCTACAAGCTTGGCACACCTGTATTTGGGGAGTTT
CTTTGTTGAAGCACCTTTAGCAGCGACTACAGCCT[C/T]GAGTCTTTTTGTGTATAACGTTACAAGCTTCGCAC
CTTGTAACGTCATATCCAAAAAGGCTGTAATCGCT[G/T]CCAAACGTACTTCAACAAAGTACTGAGTAAAGGGT
TCATACCCAAGAATACTTGAGGCTGTACTCGCTGC[C/T]TAAGTTGCTTGAACTGAGTACTGAGTAAAGGGTCT
AGCAAAAACAGTTTACAAATAAAACACTTTAATAT[C/T]ACATTTACATAAGTATTCAGACCCTTTACTCAGAA
ATGTGATATTCCAGTTTTAACTTTTCTATAAATTA[G/T]CAAAAATGGCTAAAAACCTGTTTTTGCTTTGTCAG
ATGTTTATTTATTTTTTATGAATGAGCAAATGTTC[G/T]AAAAAACTGTTTTTGCTTTGTCATTATGGGGTATT
"

for probe in $SEX_PROBES; do
    echo "Grepping $probe"
    grep -F $probe probes_124SNPs_with_revcomp.csv |
    tee -a common_probes_sex_specific.txt
done 

### get the input probes grepped in another file ###

#grep -v "#" ../../data/Stofnfiskur\ samples/SsaTrack_Annotation_for_sharing.csv | grep -v "Probe" |
#for probe in $(cut -f 7 -d ";"); do
#    echo "Grepping $probe";
#    if grep -F $probe probes_124SNPs_with_revcomp.csv; then
#        echo $probe |
#        tee -a common_probes_input.txt;
#    fi
#done
