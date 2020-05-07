library(genepop)

#option3 population differentiation
#suboption4 genotypic difference for all pairs of
#populations
test_diff(
  inputFile = "C:/Users/alexa/Matis/results/30_04_20/genepop_cor_seage.txt",
  genic = FALSE,
  pairs = TRUE,
  outputFile = "C:/Users/alexa/Matis/results/05_05_20/test_diff_cor_seage",
  settingsFile = "C:/Users/alexa/Matis/results/05_05_20/test_diff_cor_seage_settings",
  dememorization = 10000,
  batches = 100,
  iterations = 5000,
  verbose = interactive()
)

#option5 basic information
#suboption1 allele and genotype frequencies
basic_info(inputFile = "C:/Users/alexa/Matis/results/30_04_20/genepop_cor_seage.txt", 
           outputFile = "C:/Users/alexa/Matis/results/30_04_20/gt_freq_cor_seage", 
           verbose = interactive())

