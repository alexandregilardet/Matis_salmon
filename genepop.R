library(genepop)

input <- "C:/Users/alexa/Matis/results/08_05_20/genepop_old_data.txt"
output <- "C:/Users/alexa/Matis/results/08_05_20/test_diff_old_data"
setting <- "C:/Users/alexa/Matis/results/08_05_20/test_diff_old_data_settings"

#option3 population differentiation
#suboption4 genotypic difference for all pairs of
#populations
test_diff(
  inputFile = input,
  genic = FALSE,
  pairs = TRUE,
  outputFile = output,
  settingsFile = setting,
  dememorization = 10000,
  batches = 100,
  iterations = 5000,
  verbose = interactive()
)

#option5 basic information
#suboption1 allele and genotype frequencies
basic_info(inputFile = input, 
           outputFile = output, 
           verbose = interactive())

