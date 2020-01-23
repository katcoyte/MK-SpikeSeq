# Test for temporal structure in each high-frequency infant time series
# Katharine Z. Coyte, January 2020
#
# Here we follow the approach and R package developed by Faust et al, originally published in:
# [1] Faust, K. et al. Signatures of ecological processes in microbial community time series. Microbiome 6, 120 (2018). 
#
# Note, this analysis requires the first steps of python notebook SIFig_noise_analysis to have been run to generate the time series data for each infant (new_bacteria_otu_data_for_r_X.csv)

library(devtools)
library(seqtime)
library(WrightFisher)
library(stinepack)
library(ggplot2)

### BACTERIA ###

# Initialize data structures
unique_babies <- c('035', '073', '091', '127', '128', '139', '158', '175', '203','260', '307', '374', '320')
num_unique_babies <- length(unique_babies)
save_data_bacteria <- data.frame()
save_noise_bacteria <- data.frame()
save_noise_small_bacteria <- data.frame()

# As in [1], for each infant we first interpolate the timeseries for each constituent ASV / genus to ensure equidistant time intervals. 
# We then identify noise type using the identifyNoisetypes function from the seqtime package (we also perform a neutrality test, but do not analyse
# these results due to the high level of white noise observed, as advised by K. Faust, personal correspondence)

for(baby_ix in 1:num_unique_babies){
  
  load_string <- paste('new_bacteria_otu_data_for_r_', unique_babies[baby_ix],'.csv', sep="")
  
  for_faust <- read.csv(load_string)
  my_col_names <- colnames(for_faust)
  my_col_names <- my_col_names[4:length(my_col_names)]
  t_vec <- for_faust$day
  new_t_vec <- seq(t_vec[1],t_vec[length(t_vec)],by=1)
  
  num_species = length(my_col_names)
  
  output <- matrix(ncol=length(new_t_vec), nrow=num_species)
  for(i in 1:num_species){
    cur_species <- my_col_names[i]
    cur_y <- for_faust[[cur_species]]
    new_y <- stinterp(t_vec,cur_y,new_t_vec)
    output[(i),] <- new_y$y
  }
  
  output_relative <- apply(output, 2, function(x) x/sum(x))
 
  save_data_bacteria[unique_babies[baby_ix], 'neutrality'] <- -log10(NeutralCovTest(t(output), tm=t_vec))
  
  my_noise_types=identifyNoisetypes(output, smooth=TRUE)
  
  plotNoisetypes(my_noise_types)
  
  save_noise_small_bacteria[unique_babies[baby_ix], 'black'] <- length(my_noise_types$black)
  save_noise_small_bacteria[unique_babies[baby_ix], 'brown'] <- length(my_noise_types$brown)
  save_noise_small_bacteria[unique_babies[baby_ix], 'pink'] <- length(my_noise_types$pink)
  save_noise_small_bacteria[unique_babies[baby_ix], 'white'] <- length(my_noise_types$white)
  print(baby_ix)
  
  # Store proportions of different noise types
  # black
  for(noise_ix in my_noise_types$black){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_bacteria[unique_babies[baby_ix], cur_noisey_species] <- 3
  }
  # brown
  for(noise_ix in my_noise_types$brown){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_bacteria[unique_babies[baby_ix], cur_noisey_species] <- 2
  }
  # pink
  for(noise_ix in my_noise_types$pink){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_bacteria[unique_babies[baby_ix], cur_noisey_species] <- 1
  }
  # white
  for(noise_ix in my_noise_types$white){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_bacteria[unique_babies[baby_ix], cur_noisey_species] <- 0
  }
}

save_data_bacteria['n']='n'
fill <- "gold1"
line <- "goldenrod2"
p10 <- ggplot(save_data_bacteria, aes(x = n, y = neutrality)) +
  geom_boxplot(fill = fill, colour = line) +
  scale_y_continuous(name = "-log10(p-value)",
                     breaks = c(round(-log10(0.05),3), 5, 8)) +
  scale_x_discrete(name = "Neutrality") +
  ggtitle("P-vals of neutrality test") + theme(plot.title = element_text(hjust = 0.5))
p10

save_noise_bacteria[is.na(save_noise_bacteria)] <- 0
#write.csv(save_noise_bacteria, "save_noise_bacteria.csv")
#write.csv(save_noise_small_bacteria, "save_noise_small_bacteria.csv")
# write.csv(save_data_bacteria, "save_data_bacteria.csv")



### FUNGI ###

save_data_fungi <- data.frame()
save_noise_fungi <- data.frame()
save_noise_small_fungi <- data.frame()

for(baby_ix in 1:num_unique_babies){
  
  load_string <- paste('new_fungi_otu_data_for_r_', unique_babies[baby_ix],'.csv', sep="")
  
  for_faust <- read.csv(load_string)
  my_col_names <- colnames(for_faust)
  my_col_names <- my_col_names[4:length(my_col_names)]
  t_vec <- for_faust$day
  new_t_vec <- seq(t_vec[1],t_vec[length(t_vec)],by=1)
  
  num_species = length(my_col_names)
  output <- matrix(ncol=length(new_t_vec), nrow=num_species)
  
  for(i in 1:num_species){
    cur_species <- my_col_names[i]
    cur_y <- for_faust[[cur_species]]
    new_y <- stinterp(t_vec,cur_y,new_t_vec)
    output[(i),] <- new_y$y
  }
  
  output_relative <- apply(output, 2, function(x) x/sum(x))
  
  save_data_fungi[unique_babies[baby_ix], 'neutrality'] <- -log10(NeutralCovTest(t(output), tm=t_vec))
  
  my_noise_types=identifyNoisetypes(output, smooth=TRUE)
  plotNoisetypes(my_noise_types)
  
  save_noise_small_fungi[unique_babies[baby_ix], 'black'] <- length(my_noise_types$black)
  save_noise_small_fungi[unique_babies[baby_ix], 'brown'] <- length(my_noise_types$brown)
  save_noise_small_fungi[unique_babies[baby_ix], 'pink'] <- length(my_noise_types$pink)
  save_noise_small_fungi[unique_babies[baby_ix], 'white'] <- length(my_noise_types$white)
  print(baby_ix)
  
  
  # black
  for(noise_ix in my_noise_types$black){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_fungi[unique_babies[baby_ix], cur_noisey_species] <- 3
  }
  # brown
  for(noise_ix in my_noise_types$brown){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_fungi[unique_babies[baby_ix], cur_noisey_species] <- 2
  }
  # pink
  for(noise_ix in my_noise_types$pink){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_fungi[unique_babies[baby_ix], cur_noisey_species] <- 1
  }
  # white
  for(noise_ix in my_noise_types$white){
    cur_noisey_species <-my_col_names[noise_ix]
    save_noise_fungi[unique_babies[baby_ix], cur_noisey_species] <- 0
  }
}

save_data_fungi['n']='n'
fill <- "gold1"
line <- "goldenrod2"

p10 <- ggplot(save_data_fungi, aes(x = n, y = neutrality)) +
  geom_boxplot(fill = fill, colour = line) +
  scale_y_continuous(name = "-log10(p-value)",
                     breaks = c(round(-log10(0.05),3), 5, 8)) +
  scale_x_discrete(name = "Neutrality") +
  ggtitle("P-vals of neutrality test") + theme(plot.title = element_text(hjust = 0.5))
p10

save_noise_fungi[is.na(save_noise_fungi)] <- 0
#write.csv(save_noise_fungi, "save_noise_fungi.csv")
#write.csv(save_noise_small_fungi, "save_noise_small_fungi.csv")
#write.csv(save_data_fungi, "save_data_fungi.csv")








