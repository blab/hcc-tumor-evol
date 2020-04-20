#script to sample cells tumor simulations output and convert to beast output

#load packages
library(tidyverse)

#intializatin
set.seed(192837)
setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")

#read in simulation data
CC_sim_cells <- read.csv("outputs/sim_CC_cells.csv")
CSC_sim_cells <- read.csv("outputs/sim_CSC_cells.csv")


#temp global var for developing functions
cells_df <- CC_sim_cells

#calculate parameters of simulations
calc_sim_stats <- function(cells_df) {
  stats_df <- data.frame("avg_mutation_rate" = mean(cells_df$mutation_rate),
                         "avg_proliferation_rate" = mean(cells_df$proliferation_rate),
                         "avg_death_rate" = mean(cells_df$alpha))
  return(stats_df)
}

CC_stats_df <- calc_sim_stats(CC_sim_cells)
CC_stats_df <- calc_sim_stats(CC_sim_cells)

#sample surviving cells
sample_alive_cells <- function(cells_df){
  endpoint <- max(cells_df$deathdate)
  alive_cells <- 
  
}
