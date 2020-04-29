#!/usr/bin/env Rscript

#make xml file for spatial sim
setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

# print(args[1])
# print(args[2])
# print(args[3])

normalize_locs <- function(cells){
  center_x <- min(cells$locx) + (max(cells$locx) - min(cells$locx))/2
  center_y <- min(cells$locy) + (max(cells$locy) - min(cells$locy))/2
  norm_cells <- cells %>% 
    mutate("norm_locx" = locx - center_x, "norm_locy" = locy - center_y)
  return(norm_cells)
}


CC_sim_cells <- read.csv(args[1]) %>% 
  normalize_locs

final_sampled_cells <- read.csv(args[2])
# CC_sim_cells <- read.csv("outputs/sim_CC_cells.csv") %>% 
#   normalize_locs

# final_sampled_cells <- read.csv("outputs/CC_final_sampled_cells.csv")

#function to find all unique mutations

isolate_muts <- function(i, mutations_vec, mutations_column) { #recurisve function
  
  #for each row split and filter to unique mutations in that cell
  row_muts <- na.omit(
    as.numeric(
      unlist(
        strsplit(as.character(mutations_column[i]),"[^0-9]")
      )
    )
  )
  
  mutations_vec <- c(mutations_vec, row_muts)
  
  if (i == length(mutations_column)) {
    #return final mutations_vec when get to the last row
    return(sort(unique(mutations_vec)))
    
  } else {
    return(isolate_muts(i=i+1, mutations_vec, mutations_column))
  }
}

get_row_muts <- function(mutations_row) {
  
  row_muts <- na.omit(
    as.numeric(
      unlist(
        strsplit(as.character(mutations_row),"[^0-9]")
      )
    )
  )
  return(row_muts)
  
}

compare_muts <- function(row_muts, all_muts) {
  return(as.numeric(all_muts %in% row_muts))
}

#define presence absence of muts in each cell
define_mut_presence_absence <- function(sampled_cells) {
  #find all mutations in sample
  all_muts <- isolate_muts(i = 1, mutations_vec = c(),
                           mutations_column = sampled_cells$mutations)
  rowise_muts <- purrr::map(sampled_cells$mutations, get_row_muts)
  
  mut_presence_absence <- purrr::map(rowise_muts, ~compare_muts(.,all_muts)) %>% 
    do.call(rbind,.)
  
  return(mut_presence_absence)
}

all_muts <- define_mut_presence_absence(CC_sim_cells)
final_sampled_muts <- all_muts[match(final_sampled_cells$index, CC_sim_cells$index),]

final_sampled_muts <- final_sampled_muts[, colSums(final_sampled_muts != 0) > 0]

xml_file <- args[3]

write('<?xml version="1.0" standalone="yes"?>', xml_file)
write('<beast>', xml_file, append = TRUE)
write(paste0("\t", 	'<taxa id="taxa">'), xml_file, append = TRUE)
for (i in 1:nrow(final_sampled_cells)) {
  
  write(paste0("\t", "\t",	'<taxon id="', 
               final_sampled_cells$index[i],
               "|",
               final_sampled_cells$norm_locx[i],
               "|",
               final_sampled_cells$norm_locy[i],'">'),
        file = xml_file,
        append = TRUE)

  write(paste0('\t', '\t', '\t', '<attr name="x">'), file = xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", "\t", final_sampled_cells$norm_locx[i] ), xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", '</attr>'), file = xml_file, append = TRUE)
  
  write(paste0("\t", "\t", "\t", '<attr name="y">'), xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", "\t", final_sampled_cells$norm_locy[i]), xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", '</attr>'), file = xml_file, append = TRUE)
  
  write(paste0("\t", "\t", "\t", '<attr name="location">'), file = xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", "\t", final_sampled_cells$norm_locx[i], " ", final_sampled_cells$norm_locy[i]), xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", '</attr>'),xml_file, append = TRUE)
  write(paste0("\t", "\t", '</taxon>'), xml_file, append = TRUE)
  

        #   paste(final_sampled_muts[i,], collapse = '')),
        # file = fasta_file,
        # append=TRUE)
  
}

write(paste0("\t", '</taxa>'), xml_file, append = TRUE)
write(paste0("\t", '<alignment id="alignment" dataType="binary">'), xml_file, append = TRUE)

for (i in 1:nrow(final_sampled_cells)) {
  write(paste0("\t", "\t", '<sequence>'), xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", '<taxon idref="',
               final_sampled_cells$index[i],
               "|",
               final_sampled_cells$norm_locx[i],
               "|",
               final_sampled_cells$norm_locy[i],'"/>'), xml_file, append = TRUE)
  write(paste0("\t", "\t", "\t", paste(final_sampled_muts[i,], collapse = '')), xml_file, append = TRUE)
  write(paste0("\t", "\t", '</sequence>'), xml_file, append = TRUE)
}
write(paste0("\t", '</alignment>'), xml_file, append=TRUE)
