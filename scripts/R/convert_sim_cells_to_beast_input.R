#script to sample cells tumor simulations output and convert to beast output

#load packages
library(tidyverse)
library(cluster)

#intializatin
set.seed(192837)
setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")

#read in simulation data
CC_sim_cells <- read.csv("outputs/sim_CC_cells.csv")
CSC_sim_cells <- read.csv("outputs/sim_CSC_cells.csv")


#temp global var for developing functions
cells_df <- CC_sim_cells

#calculate parameters of simulations
# calc_sim_stats <- function(cells_df) {
#   stats_df <- data.frame("avg_mutation_rate" = mean(cells_df$mutation_rate),
#                          "avg_proliferation_rate" = mean(cells_df$proliferation_rate),
#                          "avg_death_rate" = mean(cells_df$alpha))
#   return(stats_df)
# }
# 
# CC_stats_df <- calc_sim_stats(CC_sim_cells)
# CC_stats_df <- calc_sim_stats(CC_sim_cells)


#filter surviving cells

filter_alive_cells <- function(cells_df) {
  endpoint <- max(cells_df$deathdate)
  
  alive_cells <- cells_df %>% 
    filter(deathdate == endpoint)
  
  return(alive_cells)
}

#sample surviving cells

sample_alive_cells <- function(alive_cells, n){
  
  sampled_indexes <- alive_cells[sample(1:nrow(alive_cells),
                            size = min(n, nrow(alive_cells))), "index"]
  sampled_alive_cells <- alive_cells %>% 
    mutate("sampled" = index %in% sampled_indexes)
  
  return(sampled_alive_cells)
}

#filter and sample alive cells

CC_sampled_cells <- CC_sim_cells %>% 
  filter_alive_cells %>% 
  sample_alive_cells(., n = 20)

CSC_sampled_cells <- CSC_sim_cells %>% 
  filter_alive_cells %>% 
  sample_alive_cells(., n = 20)

#plotting sampled and simulated cells
plot_sampled_cells <- function(sampled_cells, model = "") {
  g <- ggplot(sampled_cells, aes(x = locx,
                                 y = locy,
                                 color = sampled)) +
    geom_point(size = 2.5) +
    theme_bw() +
    xlab("X") +
    ylab("Y") +
    scale_color_manual(values = c("FALSE" = "gray36", "TRUE" = "brown3")) +
    theme(legend.position = "none") +
    coord_fixed() +
    ggtitle(model)
  g
}


plot_sampled_cells(CC_sampled_cells, "Clonal model")
plot_sampled_cells(CSC_sampled_cells, "CSC model")

#convert sampled alive cells into mutation matrix


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
  rowise_muts <- map(sampled_cells$mutations, get_row_muts)
  
  mut_presence_absence <- map(rowise_muts, ~compare_muts(.,all_muts)) %>% 
    do.call(rbind,.)
  
  return(mut_presence_absence)
}

cluster_cells <- function(sampled_cells, k = 5) {
  
  pa <- as.data.frame(define_mut_presence_absence(sampled_cells))

  clusters <- kmeans(pa, centers = k)$cluster
  clustered_cells <- sampled_cells %>% 
    add_column("cluster" = clusters)
  
  return(clustered_cells)
}

set.seed(2018)
clustered_cells <- cluster_cells(sampled_cells, k = 7)


plot_clustered_cells <- function(clustered_cells, model = "") {
  g <- ggplot(clustered_cells, aes(x = locx,
                                 y = locy,
                                 color = 
                                                as.factor(cluster)
                                                )
                                 ) +
    geom_point(size = 2) +
    theme_bw() +
    xlab("X") +
    ylab("Y") +
    #scale_color_manual(values = c("FALSE" = "gray36", "TRUE" = "brown3")) +
    theme(legend.position = "none") +
    coord_fixed() +
    ggtitle(model)
  g
}

plot_clustered_cells(clustered_cells)


#first normalize locations so that 0 is the center
normalize_locs <- function(cells){
  center_x <- min(cells$locx) + (max(cells$locx) - min(cells$locx))/2
  center_y <- min(cells$locy) + (max(cells$locy) - min(cells$locy))/2
  norm_cells <- cells %>% 
    mutate("norm_locx" = locx - center_x, "norm_locy" = locy - center_y)
  return(norm_cells)
}


normalized_cell <- normalize_locs(clustered_cells)


#generate fasta file for beast
#>index|clone|xloc|yloc
#1100---000000-00000-000-0000-0-0---
final_sampled_cells <- normalized_cell %>% 
  filter(sampled == TRUE)

final_sampled_muts <- define_mut_presence_absence(sampled_cells)[which(final_sampled_cells$sampled == TRUE),]
final_sampled_muts <- final_sampled_muts[, colSums(final_sampled_muts != 0) > 0]

#final_sampled_muts[final_sampled_muts==0]<-"A"
#final_sampled_muts[final_sampled_muts==1]<-"G"

fasta_file <- "outputs/CC_sim.fa"

write(c(paste0(">",
               final_sampled_cells$index[1],
               "|",
               final_sampled_cells$norm_locx[1],
               "|",
               final_sampled_cells$norm_locy[1]),
        paste(final_sampled_muts[i,], collapse = '')),
      file = fasta_file,
      append=FALSE)

for (i in 2:nrow(final_sampled_cells)) {
  
  write(c(paste0(">",
                      final_sampled_cells$index[i],
                      "|",
                      final_sampled_cells$norm_locx[i],
                      "|",
                      final_sampled_cells$norm_locy[i]),
               paste(final_sampled_muts[i,], collapse = '')),
             file = fasta_file,
             append=TRUE)

}



