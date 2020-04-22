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

#first normalize locations so that 0 is the center
normalize_locs <- function(cells){
  center_x <- min(cells$locx) + (max(cells$locx) - min(cells$locx))/2
  center_y <- min(cells$locy) + (max(cells$locy) - min(cells$locy))/2
  norm_cells <- cells %>% 
    mutate("norm_locx" = locx - center_x, "norm_locy" = locy - center_y)
  return(norm_cells)
}

norm_cells_df <- normalize_locs(cells_df)
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

#find MRCA of sampled cells

collect_all_ancestors <- function(index, all_cells, ancestors = c()) {
  if( index == 0) {
    return(ancestors)
  } else {
    new_index <- all_cells[all_cells$index == index,"parent_index"]
    return(collect_all_ancestors(index = new_index, all_cells, ancestors = c(ancestors, new_index)))
  }
}

test <- collect_all_ancestors(index = 42, cells_df)
test_2 <- collect_all_ancestors(index = 153, cells_df)
test[min(which(test %in% test_2))]

find_MRCA <- function(index_1, index_2, all_cells) {
  ancestors_1 <- collect_all_ancestors(index_1, all_cells)
  ancestors_2 <- collect_all_ancestors(index_2, all_cells)
  mrca <- ancestors_1[min(which(ancestors_1 %in% ancestors_2 ))]
  return(mrca)
}

collect_nodes_for_phylo <- function(orphan_cells, all_cells, n, node_vec = c(),
                                    edge_matrix = matrix(ncol=2), edge_lengths = c()) {
  all_cells <- all_cells %>% arrange(index)
  #if just starting record information for leaves
  if (length(node_vec) == 0) {
    node_vec <- orphan_cells
    n <- length(orphan_cells)
  }
  if (sum(orphan_cells == 0)/length(orphan_cells) == 1) {

    trait_data <- data.frame("node" = node_vec, "x_coord" = all_cells[node_vec+1,"norm_locx"],
                             "y_coord" = all_cells[node_vec+1,"norm_locy"])
  
    return(list("trait_data" = trait_data, "edge_matrix" = edge_matrix[-1,], "edge_lengths" = edge_lengths))
  } else {
    #make all pariwise combinations of cells without mrca
    all_combs <- as.data.frame(t(combn(orphan_cells, 2)))
    
    #find ancestors of candidate cells 
    common_ancestors <- map2_dbl(all_combs$V1, all_combs$V2, ~find_MRCA(.x, .y, all_cells))
    common_ancestors_data <- all_cells[common_ancestors+1,]
    
    #find youngest ancestor and two children
    youngest_ancestor <- common_ancestors_data$index[common_ancestors_data$birthdate == max(common_ancestors_data$birthdate)]
    
    if (length(youngest_ancestor) == 1) {

    #add new node 
    node_vec <- c(node_vec, youngest_ancestor)

    child_1 <- all_combs[which(common_ancestors == youngest_ancestor),1]
    child_2 <- all_combs[which(common_ancestors == youngest_ancestor),2]
    c1 <- which(node_vec == child_1)
    c2 <- which(node_vec == child_2)
    a <- length(node_vec)

    #update edge matrrix
    edge_matrix <- rbind(edge_matrix, c(a,c1), c(a,c2))
    
    if (c1 <= n) {
      t_c1 <- all_cells[child_1+1, "deathdate"]
    } else {
      t_c1 <- all_cells[child_1+1, "birthdate"]
    }
    
    if (c2 <= n) {
      t_c2 <- all_cells[child_2+1, "deathdate"]
    } else {
      t_c2 <- all_cells[child_2+1, "birthdate"]
    }
    
    t_a <- all_cells[youngest_ancestor+1, "birthdate"]
    
    edge_lengths <- c(edge_lengths, (t_c1 - t_a), t_c2 - t_a)
    
    orphan_cells <- orphan_cells[orphan_cells != child_1]
    orphan_cells <- orphan_cells[orphan_cells != child_2]
    orphan_cells <- c(orphan_cells, youngest_ancestor)
  
    } else {
    youngest_ancestor <- youngest_ancestor[1]
    node_vec <- c(node_vec, youngest_ancestor)
    children <- unique(unlist(all_combs[which(common_ancestors == youngest_ancestor),]))
    
    a <- length(node_vec)
    t_a <- all_cells[youngest_ancestor+1, "birthdate"]


    for (child in children) {

      c <- which(node_vec == child)
      if (c <= n) {
        t_c <- all_cells[child_1+1, "deathdate"]
      } else {
        t_c <- all_cells[child_1+1, "birthdate"]
      }
      edge_matrix <- rbind(edge_matrix, c(a,c))
      edge_lengths <- c(edge_lengths, t_c - t_a)
      orphan_cells <- orphan_cells[orphan_cells != child]
      }
    orphan_cells <- c(orphan_cells, youngest_ancestor)
    }
    return(collect_nodes_for_phylo(orphan_cells, all_cells, n, node_vec = node_vec,
                                        edge_matrix = edge_matrix, edge_lengths = edge_lengths))
    
    
  }
  
}
leaves <- final_sampled_cells$index
tree_info <- collect_nodes_for_phylo(leaves, norm_cells_df)
tr <- list(edge = tree_info$edge_matrix, tip.label = leaves, Nnode = (length(leaves) - 1),
           edge.length = tree_info$edge_lengths)
class(tr) <- "phylo"

plot(tr)
ggtree(tr)
#for tree
#nodes list
find_MRCA(42,153,cells_df)
