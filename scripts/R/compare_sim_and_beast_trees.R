#script to compare mcra spatial location between beast infered and true simulated ancestors

#load libraries
library(tidyverse)
library(phytools)
library(ape)
library(ggtree)
library(phangorn)
library(stringr)
setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")

#read trees 
CC_sim_tree <- readRDS("outputs/objects/CC_sim_tree.rds")
CC_beast_tree <- readRDS("outputs/objects/CC_beast_tree.rds")
beast_tips <- purrr::map_chr(CC_beast_tree@phylo$tip.label, function(s) str_extract(string = s, pattern = "^[0-9]+"))

sim_tips <- CC_sim_tree@phylo$tip.label

clustered_cells <- read_csv("outputs/CC_clustered_cells.csv")

clustered_cells_sampled <- clustered_cells[clustered_cells$index %in% sim_tips,]

sim_clusters <- c(clustered_cells_sampled[match(sim_tips, clustered_cells_sampled$index),]$cluster,
                  rep(NA, CC_sim_tree@phylo$Nnode))

beast_clusters <- c(clustered_cells_sampled[match(beast_tips, clustered_cells_sampled$index),]$cluster,
                  rep(NA, CC_beast_tree@phylo$Nnode))

CC_sim_tree@data <- CC_sim_tree@data %>% 
  arrange(node) %>% 
  add_column("cluster" = sim_clusters)

CC_beast_tree@data <- CC_beast_tree@data %>% 
  arrange(node) %>% 
  add_column("cluster" = beast_clusters)

g <- ggtree(CC_sim_tree) +
geom_tippoint(aes(fill=as.factor(cluster)), pch=21, color = "black", size=3)

ggsave("figures/spatial_sim/clustered_phylo.png",g)

g <- ggtree(CC_beast_tree) +
  geom_tippoint(aes(fill=as.factor(cluster)), pch=21, color = "black", size=2.5)
ggsave("figures/spatial_sim/clustered_phylo_beast.png",g)






#get corresponding node for each list of children in beast tree
beast_tree <- CC_beast_tree

beast_node_vec <- c()
for (i in 1:length(node_child_list)) {

  child_vec <- node_child_list[[i]]
  strings <- beast_tree@phylo$tip.label
  leaves <- purrr::map_chr(strings, function(s) str_extract(string = s, pattern = "^[0-9]+"))
  #patterns <- as.character(child_vec)
  child_nodes <- which(leaves %in% as.character(child_vec))
  #child_nodes <- purrr::map(patterns, function(p) str_which(p,strings))

  parent_node <- findMRCA(beast_tree@phylo, tips=child_nodes)

  
  beast_node_vec <- c(beast_node_vec, parent_node)
}

beast_x_coord_vec <- beast_tree@data[match(beast_node_vec, beast_tree@data$node),]$location1
beast_y_coord_vec <- beast_tree@data[match(beast_node_vec, beast_tree@data$node),]$location2

beast_posterior_vec <- beast_tree@data[match(beast_node_vec, beast_tree@data$node),]$posterior


sim_x_coord_vec <- tree@data[match(all_nodes, tree@data$node),]$X_coord
sim_y_coord_vec <- tree@data[match(all_nodes, tree@data$node),]$Y_coord

n <- length(tree@phylo$tip.label)
m <- n + tree@phylo$Nnode

node_conversion_df <- data.frame("sim_node" = all_nodes,
                                 "beast_node" = beast_node_vec,
                                 "sim_x_coord" =sim_x_coord_vec,
                                 "sim_y_coord" =sim_y_coord_vec,
                                 "beast_x_coord" = beast_x_coord_vec,
                                 "beast_y_coord" = beast_y_coord_vec,
                                 "beast_posterior" = beast_posterior_vec) %>% 
  mutate("euclid_dist" = sqrt((sim_x_coord - beast_x_coord)^2 + (sim_y_coord - beast_y_coord)^2),
         "sim_height" = node.height(tree@phylo)[(n+1):m]) 

max_distance = sqrt((max(node_conversion_df$sim_x_coord) - min(node_conversion_df$sim_x_coord))^2 +
  (max(node_conversion_df$sim_y_coord) - min(node_conversion_df$sim_y_coord))^2)

sim_beast_comparison_df <- node_conversion_df %>% 
  mutate("percent_euclid_dist" = euclid_dist/max_distance)


g <- ggplot(sim_beast_comparison_df, aes(x = sim_height, y = percent_euclid_dist, color = beast_posterior)) +
  geom_point(size = 3) + theme_bw() + xlab("node height") + ylab("Fraction distance") + ggtitle("Location error")

ggsave("figures/spatial_sim/node_distance_vs_height.png", g)

 




# 
# g <- ggtree(beast_tree, aes(color = as.factor(alleles_fasta_meta))) + geom_tiplab(size = 2, color = "black")
# g <- ggtree(beast_tree, aes(color = as.factor(index))) + geom_tiplab(size = 2, color = "black")
# ggsave("figures/spatial_sim/beast_labeled_tips.png", g)
# g <- ggtree(tree, aes(color = as.factor(index))) + geom_tiplab(size = 2, color = "black")
# ggsave("figures/spatial_sim/sim_labeled_tips.png", g)


