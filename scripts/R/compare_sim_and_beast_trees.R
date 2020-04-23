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


#find equivalent nodes
#nodes can be defined by unique descendents 
#can't assume that beast correctly identified all nodes
#so start with true simulated node and find common mrca with all of those descendents. 
tree <- CC_sim_tree

#get parameters of leaves and nodes
n <- length(tree@phylo$tip.label)
m <- n + tree@phylo$Nnode

#make vector of all nodes
all_nodes <- (n+1):m

#get list of children for each node
node_child_list <- list()

for (n in all_nodes) {
  children <- unlist(Descendants(tree@phylo, n ,"tips"))
  children_indexes <- tree@data[children, ]$index
  node_child_list[[as.character(n)]] <- children_indexes
  
}

#get corresponding node for each list of children in beast tree
beast_tree <- CC_beast_tree

beast_node_vec <- c()
for (i in 1:length(node_child_list)) {

  child_vec <- node_child_list[[i]]
  strings <- beast_tree@phylo$tip.label
  patterns <- as.character(child_vec)

  child_nodes <- purrr::map_int(patterns, function(p) str_which(p,strings))
  
  parent_node <- findMRCA(beast_tree@phylo, tips=child_nodes )

  
  beast_node_vec <- c(beast_node_vec, parent_node)
}

beast_x_coord_vec <- beast_tree@data[match(beast_node_vec, beast_tree@data$node),]$X_coord
beast_y_coord_vec <- beast_tree@data[match(beast_node_vec, beast_tree@data$node),]$Y_coord

beast_posterior_vec <- beast_tree@data[match(beast_node_vec, beast_tree@data$node),]$posterior


sim_x_coord_vec <- tree@data[match(all_nodes, tree@data$node),]$X_coord
sim_y_coord_vec <- tree@data[match(all_nodes, tree@data$node),]$Y_coord


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
  geom_point(size = 3) + theme_bw() + xlab("node height") + ylab("Fraction distance") + ggtitle("Distance between nodes")

ggsave("figures/spatial_sim/node_distance_vs_height.png", g)

 
ggplot(sim_beast_comparison_df, aes(x = beast_posterior, y = percent_euclid_dist)) +
  geom_point() + theme_bw()


g <- ggtree(beast_tree, aes(color = X_coord)) +
  geom_text2(aes(label=round(as.numeric(posterior), 2)), color = "black")

ggsave("figures/spatial_sim/beast_phylo_w_posterior_labels.png", g)


