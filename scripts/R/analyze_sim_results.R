####### Phylogenetic analysis of ground truth of simulation
library(tidyverse)
library(phytools)
library(ape)
library(ggtree)
library(cowplot)
library(biwavelet)
library(autoimage)
library(gridGraphics)
setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")


#####FUNCTIONS######

#first normalize locations so that 0 is the center
normalize_locs <- function(cells){
  center_x <- min(cells$locx) + (max(cells$locx) - min(cells$locx))/2
  center_y <- min(cells$locy) + (max(cells$locy) - min(cells$locy))/2
  norm_cells <- cells %>% 
    mutate("norm_locx" = locx - center_x, "norm_locy" = locy - center_y)
  return(norm_cells)
}


collect_all_ancestors <- function(index, all_cells, ancestors = c()) {
  if( index == 0) {
    return(ancestors)
  } else {
    new_index <- all_cells[all_cells$index == index,"parent_index"]
    return(collect_all_ancestors(index = new_index, all_cells, ancestors = c(ancestors, new_index)))
  }
}

#function to find MRCA between two cells given indexes

find_MRCA <- function(index_1, index_2, all_cells) {
  ancestors_1 <- collect_all_ancestors(index_1, all_cells)
  ancestors_2 <- collect_all_ancestors(index_2, all_cells)
  mrca <- ancestors_1[min(which(ancestors_1 %in% ancestors_2 ))]
  return(mrca)
}


#function to take df of cells in simulation and convert to phylo class object
#output object crashes plot
# collect_nodes_for_phylo <- function(orphan_cells, all_cells, n, node_vec = c(),
#                                     edge_matrix = matrix(ncol=2), edge_lengths = c()) {
#   
#   all_cells <- all_cells %>% arrange(index)
#   
#   #if just starting record information for leaves
#   if (length(node_vec) == 0) {
#     node_vec <- orphan_cells
#     n <- length(orphan_cells)
#   }
#   if (sum(orphan_cells == 0)/length(orphan_cells) == 1) {
#     
#     trait_data <- data.frame("node" = node_vec, "x_coord" = all_cells[node_vec+1,"norm_locx"],
#                              "y_coord" = all_cells[node_vec+1,"norm_locy"])
#     
#     return(list("trait_data" = trait_data, "edge_matrix" = edge_matrix[-1,], "edge_lengths" = edge_lengths))
#     
#   }else {
#     
#     #make all pariwise combinations of cells without mrca
#     all_combs <- as.data.frame(t(combn(orphan_cells, 2)))
#     
#     #find ancestors of candidate cells 
#     common_ancestors <- map2_dbl(all_combs$V1, all_combs$V2, ~find_MRCA(.x, .y, all_cells))
#     common_ancestors_data <- all_cells[common_ancestors+1,]
#     
#     #find youngest ancestor and two children
#     youngest_ancestor <- common_ancestors_data$index[common_ancestors_data$birthdate == max(common_ancestors_data$birthdate)]
#     
#     #in case more than 2 cells have same ancestor
#     youngest_ancestor <- youngest_ancestor[1]
#     
#     #add ancestor to nodes on tree
#     node_vec <- c(node_vec, youngest_ancestor)
#     
#     #keep track of children descended from that ancestor
#     children <- unique(unlist(all_combs[which(common_ancestors == youngest_ancestor),]))
#     
#     #node number
#     a <- length(node_vec)
#     
#     #get time info for ancestor
#     t_a <- all_cells[youngest_ancestor+1, "birthdate"]
#     
#     #make edges for each child
#     
#     for (child in children) {
#       
#       c <- which(node_vec == child)
#       if (c <= n) { #if leaf
#         
#         t_c <- all_cells[child + 1, "deathdate"] #use deathdate
#         
#       } else { #if node
#         
#         t_c <- all_cells[child + 1, "birthdate"] #use birthday
#       }
#       
#       edge_matrix <- rbind(edge_matrix, c(a,c)) #add edge
#       edge_lengths <- c(edge_lengths, t_c - t_a) #add lengths
#       orphan_cells <- orphan_cells[orphan_cells != child] #remove from orphan cells
#     }
#     
#     orphan_cells <- c(orphan_cells, youngest_ancestor) #add ancestor as orphan
#     
#     #recursive function
#     return(collect_nodes_for_phylo(orphan_cells, all_cells, n, node_vec = node_vec,
#                                    edge_matrix = edge_matrix, edge_lengths = edge_lengths))
#     
#   }
#   
# }


#function to take df of cells in simulation and convert newick string

convert_nodes_to_string <- function(orphan_cells, all_cells, subtrees = c(), node_vec = c(), n = 0){
  
  all_cells <- all_cells %>% arrange(index)
  
  if (length(subtrees)== 0) {
    subtrees <- as.character(orphan_cells)
    node_vec <- orphan_cells
    n <- length(orphan_cells)
  }
  
  if (sum(orphan_cells == 0)/length(orphan_cells) == 1) {
    
    return(list("tree.text"= paste0(last(subtrees),";"), "node_vec" = node_vec))
    
  } else {
    
    #make all pariwise combinations of cells without mrca
    all_combs <- as.data.frame(t(combn(orphan_cells, 2)))
    
    #find ancestors of candidate cells 
    common_ancestors <- map2_dbl(all_combs$V1, all_combs$V2, ~find_MRCA(.x, .y, all_cells))
    common_ancestors_data <- all_cells[common_ancestors+1,]
    
    #find youngest ancestor and two children
    youngest_ancestor <- common_ancestors_data$index[common_ancestors_data$birthdate == max(common_ancestors_data$birthdate)]
    
    #in case more than 2 cells have same ancestor
    youngest_ancestor <- youngest_ancestor[1]
    
    #keep track of children descended from that ancestor
    children <- unique(unlist(all_combs[which(common_ancestors == youngest_ancestor),]))
    t_a <- all_cells[youngest_ancestor +1, "birthdate"]
    
    branch_lengths <- c()
    for (child in children) {
      
      c <- which(node_vec == child)
      if (c <= n) { #if leaf
        
        t_c <- all_cells[child + 1, "deathdate"] #use deathdate
        
      } else { #if node
        
        t_c <- all_cells[child + 1, "birthdate"] #use birthday
      }
      
      branch_lengths <- c(branch_lengths, as.character(round(t_c-t_a,2)))
      
    }
    subtrees <- c(subtrees, paste0("(",paste(paste(subtrees[which(orphan_cells %in% children)],branch_lengths,sep=":"), collapse = ","),")",as.character(youngest_ancestor)))
    
    subtrees <- subtrees[-which(orphan_cells %in% children)] #remove child trees after used
    
    orphan_cells <- c(orphan_cells, youngest_ancestor) #add ancestor as orphan
    node_vec <- c(node_vec, youngest_ancestor)
    orphan_cells <- orphan_cells[-which(orphan_cells %in% children)] #remove children
    
    #recursive function
    return(convert_nodes_to_string(orphan_cells, all_cells, subtrees, node_vec, n = n))
    
  }
  
}


#function to take newick output and convert to class treedata for ggtree

convert_nwk_to_treedata <- function(nwk_list, all_cells) {
  
  #read newick string into tree and convert into correct class
  tree <- as.treedata(read.tree(text = nwk_list$tree.text))
  
  all_cells <- all_cells %>% 
    arrange(index)
  
  #get features for each node
  tree_data <- data.frame("X_coord" = all_cells[nwk_list$node_vec+1,"norm_locx"],
                               "Y_coord" = all_cells[nwk_list$node_vec+1,"norm_locy"],
                               "index" = as.character(nwk_list$node_vec) )
  tree_data <- tree_data[match(c(tree@phylo$tip.label,tree@phylo$node.label), tree_data$index),]
  tree_data$node <- 1:nrow(tree_data)
  tree@treetext <- nwk_list$tree.text
  tree@data <- as_tibble(tree_data)
  
  return(tree)
}



#plot spatial and time phylogenies

plot_space_and_time_phylos <- function(tree,
                               time_file = "phylo.png",
                               comb_file = "phylo_spatial.png",
                               spatial_file = "spatial.png") {

  n <- length(tree@phylo$tip.label)
  m <- tree@phylo$Nnode

  #matrix with trait info for leaves
  X <- get.data(tree) %>% 
    arrange(node) %>% 
    filter(node %in% 1:n) %>% 
    select(X_coord, Y_coord) %>% 
    as.matrix

  rownames(X) <- tree@phylo$tip.label

  #matrix with trait info for nodes
  A <- get.data(tree) %>% 
    arrange(node) %>% 
    filter(node %in% (n+1):(n+m)) %>% 
    select(X_coord, Y_coord) %>% 
    as.matrix

  rownames(A) <- as.character((n+1):(n+m))

  #make continuous color scale 

  color.gradient <- function(x, colors=c("red","yellow","springgreen","royalblue"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }

  #fit color scale to Euclidian distance from center
  col_vec <- c(color.gradient(sqrt(X[,1]^2 + X[,2]^2)), color.gradient(sqrt(A[,1]^2 + A[,2]^2)))
  names(col_vec) <- 1:(n+m)
  
  #make spatial only figure
  plot.new()
  png(filename = spatial_file)
  phylomorphospace(as.phylo(tree),
                   X,
                   A=A,
                   control=list(col.node=col_vec),
                   xlab = "X_coord",
                   ylab = "Y_coord",
                   label='off')
  dev.off()
  #make combined figure
  plot.new()

  reset.par()

  #dev.copy(png,"test.png")
  dev.new()
  par(xpd = NA, # switch off clipping, necessary to always see axis labels
      bg = "transparent", # switch off background to avoid obscuring adjacent plots
      oma = c(2, 2, 0, 0)) # move plot to the right and up
 
  phylomorphospace(as.phylo(tree),
                 X,
                 A=A,
                 control=list(col.node=col_vec),
                 xlab = "X_coord",
                 ylab = "Y_coord",
                 label='off')
  
  p1 <- recordPlot()  

  dev.off()
  p2 <- ggtree(tree) +
    geom_tippoint(aes(fill=as.factor(node)), pch=21, color = "black", size=2.5) +
    geom_point2(aes(fill=as.factor(node)),pch=21, color = "black", size=2.5) +
    scale_fill_manual(values = col_vec)
  
  comb_plot <- plot_grid(p2, p1,
                             labels = 'AUTO')
  
  
  save_plot(time_file, p2)
  save_plot(comb_file, comb_plot, base_asp = 2)
}

####DO ANALYSIS AND PLOTTING####

#function to collect all ancestors of cell
final_sampled_cells <- read.csv(file="outputs/CC_final_sampled_cells.csv")
#read in simulation data
CC_sim_cells <- read.csv("outputs/sim_CC_cells.csv") %>% 
  normalize_locs

#input clonal simulation results

CC_sim_nwk_list <- convert_nodes_to_string(final_sampled_cells$index, CC_sim_cells)
#convert nwk to tree
CC_sim_tree <- convert_nwk_to_treedata(CC_sim_nwk_list, CC_sim_cells)

#plotting
plot_space_and_time_phylos(CC_tree,
                           comb_file = "figures/spatial_sim/sim_phylo_spatial.png",
                           time_file =  "figures/spatial_sim/sim_phylo.png",
                           spatial_file = "figures/spatial_sim/sim_spatial.png")

#save tree R object

saveRDS(CC_sim_tree, file = "outputs/objects/CC_sim_tree.rds")
