####### Phylogenetic analysis of ground truth of simulation
library(tidyverse)
library(phytools)
library(ape)
library(ggtree)
setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")
#function to collect all ancestors of cell
final_sampled_cells <- read.csv(file="outputs/CC_final_sampled_cells.csv")

#first normalize locations so that 0 is the center
normalize_locs <- function(cells){
  center_x <- min(cells$locx) + (max(cells$locx) - min(cells$locx))/2
  center_y <- min(cells$locy) + (max(cells$locy) - min(cells$locy))/2
  norm_cells <- cells %>% 
    mutate("norm_locx" = locx - center_x, "norm_locy" = locy - center_y)
  return(norm_cells)
}
#read in simulation data
CC_sim_cells <- read.csv("outputs/sim_CC_cells.csv") %>% 
  normalize_locs


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
    
  }else {
    
    #make all pariwise combinations of cells without mrca
    all_combs <- as.data.frame(t(combn(orphan_cells, 2)))
    
    #find ancestors of candidate cells 
    common_ancestors <- map2_dbl(all_combs$V1, all_combs$V2, ~find_MRCA(.x, .y, all_cells))
    common_ancestors_data <- all_cells[common_ancestors+1,]
    
    #find youngest ancestor and two children
    youngest_ancestor <- common_ancestors_data$index[common_ancestors_data$birthdate == max(common_ancestors_data$birthdate)]
    
    #in case more than 2 cells have same ancestor
    youngest_ancestor <- youngest_ancestor[1]
    
    #add ancestor to nodes on tree
    node_vec <- c(node_vec, youngest_ancestor)
    
    #keep track of children descended from that ancestor
    children <- unique(unlist(all_combs[which(common_ancestors == youngest_ancestor),]))
    
    #node number
    a <- length(node_vec)
    
    #get time info for ancestor
    t_a <- all_cells[youngest_ancestor+1, "birthdate"]
    
    #make edges for each child
    
    for (child in children) {
      
      c <- which(node_vec == child)
      if (c <= n) { #if leaf
        
        t_c <- all_cells[child + 1, "deathdate"] #use deathdate
        
      } else { #if node
        
        t_c <- all_cells[child + 1, "birthdate"] #use birthday
      }
      
      edge_matrix <- rbind(edge_matrix, c(a,c)) #add edge
      edge_lengths <- c(edge_lengths, t_c - t_a) #add lengths
      orphan_cells <- orphan_cells[orphan_cells != child] #remove from orphan cells
    }
    
    orphan_cells <- c(orphan_cells, youngest_ancestor) #add ancestor as orphan
    
    #recursive function
    return(collect_nodes_for_phylo(orphan_cells, all_cells, n, node_vec = node_vec,
                                   edge_matrix = edge_matrix, edge_lengths = edge_lengths))
    
  }
  
}


#function to take df of cells in simulation and convert to phylo class object

convert_nodes_to_string <- function(orphan_cells, all_cells, subtrees = c()){
  
  all_cells <- all_cells %>% arrange(index)
  
  if (length(subtrees)== 0) {
    subtrees <- as.character(orphan_cells)
  }
  
  if (sum(orphan_cells == 0)/length(orphan_cells) == 1) {
    
    return(paste0(last(subtrees),";"))
    
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
  
    subtrees <- c(subtrees, paste0("(",paste(subtrees[which(orphan_cells %in% children)], collapse = ","),")"))
    
    subtrees <- subtrees[-which(orphan_cells %in% children)] #remove child trees after used
    
    orphan_cells <- c(orphan_cells, youngest_ancestor) #add ancestor as orphan
    orphan_cells <- orphan_cells[-which(orphan_cells %in% children)] #remove children
    
    print(last(subtrees))
    #recursive function
    return(convert_nodes_to_string(orphan_cells, all_cells, subtrees))
    
  }
  
}

#test for clonal simulation results
leaves <- final_sampled_cells$index


tree_text <- convert_nodes_to_string(leaves, CC_sim_cells)

test.tree <- read.tree(text = tree_text)

tree_info <- collect_nodes_for_phylo(leaves, CC_sim_cells)
tr <- list(edge = tree_info$edge_matrix, tip.label = as.character(leaves), Nnode = (length(leaves) - 1),
           edge.length = tree_info$edge_lengths)
class(tr) <- "phylo"

tr <- list(edge = matrix(c(8,1,8,2,7,8,6,7,6,4,5,6), 6, 2), tip.label = c("a","b","c","d"), Nnode = 4L)
class(tr) <- "phylo"

plot(tr)
tr <- list(edge = as.matrix(), tip.label = as.character(leaves), Nnode = (length(leaves) - 1),
           edge.length = tree_info$edge_lengths)
class(tr) <- "phylo"
ggtree(tr)
       
       
 ## read tree from string
text.string<-"(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark, test);"
test.text <- "(((1,2),(3,4)),((5,6),(7,8)));"
test.tree <- read.tree(text = test.text)
vert.tree<-read.tree(text=text.string)
plot(vert.tree)
