#plot Beast output

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("treeio")

setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")
beast_file <- system.file("outputs/CC_output.tree", package="ggtree")

CC_sim_tree <- read.beast("outputs/CC_mcc.tree")
  
ggtree(CC_sim_tree, aes(color = X_coord))

ggtree(CC_sim_tree, 
       yscale = "Y_coord")

ggplot(CC_sim_tree, aes(x = X_coord, y = Y_coord)) + geom_point()

str(CC_sim_tree)
str(as.phylo(CC_sim_tree))
tree_data <- get.data(CC_sim_tree)

ggtree(CC_sim_tree)
#+ 
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3) + 
  geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9, 
                 x=branch), vjust=0) 
  
  file <- system.file("extdata/BEAST", "beast_mcc.tree", 
                      package="treeio")
  beast <- read.beast(file)
  ggtree(beast) + 
    geom_tiplab(align=TRUE, linetype='dashed', linesize=.3) + 
    geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + 
    geom_text2(aes(label=round(as.numeric(posterior), 2), 
                   subset=as.numeric(posterior)> 0.9, 
                   x=branch), vjust=0) 

   