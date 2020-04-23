#plot Beast output

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("treeio")
library("phytools")
library(cowplot)
library(biwavelet)
library(autoimage)
library(gridGraphics)


setwd("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")


CC_beast_tree <- read.beast("outputs/CC_sim_2.mcc.tree")


#plot spatial and time phylogenies

plot_space_and_time_phylos <- function(tree,
                                       time_file = "temporal.png",
                                       comb_file = "spatial_temporal.png",
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


plot_space_and_time_phylos(CC_beast_tree,
                           comb_file = "figures/spatial_sim/beast_phylo_spatial_2.png",
                           time_file =  "figures/spatial_sim/beast_phylo_2.png",
                           spatial_file ="figures/spatial_sim/beast_spatial_2.png")

saveRDS(CC_beast_tree, file = "outputs/objects/CC_beast_tree_2.rds")



