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
beast_file <- system.file("outputs/CC_output.tree", package="ggtree")

CC_sim_tree <- read.beast("outputs/CC_mcc.tree")


#plot spatial and time phylogenies

n <- length(CC_sim_tree@phylo$tip.label)
m <- CC_sim_tree@phylo$Nnode

#matrix with trait info for leaves
X <- get.data(CC_sim_tree) %>% 
  arrange(node) %>% 
  filter(node %in% 1:n) %>% 
  select(X_coord, Y_coord) %>% 
  as.matrix

rownames(X) <- CC_sim_tree@phylo$tip.label

#matrix with trait info for nodes
A <- get.data(CC_sim_tree) %>% 
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

#make combined figure
plot.new()

reset.par()
#dev.new()
dev.copy(png,"figures/spatial_sim/beast_spatial.png")
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(2, 2, 0, 0)) # move plot to the right and up
phylomorphospace(as.phylo(CC_sim_tree),
                 X,
                 A=A,
                 control=list(col.node=col_vec),
                 xlab = "X_coord",
                 ylab = "Y_coord",
                 label='off')
p1 <- recordPlot()  
#p3 <- plot_grid(p1)
#fancyTree(as.phylo(CC_sim_tree),type="scattergram",X=X,A=A,label = "off")
dev.off()
p2 <- ggtree(CC_sim_tree) +
  geom_tippoint(aes(fill=as.factor(node)), pch=21, color = "black", size=2.5) +
  geom_point2(aes(fill=as.factor(node)),pch=21, color = "black", size=2.5) +
  scale_fill_manual(values = col_vec)

beast_results <- plot_grid(p2, p1,
          labels = 'AUTO')


save_plot("figures/spatial_sim/beast_phylo.png", p2)
save_plot("figures/spatial_sim/beast_phylo_spatial.png", beast_results)


tree_data <- get.data(CC_sim_tree)

