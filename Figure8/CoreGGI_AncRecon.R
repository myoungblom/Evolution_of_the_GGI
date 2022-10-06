require(ape)
require(adegenet)
require(viridis)
require(factoextra)
require(phytools)
require(dplyr)
require(ggtree)
require(ggimage)
require(ggplot2)
require(treeio)
require(cowplot)
require(scales)
require(splitstackshape)
require(reshape2)
require(phangorn)
library(tidyverse)

###############
##### PCA #####
###############

# read core ggi alignment fasta
ggi <- fasta2DNAbin("CoreGGI_aln.fasta")

# convert to genid object
ggi_obj <- DNAbin2genind(ggi)

# replace missing data
X <- scaleGen(ggi_obj, NA.method="zero")

# pca analysis
pca1 <- dudi.pca(X, cent=F, scale=F, scannf=F, nf=3)

# assigning PCA groups based on clustering
pca_groups <- pca1$l1
pca_groups$isolate <- str_split_fixed(rownames(pca1$l1),"_",2)[,1]
pca_groups <- pca_groups %>% mutate(group=case_when(RS1 < 0 & RS2 > 0 ~ "I",
                                                 RS1 > 0 & RS2 > 0 ~ "II",
                                                 RS1 > 0 & RS2 < 0 ~ "IV",
                                                 RS1 < 0 & RS2 < 0 ~ "III"))
pca_groups$group <- as.factor(pca_groups$group)

# moving some group VI's to group I bc they cluster better
move <- c(107,15,97)
pca_groups$group[move] <- "I"


custom_colors <- viridis(n=4)

# plot with ggplot
PCAplot <- ggplot()+geom_point(data=pca_groups,aes(x=RS1,y=RS2,color=group),size=4, alpha=0.7)+
  theme_minimal()+xlab("PC1: 16.2%")+ylab("PC2: 9.3%")+ 
  scale_color_manual(values=custom_colors,name="Group",labels=levels(pca_groups$group),breaks=levels(pca_groups$group))+
  geom_hline(yintercept=0,lty=2,color="grey") + geom_vline(xintercept=0,lty=2,color="grey")

PCAplot

############################################
##### PLOT GGI GROUPS ON CORE GGI TREE #####
############################################

ggi_tree <- read.tree("GGIcore_RAxML.newick")
ggi_data <- pca_groups[,c(4,5)]

#Order sources info dataframe by tree tip label
ggi_data <- ggi_data[match(ggi_tree$tip.label,ggi_data$isolate),]

#save tip labels before adding the source
OG_tipLabel <- ggi_tree$tip.label

#Add meta info to tree tip labels
ggi_tree$tip.label <- paste(ggi_data$group, "_",ggi_tree$tip.label, sep = "")

#group taxa by meta info
groupInfo <- split(ggi_tree$tip.label,gsub("_\\w+","",ggi_tree$tip.label))

#group branches/tips by groupInfo
ggi_tree <- groupOTU(ggi_tree,groupInfo)

#replace tip labels with originals
ggi_tree$tip.label <- OG_tipLabel

#Get data to color branch by group
dat <- fortify(ggi_tree)
index.pos <- dat$isTip == TRUE
index.neg <- dat$isTip == FALSE
dat$ggi[index.neg] <- "aaaaINTERNAL"
dat$ggi[index.pos] <- ggi_data$group
dat$ggi <- as.factor(dat$ggi)

#colors
custom_colors = c(viridis(n=4),"black")

#plot tree
r <- ggtree(ggi_tree,layout = "circular", size = 0.75)
r

#add tips
s <-  r %<+% dat + aes(color=(ggi)) + 
  geom_tippoint(aes(color=ggi),size=2) + 
  scale_color_manual(name = "Group", values=c(custom_colors), labels = levels(ggi_data$group), breaks=c("1","2","3","4"))+
  theme(legend.position="none")  
s

####################################
##### ANCESTRAL RECONSTRUCTION #####
####################################

# reading tree and metadata
tree <- read.tree("GCcore_RAxML.newick")
group_data <- pca_groups[,c(4,5)]
ggineg <- tree$tip.label[!(tree$tip.label %in% group_data$isolate)]
group_data <- rbind(group_data,
              data.frame(isolate=ggineg,group="Absent"))

# putting GGI groups in order of tree tip labels
tip.label.order <- rev(group_data$isolate[match(tree$tip.label,group_data$isolate)])
anc_data <- group_data %>% mutate(isolate = factor(isolate,levels=tip.label.order)) %>% arrange(isolate)
ggi.groups <- anc_data$group
names(ggi.groups) <- anc_data$isolate
ggi.groups <- as.factor(ggi.groups)

# running ancestral reconstruction
ERrecon <- ace(ggi.groups, tree, type="discrete", model="ER")

####################################
#### PLOTTING ANC RECON ON TREE ####
####################################

# plot the tree first

# ggi groups
tree_data <- group_data

# reorder ggi groups to match tree
tree_data <- tree_data[match(tree$tip.label,tree_data$isolate),]

# save tip labels before adding the source
OG_tipLabel <- tree$tip.label

# add meta info to tree tip labels
tree$tip.label <- paste(tree_data$group, "_",tree$tip.label, sep = "")

# group taxa by meta info
groupInfo <- split(tree$tip.label,gsub("_\\w+","",tree$tip.label))

# group branches/tips by groupInfo
tree <- groupOTU(tree,groupInfo)

# replace tip labels with originals
tree$tip.label <- OG_tipLabel

# get data to color branch by group
dat <- fortify(tree)
index.pos <- dat$isTip == TRUE
index.neg <- dat$isTip == FALSE
dat$ggi[index.neg] <- "aaaaINTERNAL"
dat$ggi[index.pos] <- tree_data$group
dat$ggi <- as.factor(dat$ggi)

#colors
custom_colors <- c("black",viridis(n=4),"black")

#plot tree
r <- ggtree(tree,layout = "rectangular", size = 0.75)
r

#add tips
s <-  r %<+% dat + aes(color=ggi) +
  geom_tippoint(aes(color=ggi),size=1.5) +
  scale_color_manual(name="GGI Group",values=c(custom_colors), 
                     labels = c("Absent","I","II","III","IV"),
                     breaks=c(5,1,2,3,4))+
  theme(legend.position = c(0.1,0.9), text = element_text(size=10, face = "bold"))

s

# ancestral reconstruction data for plotting at each node
pie.data <- as.data.frame(ERrecon$lik.anc)
pie.data$node <- 1:tree$Nnode+Ntip(tree)
pies <- nodepie(pie.data,cols=c(1:5), alpha=0.8)
pies <- lapply(pies, function(g) g+scale_fill_manual(values=c("black",viridis(n=4)),name="GGI Group",
                                                     breaks=c("Absent","I","II","III","IV"), labels=c("Absent","I","II","III","IV")))

# add pie charts to node of tree
pie.tree <- s + geom_inset(pies, width=0.065, height=0.065)
pie.tree
