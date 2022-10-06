require(ggplot2)
require(scales)
require(dplyr)

# maximum of null distribution = significance cutoff
fst_sig <- read.delim("ggi_nullWcFst.txt",header=T)$Max[1]

# read fst data for all SNPs in core genome
fst <- read.delim("ggi_wcFst.txt", header = FALSE)
colnames(fst) <- c("seqid","position","targetFreq","backgroundFreq","wcFst")

# read homoplasy data from all core genome homoplasies (ie mutations that occured >= 2 on the phylogeny)
homo <- read.delim("homoplasy/all_gc_homoplasies.txt")

# make separate df for homoplasies in or above the 95th percentile
homo_outliers <- subset(homo, homo$mutationMultiplicity > quantile(homo$mutationMultiplicity,0.95))

# dataframe combining fst values for homoplasy outliers
# arranging dataframe so homoplasies are last (so they are plotted last and stand out on top of fst outliers)
fst <- fst %>% mutate(sig = case_when(position %in% homo_outliers$position ~ "zhomo",
                                      wcFst > fst_sig ~ "sig",
                                      wcFst < fst_sig ~ "NS")) %>%
  arrange(sig)

# plot colors
fst_colors <- c("lightgrey","black","red")

# manhattan plot with significant fst ouliers colored in black and homoplasies colored in red
fstPlot <- ggplot(fst, aes(x=position, y = wcFst)) + geom_point(aes(color=sig), alpha=0.7, size=2)+
  xlab("\nGenome Position") + theme_bw()+
  geom_hline(yintercept = fst_sig, linetype = "dashed", color = "red", size = 1)+
  theme(plot.title = element_text(size = 18, face = "bold", color = "grey27", family = "sans"), 
        axis.text = element_text(color = "grey27", face = "bold", family = "sans", size = 12), 
        axis.title = element_text(color = "grey27", face = "bold", size = 16, family = "sans"),
        legend.position = "bottom", legend.text = element_text(size = 12, color = "grey27", family = "sans"))+
  scale_color_manual(name = NULL, values = fst_colors, labels = c("NS Fst", "Sig Fst", "95th percentile Homoplasy"),
                     breaks = c("NS", "sig","zhomo"))+
  scale_x_continuous(labels=comma)
fstPlot
