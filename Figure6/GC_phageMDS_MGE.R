require(gplots)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
require(dplyr)
require(tidyr)
require(pheatmap)
require(stringr)

# phage analyses

ggipheno <- read.csv("ggi_scoary_traits.csv")
ggineg <- ggipheno$X[ggipheno$ggi == "0"]
ggipos <- ggipheno$X[ggipheno$ggi == "1"]


mash <- read.delim("mash_results.txt",header=F)
mash <- mash[c(1,2,3)]
colnames(mash) <- c("binA","binB","distance")

matrix.results <- acast(mash, binA~binB, value.var="distance")

fit <- cmdscale(matrix.results, eig=T, k=2)
mds <- as.data.frame(fit$points)
colnames(mds) <- c("coord1","coord2")
mds$file <- rownames(mds)
mds$isolate <- str_split_fixed(mds$file,"[.]",2)[,1]
mds <- points %>% mutate(ggi=case_when(isolate %in% ggineg ~ "neg",
                                          isolate %in% ggipos ~ "pos"))
mds$ggi <- as.factor(mds$ggi)

plot.ggi <- ggplot(data=mds, aes(x=coord1, y=coord2, color=ggi)) + geom_jitter(size=3, alpha=0.6, width=.15,height=.15)+
  xlim(-1,1) + ylim(-1,1) + scale_color_manual(name=NULL,values=c("black","#0066FF"),breaks=c("neg","pos"),labels=c("GGI-","GGI+"))+
  theme_minimal() + xlab("Coordinate 1") + ylab("Coordinate 2")
plot.ggi

# mge

data <- read.delim("mge_summary.txt")

plot <- ggplot(data = data, aes(x=ggi, y=total)) + geom_boxplot(aes(fill=ggi)) + theme_minimal()+
  scale_x_discrete(labels = c("GGI-", "GGI+")) + xlab(NULL) + ylab("Number MGE per strain")+
  scale_fill_manual(name=NULL,values=alpha(c("black","#0066FF"),0.75),breaks=c("neg","pos"),labels=c("GGI-","GGI+"))
plot

stat.test <- compare_means(total ~ ggi, data=data, method="wilcox.test", p.adjust.method = "BH")
stat.test <- stat.test %>%
  mutate(y.position=20)

stat.test$p.adj <- formatC(stat.test$p.adj, format="e", digits = 1)

plot.stats <- plot + stat_pvalue_manual(stat.test, label="p.adj", label.size = 4.5)
plot.stats
