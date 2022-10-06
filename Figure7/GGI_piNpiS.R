require(ggplot2)
require(dplyr)

setwd("~/Desktop/2022.10.05_GGI_Revisions/GithubRepo/Figure7/")

#### piN/piS per gene (boxplot for each gene) ####
ggi_genes <- read.delim("ggi_gene_order.txt",header=T)
gene_order <- ggi_genes$gene

pinpis <- read.delim("GGI_piNpiS.txt",sep="\t",header=T,na.strings = "None")
pinpis <- pinpis %>% mutate(gene=factor(gene,levels=gene_order)) %>% arrange(gene)
pinpis$gene <- as.factor(pinpis$gene)

pinpis[is.na(pinpis)] <- 0

g <- ggplot() + geom_boxplot(data=pinpis,aes(x=gene,y=PiNPiS,fill=core)) + theme_minimal() +
  theme(axis.text.x = element_text(angle=45)) + xlab("\nGGI gene") +
  scale_fill_manual(name="Gene Type",values=c("darkblue","dodgerblue","grey","white"), 
                    breaks=c("core","soft core","shell","cloud"), labels=c("Core","Soft Core","Shell","Cloud"))+
  ylab(expression(pi~"N/"~pi~"S")) + theme(axis.text.x = element_text(face="italic"))

g
