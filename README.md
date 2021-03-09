# MyRplots
A work flow of creating plots for genomic data

### Creating similarity distance triangle heat map
library(reshape2)
library(ggplot2)

Betapartiti <- read.table("Betapartitivirus3_0506.txt", fill=TRUE, col.names=paste("V",1:35), row.names=1)
colnames(Betapartiti) <- rownames(Betapartiti)
Betapartiti <- apply(Betapartiti,2, function(x){(1/(1+x))*100})
Betapartiti <- melt(Betapartiti, na.rm=TRUE)

pdf("Betapartiti.pdf")
ggplot(data = Betapartiti, aes(Var1, Var2, fill = value))+
geom_tile(color = "white")+
scale_fill_gradient2(low = "blue", high = "red", mid = "yellow",
midpoint = 70, limit = c(40,100), space = "Lab",
name="Percentage\nSimilarity") +
geom_text(aes(label = round(value, 0)), size=1)+
theme_minimal()+
theme(panel.grid.major.x = element_line(size=.1, color="white") , panel.grid.major.y = element_line(size=.1, color="white"))+
theme(axis.text.x = element_text(angle = 90, vjust = 1,
size = 6, hjust = 1))+
theme(axis.text.y = element_text(angle = 0, vjust = 1,
size = 6, hjust = 1))+
coord_fixed()+
labs(x="Virus Names", y="Virus Names")
dev.off()
![Betapartiti](https://user-images.githubusercontent.com/5850834/110440473-5e0ebe00-80c1-11eb-829a-0c36e91ab09a.jpg)

### Creating OTU taxonomy tree plot
mafft --retree 2 --maxiterate 1000 --thread 10 --leavegappyregion otu.fasta > otu.aln
FastTree -pseudo otu.aln > psudu_otu.tree
FastTree -intree psudu_otu.tree -pseudo otu.aln > otu.tree
# in R 
library(ape)
library(ggtree)
info <- read.table("info.txt", header=TRUE, row.names=1)
tree<-read.tree("otu.tree")
p <- ggtree(tree, layout='circular') %<+% info +
geom_tippoint(aes(color=Division), size=1)+
scale_color_manual(values=c("green", "purple","red"))+
geom_tiplab(aes(color=Division, angle=angle), align=T, linetype="dotted", size=2, hjust=1, linesize = 0.2, offset=0.75)+
theme(legend.position = "top")
![otu](https://user-images.githubusercontent.com/5850834/110444773-13dc0b80-80c6-11eb-9315-3ec52b6177b9.jpg)

