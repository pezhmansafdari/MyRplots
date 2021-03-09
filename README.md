# MyRplots
A work flow of creating plots for genomic data

### Creating similarity distance triangle heat map in R
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

### creating animated 3D PCA plot with PCA data extracted from DeSeq2 in R
library(rgl)
pcaData <- read.table("pcaData.txt", sep="\t", header=TRUE, row.names=1)
colvec <- c("#FF0000", "#FF0000", "#FF0000", "#0000FF", "#0000FF", "#0000FF", "#000000", "#000000", "#808080", "#006600", "#006600", "#006600", "#FF8000", "#FF8000", "#FF8000", "#FF99FF", "#FF99FF", "#FF99FF", "#33FFFF", "#33FFFF", "#33FFFF", "#808080", "#000000", "#808080", "#66FF66", "#66FF66","#66FF66", "#FFFF00", "#FFFF00","#FFFF00")
legendCol <- c("#FF0000", "#0000FF","#000000","#006600","#FF8000","#FF99FF","#33FFFF","#808080","#66FF66","#FFFF00")
groups <- factor(pcaData$Gen)
levs <- levels(groups)
x=pcaData$PC1
y=pcaData$PC2
z=pcaData$PC3
groupCol <- c("#FF99FF","#33FFFF","#808080","#66FF66","#FFFF66")
# percentage variations for first 3 PCs have been taken from the Deseq data as 23, 16 and 15
percentVar <- c(23, 16, 15)
plot3d(pcaData[,1:3], col=colvec, size=2, type="s", lit=FALSE, xlim = c(-80,80), ylim=c(-50,50), zlim=c(-50,50), xlab=paste("PC1: 23% variance"), ylab=paste("PC2: 16% variance"), zlab=paste("PC3: 15% variance"))
#
for (i in 1:length(levs)) {
    group <- levs[i]
    selected <- groups == group
    xx <- x[selected]; yy <- y[selected]; zz <- z[selected]
    ellips <- ellipse3d(cov(cbind(xx,yy,zz)),
              centre=c(mean(xx), mean(yy), mean(zz)), level = 0.95)
    shade3d(ellips, col = groupCol[i], alpha = 0.1, lit = FALSE)}
# rsize manually to the rigth size and then create legen
legend3d("right", legend =legend, pch=16, col=legendCol, cex=2, inset=c(0.15))
dir.create("animation_merge")

for (i in 1:720) {
  view3d(userMatrix=rotationMatrix(2*pi * i/360, 0, 1, 0))
  rgl.snapshot(filename=paste("animation_merge/frame-",
                              sprintf("%03d", i), ".png", sep=""))}
# 
/SoftWare/Fiji.app/ImageJ-linux64
File > Import > Image Sequence
Image > type > 8-bit color
File >Save as > Animated gif


