# MyRplots
A work flow of creating plots for gene expression data
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
