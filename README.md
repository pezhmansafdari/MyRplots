# MyRplots
A work flow of creating plots for genomic data

Similarity distance triangle heat map viruses. Betapartiti virus is the base line and distance matrix has been generated by MEGAN software. The data for creating this plot can be found in the data directory (Betapartitivirus3_0506.txt).
![Betapartiti](https://user-images.githubusercontent.com/5850834/110478112-aa232800-80ec-11eb-938b-9004dcee0879.jpg)

Phylogenetic tree of fungal OTUs. The OTU contigs have been assembled by pipits software. The multiple alignment has been done by MAFFT and the tree generaed by FastTree. 
![otu](https://user-images.githubusercontent.com/5850834/110479327-f15de880-80ed-11eb-9c3f-04ee99c90f19.jpg)

Animated 3D PCA plot with PCA data extracted from DeSeq2 after conducting DE Analysis. Due to the big size of the file, it is uploaded as a separate file(PCA.gif)

An ordination plot with vegan package and rda function. We are using raw counts mapped of RNAseq reads to the reference and the meta data file. Different combination of teh model like Phen or Gen * Phen can be modeled and plotted in a similar fashion. 
![RDA_genotype_1_2](https://user-images.githubusercontent.com/5850834/110600540-b6f75880-818c-11eb-8f07-77de21a2b685.jpg)
![RDA_genotype_1_3](https://user-images.githubusercontent.com/5850834/110600549-b959b280-818c-11eb-99a7-7cc994f9264c.jpg)


