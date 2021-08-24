
#Step 1 - Open new project with R


#Getting and setting working directory
getwd()
setwd("~/path/to/dir")

#Installing "packrat" for reproducablity
#packrat files need to be in current working directory
install.packages("packrat")
library("packrat")

#Initiating "packrat
packrat::init("~/path/to/dir")

#Loading Libraries
library(poppr)
library(BiocManager)
library(ape)
library(ggtree)
library(genepop)
library(factoextra)
library(LEA, lib.loc="~/packrat/lib/")
library(grid)
library(data.table)
library(cluster)
library(NbClust)
library(phytools)



#Refer to packrat.loc file
#Use to restore to older packages
packrat::restore()
#Use to Update packages and store version in "packrat"
packrat::snapshot()

###ANALYSIS###
##Bi-parental population##
mydat_table1 = read.table("~/Bi_parental_pop1.txt", header = TRUE, stringsAsFactors = FALSE)
mydat_table2 = read.table("~/Bi_parental_pop2.txt", header = TRUE, stringsAsFactors = FALSE)

#Creating Variables
accession_names1 = as.character(mydat_table1[,1])
accession_names2 = as.character(mydat_table2[,1])

ploidy1 = as.integer(mydat_table1[,3])
ploidy2 = as.integer(mydat_table2[,3])

allele_dat1 = mydat_table1[4:ncol(mydat_table1)]
allele_dat2 = mydat_table2[4:ncol(mydat_table2)]

mydat_genind1 = df2genind(allele_dat1, sep = ":", ind.names = accession_names1, ploidy = ploidy1)
mydat_genind2 = df2genind(allele_dat2, sep = ":", ind.names = accession_names2, ploidy = ploidy2)

recode_dat1 = recode_polyploids(mydat_genind1)
recode_dat2 = recode_polyploids(mydat_genind2)

ssr_replen1 = rep(0.0001, length(recode_dat1@all.names))
ssr_replen2 = rep(0.0001, length(recode_dat2@all.names))

#Creating dendrograms
mytree1 = bruvo.boot(recode_dat1, 
                     replen = ssr_replen1, 
                     tree = "njs",  
                     add = FALSE, 
                     loss = FALSE, 
                     sample = 2000, 
                     cutoff = 75, 
                     quiet = FALSE, 
                     showtree = FALSE)

mytree2 = bruvo.boot(recode_dat2, 
                     replen = ssr_replen1, 
                     tree = "njs",  
                     add = FALSE, 
                     loss = FALSE, 
                     sample = 2000, 
                     cutoff = 75, 
                     quiet = FALSE, 
                     showtree = FALSE)

#Plotting denrograms
den1<-ggtree(mytree1, size=.25,  branch.length = "none", layout="radial")+
  geom_tiplab2(size=3, hjust=-.05, fontface = "bold")+
  labs(title=("A"))+
  theme(plot.title = element_text(hjust = 0.05, vjust = -10, size = 15, face = "bold"))+
  xlim(NA, 10)

rmytree2<-phytools::reroot(mytree2, 22)


den2<-ggtree(rmytree2, size=.25,  branch.length = "none", layout="radial")+
  geom_tiplab2(size=3, hjust=-.05, fontface="bold")+
  labs(title=("B"))+
  theme(plot.title = element_text(hjust = 0.05, vjust = -10, size = 15, face = "bold"))+
  xlim(NA, 10)


biparental_dendros<-gridExtra::grid.arrange(den1, den2, top = textGrob("", vjust=.5, gp = gpar(fontsize = 15, fontface = 'bold')), layout_matrix = matrix(c(1,2), ncol=2, byrow = TRUE))

#writng file
ggsave(file="~/Figure_2.pdf", plot = biparental_dendros, device = "pdf")

##9-SSR fingerprintset in 629 samples##
#Importing data
Mydat = read.table("~/Raw_allele_data_grouped_pub.txt", header = TRUE, stringsAsFactors = FALSE)

#Creating variables 
#Selecting markers
Marker_names <- colnames(Mydat)
Marker_names<- Marker_names[4:length(Marker_names)]
#Selecting accessions
Accession_names <- Mydat[,1]
#Selecting ploidy
Ploidy <- as.vector(Mydat$Ploidy)
#Selecting allele data
Allele_dat <- Mydat[,4:ncol(Mydat)]
#Selecting populations
population <- Mydat[,2]
population = as.factor(population)
#Creating a genind object 
Mydat_genind = df2genind(Allele_dat, sep = ":", ind.names = Accession_names, ploidy = Ploidy, NA.char = "NA", pop = population)
#Recoding ploidy
Mydat_genind_recoded = recode_polyploids(Mydat_genind)
#Recoding SSR repeat length according to Metzger et al., 2016
ssr_replen = rep(0.0001, length(Mydat_genind@all.names))

#Creating dendrogram
Mytree = bruvo.boot(Mydat_genind_recoded, 
                     replen = ssr_replen, 
                     tree = "njs", 
                     add = F, 
                     loss = F, 
                     sample = 2000, 
                     cutoff = 75, 
                     quiet = FALSE, 
                     showtree = FALSE)

#Creating a distance matrix
Mydat_bruvo <- as.matrix(bruvo.dist(pop=Mydat_genind_recoded, 
                                    replen = ssr_replen, 
                                    add = F, 
                                    loss = F))

#Principal coodinate analysis (PCoA)
mds = cmdscale(as.matrix(Mydat_bruvo), eig=TRUE, x.ret=TRUE)

#Creating PCoA Variables
mds.values = as.vector(rownames(mds$points))
X=as.vector(mds$points[,1])
Y=as.vector(mds$points[,2])

#Creating PCoA data frame
mds.data = data.frame(Sample=mds.values, X=X, Y=Y)
rownames(mds.data) <- mds.data[,1]
mds.data[,1] <- NULL

#Geting variance explained
mds.var.per = round(mds$eig/sum(mds$eig)*100, 1)
scree_plot_dat = data.frame(PCoord = seq(1, length(mds.var.per)), per.var.explained = mds.var.per)

#Creating total within sum of squares plot 
fviz_nbclust(mds.data, FUNcluster = cluster::pam, method = c("wss"), k.max = 8, nboot = 100)

#Creating clustering based on the "majority rule"
res.nbclust <- NbClust(mds.data,
                       min.nc = 2, max.nc = 9, 
                       method = "centroid", index ="all")

#Creating clustering based of mediods
PCoA_Clusters = cluster::pam(mds.data,  k=3, diss = inherits(mds.data, "dist"))

#Creting cluster variable
my_clusters<-PCoA_Clusters$clustering

#Creating mediods clustering dataframe
myclusterdf<-as.data.frame(PCoA_Clusters$clustering, rownames=T)
names(myclusterdf)[1] <- "clusters"

setDT(myclusterdf, keep.rownames = TRUE)[]
names(myclusterdf)[1] <- "samples"

#Creating plotting variables
my_clusters<-as.factor(my_clusters)

pop_colors = c("grey20", "steelblue4", "red4")
xlab_name = paste("Principal Coord1 (", mds.var.per[1], "%)", sep = "")
ylab_name = paste("Principal Coord2 (", mds.var.per[2], "%)", sep = "")

#Creating PCoA and variance explained plots
p1<-ggplot(mds.data, aes(x=X ,y=Y))+
  theme_classic() +
  geom_point(aes(shape = population, color = my_clusters), size = 3) +
  scale_color_manual(values = pop_colors) +
  scale_shape_manual(values=c(20,3)) +
  labs(shape="Population", colour="K-mediods Cluster") +
  guides(shape = guide_legend(order = 1),  color = guide_legend(order = 2)) +
  xlab(xlab_name) +
  ylab(ylab_name) + 
  theme(axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14), 
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        legend.title = element_text(size = 14)) +
  scale_y_continuous(limits = c(-0.45, 0.50)) +
  scale_x_continuous(limits = c(-0.45, 0.8))

p2<-ggplot(head(scree_plot_dat,30), aes(x=PCoord, y=per.var.explained)) + 
  theme_classic() +
  geom_bar(stat = "identity",fill="grey") + 
  scale_y_continuous(limits = c(0, 50)) +
  geom_line(data=head(scree_plot_dat,30), aes(x=PCoord, y=per.var.explained)) +
  labs(x = "Principal Coordinates", y = "Variance Explained (%)") + 
  theme(axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14), 
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14))


g1<-gridExtra::grid.arrange(p1,p2, layout_matrix = matrix(c(1,2), ncol=2, byrow = TRUE))
gridExtra::grid.arrange(p1,p2, layout_matrix = matrix(c(1,2), ncol=2, byrow = TRUE))

ggsave(file="~/Figure_4.pdf", plot= g1, device = "pdf", dpi = 300)


#Craeting and manipulating full dendrogram
reroottree1 <- phytools::reroot(Mytree, 432)

t1<-ggtree(reroottree1 , branch.length='none', layout='radial', size = .05)

t2<-t1 %<+% myclusterdf +
  geom_tiplab(aes(color=as.factor(clusters), group=as.factor(clusters), fill = as.factor(clusters)),
              size=.6,  fontface="bold") + 
  theme()+
  labs(colour="K-means Clusters") +
  geom_treescale(fontsize =2, linesize = .25, offset = 1)


t3<-t2+ scale_color_manual(labels = c("1", "2", "3"), values = c("grey20", "steelblue4", "red4"))+
  theme(plot.margin=grid::unit(c(25,25,25,25), "mm"), legend.title = element_blank(), legend.position = "none")

t4<-flip(t3, 250, 669)

ggsave(file="~/S1_figure.pdf", plot= g1, device = "pdf", dpi = 300)

#Craeting and manipulating collapsed dendrogram
reroottree1 <- phytools::reroot(Mytree , 432)

ggtree(reroottree1, branch.length='none', layout='radial', size = .05)+geom_text(aes(label=node))



d1<-ggtree(reroottree1, branch.length='none', layout='radial', size = .04)%>% collapse(node = 614)+
  geom_point2(aes(subset=(node== 710)), shape=96, size=2, fill='steelblue4')+geom_treescale(fontsize =2, linesize = .25, offset = 1)

d2<-d1 %<+% myclusterdf +
  geom_tiplab(aes(color=as.factor(clusters), group=as.factor(clusters), fill = as.factor(clusters)),
              size=.6,  fontface="bold") + 
  theme()+
  labs(colour="K-means Clusters")

d3<-d2+ scale_color_manual(labels = c("1", "2", "3"),values = c("grey20", "steelblue4", "red4"))+
  theme(plot.margin=grid::unit(c(25,25,25,25), "mm"), legend.title = element_blank(), legend.position = "none")


ggsave(file="~/Figure_3.pdf", plot= t4, device = "pdf", dpi = 300)

#Creating sNMF plot 
Allele_dat = t(Mydat_genind_recoded$tab)

Allele_dat[is.na(Allele_dat)] = 9

write.table(Allele_dat, "/Users/mandiedriskill/Desktop/USDA/Hops/Final_scientific_documents/Data_availability/Hum_data.geno", row.names = FALSE, col.names = FALSE, sep = "")

obj.snmf <- snmf("/Users/mandiedriskill/Desktop/USDA/Hops/Final_scientific_documents/Data_availability/Hum_data.geno", K =1:10, project = "new", repetitions = 10, entropy = TRUE)


plot(obj.snmf, col = "black", pch = 19, cex = 1.2, font = 2, font.lab = 2, type = "b")

best = which.min(cross.entropy(obj.snmf, K = 3))
best

k_colors <- c("steelblue4", "red4", "grey20")

pdf("~/Figure_5.pdf", width = 11, height = 5) 
barchart(obj.snmf, K = 3, run = best, border = NA, space = 0, col = k_colors, xlab = "Individuals", ylab = "Ancestry Proportions", main = "Humulus Data (K = 3)", font.lab = 2, font = 2, las = 1) -> chart_order
dev.off() 

#Getting diversity measures (allelic and genotypic)
locus_table(x = Mydat_genind_recoded , index = "simpson", lev = "genotype")

locus_table(Mydat_genind_recoded)

poppr(Mydat_genind_recoded)

#Getting null allele frequencies
nulls(inputFile = "~/Raw_allele_data_grouped_pub.gen",
  outputFile = "",
  settingsFile = "",
  nullAlleleMethod = "",
  CIcoverage = 0.95,
  verbose = interactive()
)


##9-SSR and 25 KASP SNP comaparison in 190 samples##
#Importing 9-SSR and 25 KASP SNP data
mydatSSR = read.table("~/9SSR_190samples_withdiscrepancies.txt", header = TRUE, stringsAsFactors = FALSE)
mydatSNP = read.table("~/25_KASP_SNP.txt", header = TRUE, stringsAsFactors = FALSE)

#Selecting Marker names
marker_namesSSR <- colnames(mydatSSR)
marker_namesSNP <- colnames(mydatSNP)

marker_namesSSR<- marker_namesSSR[4:length(marker_namesSSR)]
marker_namesSNP<- marker_namesSNP[3:length(marker_namesSNP)]

#Selecting the accession names
accession_namesSSR <- mydatSSR[,1]
accession_namesSNP <- mydatSNP[,1]

#Selecting the ploidy values
ploidySSR <- as.vector(mydatSSR$Ploidy)
ploidySNP <- as.vector(mydatSNP$Ploidy)

#Selecting the SSR allele data
allele_datSSR <- mydatSSR[,4:ncol(mydatSSR)]
allele_datSNP <- mydatSNP[,3:ncol(mydatSNP)]

#populationSSR <- mydatSSR[,2]
#populationSNP <- mydatSNP[,2]

#Creating a genind object 
mydat_genindSSR = df2genind(allele_datSSR, sep = ":", ind.names = accession_namesSSR, ploidy = ploidySSR, NA.char = "NA")
mydat_genindSNP = df2genind(allele_datSNP, sep = ":", ind.names = accession_namesSNP, ploidy = ploidySNP, NA.char = "NA")

#Generating dendrogram with prevosti's genetic distance
#Opened case (#239) with "poppr" on the difference in graphs when updating to version >2.9.0.
#"poppr" requested to open case with "ape".
set.seed(9999)
mytreeSSR = aboot(mydat_genindSSR, 
                  tree = "njs", 
                  sample = 2000, 
                  cutoff = 75, 
                  distance = prevosti.dist, 
                  quiet = FALSE, 
                  showtree = FALSE)

set.seed(9999)
mytreeSNP = aboot(mydat_genindSNP, 
                  tree = "njs", 
                  sample = 2000, 
                  cutoff = 75, 
                  distance = prevosti.dist, 
                  quiet = FALSE, 
                  showtree = FALSE, 
                  missing = "ignore",
                  root = FALSE)

#Plotting dendogram graph
denSSR<-ggtree(mytreeSSR, size=.05)+
  geom_tiplab(size=.5)+
  ggplot2::xlim(NA,1)

denSSR

denSNP<-ggtree(mytreeSNP, size=.05)+
  geom_tiplab(size=.5)+
  ggplot2::xlim(NA,1)

denSNP


SSR_KASP_compare<-gridExtra::grid.arrange(denSSR, denSNP, top = textGrob("", vjust=.5, gp = gpar(fontsize = 15, fontface = 'bold')), layout_matrix = matrix(c(2,1), ncol=2, byrow = TRUE))

ggsave(file="/~/Fig1.eps", plot= SSR_KASP_compare, device = "pdf", dpi = 300)

#Importing 9-SSR and 25 KASP allele data - only unique genotypes
mydatSSR_grouped = read.table("~/9SSR_190samples_grouped.txt", header = TRUE, stringsAsFactors = FALSE)
mydatSNP_grouped = read.table("~/25_KASP_SNP_grouped.txt", header = TRUE, stringsAsFactors = FALSE)

#Selecting Marker names
marker_namesSSR_grouped <- colnames(mydatSSR_grouped)
marker_namesSNP_grouped <- colnames(mydatSNP_grouped)

marker_namesSSR_grouped<- marker_namesSSR_grouped[4:length(marker_namesSSR_grouped)]
marker_namesSNP_grouped<- marker_namesSNP_grouped[4:length(marker_namesSNP_grouped)]

#Selecting the accession names
accession_namesSSR_grouped <- mydatSSR_grouped[,1]
accession_namesSNP_grouped <- mydatSNP_grouped[,1]

#Selecting the ploidy values
ploidySSR_grouped <- as.vector(mydatSSR_grouped$Ploidy)
ploidySNP_grouped <- as.vector(mydatSNP_grouped$Ploidy)

#Selecting the SSR allele data
allele_datSSR_grouped <- mydatSSR_grouped[,4:ncol(mydatSSR_grouped)]
allele_datSNP_grouped <- mydatSNP_grouped[,4:ncol(mydatSNP_grouped)]

#selecting population data
populationSSR_grouped <- mydatSSR_grouped[,2]
populationSNP_grouped <- mydatSNP_grouped[,2]

#Creating a genind object 
mydat_genindSSR_grouped = df2genind(allele_datSSR_grouped, sep = ":", ind.names = accession_namesSSR_grouped, ploidy = ploidySSR_grouped, NA.char = "NA", pop = populationSSR_grouped)
mydat_genindSNP_grouped = df2genind(allele_datSNP_grouped, sep = ":", ind.names = accession_namesSNP_grouped, ploidy = ploidySNP_grouped, NA.char = "NA", pop = populationSNP_grouped)

mydat_genind_recodedSSR_grouped = recode_polyploids(mydat_genindSSR_grouped)

#Diversity measures in the 9-SSR and 25 KASP SNP samples
print("SSR")
locus_table(x = mydat_genind_recodedSSR_grouped , index = "simpson", lev = "genotype")
print("SNP")
locus_table(x = mydat_genindSNP_grouped , index = "simpson", lev = "genotype")

print("SSR")
locus_table(mydat_genind_recodedSSR_grouped)
print("SNP")
locus_table(mydat_genindSNP_grouped)


print("SSR")
poppr(mydat_genind_recodedSSR_grouped)
print("SNP")
poppr(mydat_genindSNP_grouped)

writeLines(capture.output(sessionInfo()), "~/session_info.txt")
 

