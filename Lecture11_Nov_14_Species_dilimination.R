###### Lecture 11. Species delimination using DAPC #######

### load packages for DAPC analysis
library(adegenet)
library(vcfR)

### load up the example data from the package 
data(dapcIllus)
class(dapcIllus)
names(dapcIllus)
x <- dapcIllus$a
x
grp <- find.clusters(x, max.n.clust=40) 
names(grp)
head(grp$grp, 10)
table(pop(x), grp$grp)
### visualization of the grouping of individuals 
table.value(table(pop(x), grp$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:6))
### applying DAPC function 
dapc1 <- dapc(x, grp$grp)
dapc1
scatter(dapc1)
#################################################
### now start with the real dataset
#### read the vcf data
vcfR <- read.vcfR("microsternarchus_unlinked.vcf")
popdata <- read.csv("microsternarchus.60missing.csv")

fpop <- popdata$INDV ### assign individual ID to each individual
fishdata<-vcfR2genind(vcfR,pop=fpop)
fishdata
set.seed(36)
fishgrp <-find.clusters(fishdata, max.n.clust=40,
                        n.start = 50,
                        n.iter = 10000)
table.value(table(pop(fishdata), fishgrp$grp), col.lab=paste("inf", 1:11), row.lab=paste("ori", 1:72))
fishdapc1 <- dapc(fishdata, fishgrp$grp)
scatter(fishdapc1) ### all data points are clustered together ###

fishdapc1$ind.coord

plot.df<-as.data.frame(fishdapc1$ind.coord)
ggplot(data=plot.df, aes(x=LD1, y=LD2))+
  geom_point(cex=1,alpha=0.2)+
  theme_classic()+
  labs(x="Linear discriminant 1",y="Linear discriminant 2")

##### let clean the vcf file a bit more, and put more prior grouping information using PCA ####

####################
summary(fishdata$tab)
fishgenetics<-fishdata$tab ### in this dataset, we have allele counts for each locus and each individual 

### we want to remove loci that has missing data (NA) for more than 7 individuals 
missingind<-c() ### count the total number of missing data (NAs) per loci 
for (j in 1:length(fishgenetics[1,])){
  missingind[j]<-length(fishgenetics[is.na(fishgenetics[,j]),j])  
}

fishgenetics2<-fishgenetics[,missingind<8] ### remove loci with missing data more than 7 individuals

### what if we want to remove individuals that has more than 20% of missing data across loci?

#### your code here 

### what if we want to remove loci that has allele frequency larger than 95% or lower than 5%? 

#### your code here 

#### impute the missing data using population average
fishdata_imputed <- fishgenetics2
for (i in 1:ncol(fishdata_imputed)) {
  fishdata_imputed[is.na(fishdata_imputed[, i]), i] <- mean(fishdata_imputed[, i], na.rm = TRUE)
}

pca_result <- prcomp(fishdata_imputed, scale. = TRUE, center = TRUE)

plot(pca_result$x[, 1], pca_result$x[, 2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCA of Allele Counts")

pca_result$x[1:5,1:5] ### take a look at the structure of PCA result, each column stands for each principle component 
summary(pca_result$x[,1])

### assign individuals with PC1<0, as group1, for other individuals, keep them individual ID, put this new grouping information in fpop

#### your code here

### using the new grouping information to run DAPC. 

#### your code here