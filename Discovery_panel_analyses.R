##### 
##### 
##### 
# Analysis of oat discovery panel for doi:: 
# Genomic heritability and GWAS


#1. Load libraries
############
library(rrBLUP)
library(statgenGWAS)
#############

#2. Read in genotype and phenotype input files
#############
## phenotype data available https://doi.org/10.1093/genetics/iyaa043 in File S1
## genotype data available https://datacommons.cyverse.org/browse/iplant/home/shared/GoreLab/dataFromPubs/Brzozowski_OatMetabolome_2021 as "DiversityPanel_GBS.csv"
blups <- read.csv("./DiscoveryPanel_drBLUPs.csv")
snps <- read.csv("/Users/Lauren/Dropbox/Projects/Oat_postdoc/specialized_met_manuscript/DataFiles_forPublication/DiversityPanel_GBS.csv")
snps[1:10,1:10]
##############

#3. Calculate principal components
###################
gbs_pcs <- prcomp(scale(as.matrix(snps)), center = T, scale = T)
gbs_pcs <- gbs_pcs$x
################

#4. Calculate relationship matrix
################
panel_amat <- rrBLUP::A.mat(as.matrix(snps), min.MAF = 0.05)
###############

#5. Calculate metabolite heritability 
##############
met_h2 <- matrix(nrow = (ncol(blups)-1), ncol = 2); met_h2 <- as.data.frame(met_h2)
colnames(met_h2) <- c("met", "h2"); met_h2$met <- colnames(blups)[2:ncol(blups)]

for (i in 2:ncol(blups)) {
  ans_met <- rrBLUP::kin.blup(data=blups, geno="line", pheno=colnames(blups)[i-1], K=panel_amat)
  met_h2[i,2] <- ans_met$Vg / (ans_met$Vg + ans_met$Ve) }
###############


#6. Conduct GWAS
#################
#follow #https://cran.r-project.org/web/packages/statgenGWAS/vignettes/GWAS.html
## marker info available https://datacommons.cyverse.org/browse/iplant/home/shared/GoreLab/dataFromPubs/Brzozowski_OatMetabolome_2021 as "marker_genome_position.csv"

#6.1 marker map = chr pos are columns, makers are row names
marker_map_forSGGWAS  <- read.csv("./marker_genome_position.csv")
rownames(marker_map_forSGGWAS) <- marker_map_forSGGWAS $Marker
marker_map_forSGGWAS  <- marker_map_forSGGWAS [,c(4,3)]; colnames(marker_map_forSGGWAS) <- c("chr", "pos")

#6.2 marker matrix names of the markers in its column names and the genotypes in its row names. 
marker_matrix_forSGGWAS <- snps

#6.3 - pheno first column of all elements of pheno should be genotype and all the other columns should represent different traits. 
pheno_forSGGWAS  <- combined_drblups
colnames(pheno_forSGGWAS)[1] <- "genotype"
rownames(pheno_forSGGWAS) <- pheno_forSGGWAS$genotype

#6.4 kinship - from above
#panel_amat

#6.5 covar This data.frame has genotypes in its row names and the covariates in the column names
covar_forSGGWAS <- gbs_pcs[,c(1:5)]

#6.6 discovery panel GWAS
panel_gwas_res <- list()
for (i in 2:ncol(pheno_forSGGWAS)) {
  tempGData <- createGData(geno = marker_matrix_forSGGWAS, 
                           map = marker_map_forSGGWAS, 
                           kin = panel_amat, 
                           pheno = pheno_forSGGWAS[,c(1,i)], 
                           covar = covar_forSGGWAS)
  tempResults <- runSingleTraitGwas(gData = tempGData)
  panel_gwas_res[[i]] <- tempResults
  names(panel_gwas_res)[i] <- colnames(pheno_forSGGWAS)[i]
}
#############

