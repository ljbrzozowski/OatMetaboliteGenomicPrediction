##### 
##### 
##### 
# Analysis of oat validation panel for doi:: 
# Multi-kernel genomic prediction


#1. Load libraries
############
library(rrBLUP)
library(BGLR)
#############

#2. Read in genotype and phenotype input files
#############
## phenotype data available {doi} as Supplemental Data 1 
## genotype data available https://datacommons.cyverse.org/browse/iplant/home/shared/GoreLab/dataFromPubs/Brzozowski_OatMetabolome_2021 as "ElitePanel_GBS.csv"

blups_MN <- read.csv("./DiscoveryPanel_drBLUPs_MN.csv")
blups_SD <- read.csv("./DiscoveryPanel_drBLUPs_SD.csv")
blups_WI <- read.csv("./DiscoveryPanel_drBLUPs_WI.csv")
snps <- read.csv("./ValidationPanel_GBS.csv")
##############


#3. Define k-fold cross validation sets
##############
cv_sets <- list()
for (i in 1:50) { cv_sets[[i]] <- sample(1:189, size=38, replace = F)} #189 individuals
###############

#4. Construct kernels
##############
computeGK <- function(X, bw){
  X <- scale(X, center = TRUE, scale = TRUE)
  D<- as.matrix(dist(X, method="euclidean"))^2
  D <- D/mean(D)
  GK <- exp(-bw * D)
  return(GK)
}

#whole genome, for GBLUP
GK_whole_genome <- computeGK(snps, bw=1)
ETA_gblup <- list(K1=list(K = GK_whole_genome, model="RKHS"))

#define kernels by SNP lists (not run)
GK_metabolite <- computeGK(snps[,which(colnames(snps) %in% kernel_snps)], bw=1) 
GK_not_metabolite <- computeGK(snps[,-c(which(colnames(snps) %in% kernel_snps))], bw=1) 
ETA_metabolite <- list(K1=list(K = GK_metabolite, model="RKHS"), K2=list(K = GK_not_metabolite, model="RKHS"))
###############


#5. Do k-fold cross validation
##############
#xval_rep = 50 # resampling runs
#niter = 1200 # iterations
#bi = 200 #burn in

for (j in 1:ncol(blups_MN)) {
  temp_resMN <- as.data.frame(matrix(nrow=xval_rep, ncol = 2)) ;colnames(temp_resMN) <- c("GBLUP", "Metabolite")
  temp_resSD <- as.data.frame(matrix(nrow=xval_rep, ncol = 2)); colnames(temp_resSD) <- c("GBLUP", "Metabolite")
  temp_resWI <- as.data.frame(matrix(nrow=xval_rep, ncol = 2)); colnames(temp_resWI) <- c("GBLUP", "Metabolite")
  tempFileName_MN <- paste0("./5foldCV_MN_", colnames(blups_MN)[j], ".csv")
  tempFileName_SD <- paste0("./5foldCV_SD_", colnames(blups_SD)[j], ".csv")
  tempFileName_WI <- paste0("./5foldCV_WI_", colnames(blups_WI)[j], ".csv")

  for (i in 1:xval_rep) {
    # MN
    y.na <- blups_MN; tst <- cv_sets[[i]]; y.na[tst,] <- NA
    gblub_cv<-BGLR(y=y.na[,j],ETA=ETA_gblup,nIter=niter, burnIn=bi,saveAt='temp'); temp_resMN[i,1] <- cor(blups_MN[tst,j], gblup_cv$yHat[tst], use="complete.obs")
    metabolite_cv<-BGLR(y=y.na[,j],ETA=ETA_metabolite,nIter=niter, burnIn=bi,saveAt='temp'); temp_resMN[i,2] <- cor(blups_MN[tst,j], metabolite_cv$yHat[tst], use="complete.obs")
    write.csv(temp_resMN,tempFileName_MN)
    
    # SD
    y.na <- blups_SD; tst <- cv_sets[[i]]; y.na[tst,] <- NA
    gblub_cv<-BGLR(y=y.na[,j],ETA=ETA_gblup,nIter=niter, burnIn=bi,saveAt='temp'); temp_resSD[i,1] <- cor(blups_SD[tst,j], gblup_cv$yHat[tst], use="complete.obs")
    metabolite_cv<-BGLR(y=y.na[,j],ETA=ETA_metabolite,nIter=niter, burnIn=bi,saveAt='temp'); temp_resSD[i,2] <- cor(blups_SD[tst,j], metabolite_cv$yHat[tst], use="complete.obs")
    write.csv(temp_resSD,tempFileName_SD)
    
    # WI
    y.na <- blups_WI; tst <- cv_sets[[i]]; y.na[tst,] <- NA
    gblub_cv<-BGLR(y=y.na[,j],ETA=ETA_gblup,nIter=niter, burnIn=bi,saveAt='temp'); temp_resWI[i,1] <- cor(blups_WI[tst,j], gblup_cv$yHat[tst], use="complete.obs")
    metabolite_cv<-BGLR(y=y.na[,j],ETA=ETA_metabolite,nIter=niter, burnIn=bi,saveAt='temp'); temp_resWI[i,2] <- cor(blups_WI[tst,j], metabolite_cv$yHat[tst], use="complete.obs")
    write.csv(temp_resWI,tempFileName_WI)
  }}
##################

#6. Estimate variance componenets
###########
#niter = 1200 # iterations
#bi = 200 #burn in
var_comps <- as.data.frame(matrix(nrow = ncol(blups_MN), ncol = 5)); colnames(var_comps) <- c("met", "var_MetKernel", "varRest", "varE", "DIC"); var_comps$met <- colnames(blups_MN)
MN_var_comps_list <- list(gblup=var_comps, metabolite=var_comps); SD_var_comps_list <- list(gblup=var_comps, metabolite=var_comps); WI_var_comps_list <- list(gblup=var_comps, metabolite=var_comps)


for (i in 1:ncol(blups_MN)) {

  fm_gblup=BGLR(y=blups_MN[,i],ETA=ETA_gblup,saveAt='temp_')
  MN_var_comps_list$gblup[i,3] <- fm_gblup$ETA$K1$varU
  MN_var_comps_list$gblup[i,4] <- fm_gblup$varE
  MN_var_comps_list$gblup[i,5] <- fm_gblup$fit$DIC
  
  fm_gblup=BGLR(y=blups_SD[,i],ETA=ETA_gblup,nIter=niter, burnIn=bi,saveAt='temp_')
  SD_var_comps_list$gblup[i,3] <- fm_gblup$ETA$K1$varU
  SD_var_comps_list$gblup[i,4] <- fm_gblup$varE
  SD_var_comps_list$gblup[i,5] <- fm_gblup$fit$DIC
  
  fm_gblup=BGLR(y=blups_WI[,i],ETA=ETA_gblup,nIter=niter, burnIn=bi,saveAt='temp_')
  WI_var_comps_list$gblup[i,3] <- fm_gblup$ETA$K1$varU
  WI_var_comps_list$gblup[i,4] <- fm_gblup$varE
  WI_var_comps_list$gblup[i,5] <- fm_gblup$fit$DIC
  
  fm=BGLR(y=blups_MN[,i],ETA=ETA_metabolite,nIter=niter, burnIn=bi,saveAt='temp_')
  MN_var_comps_list$metabolite[i,2] <- fm$ETA$K1$varU
  MN_var_comps_list$metabolite[i,3] <- fm$ETA$K2$varU
  MN_var_comps_list$metabolite[i,4] <- fm$varE
  MN_var_comps_list$metabolite[i,5] <- fm$fit$DIC
  
  fm=BGLR(y=blups_SD[,i],ETA=ETA_metabolite,nIter=niter, burnIn=bi,saveAt='temp_')
  SD_var_comps_list$metabolite[i,2] <- fm$ETA$K1$varU
  SD_var_comps_list$metabolite[i,3] <- fm$ETA$K2$varU
  SD_var_comps_list$metabolite[i,4] <- fm$varE
  SD_var_comps_list$metabolite[i,5] <- fm$fit$DIC
    
  fm=BGLR(y=blups_WI[,i],ETA=ETA_metabolite,nIter=niter, burnIn=bi,saveAt='temp_')
  WI_var_comps_list$metabolite[i,2] <- fm$ETA$K1$varU
  WI_var_comps_list$metabolite[i,3] <- fm$ETA$K2$varU
  WI_var_comps_list$metabolite[i,4] <- fm$varE
  WI_var_comps_list$metabolite[i,5] <- fm$fit$DIC
  
  saveRDS(MN_var_comps_list, "./MN_var_comps_list.RDS")
  saveRDS(SD_var_comps_list, "./SD_var_comps_list.RDS")
  saveRDS(WI_var_comps_list, "./WI_var_comps_list.RDS")
}
##########
