
print("Loading requiered libraries...")

suppressMessages(library(data.table))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))

args = commandArgs(T)

beta_file <- args[1]
pheno_file <- args[1]
variables_numeric <- args[3]
variables_factor <- args[4]
cross_hyd_ref <- args[2]
SNPPprob_ref <- args[3]
EPIC_manifest <- args[4]
out_prefix <- args[5]

##########################################################################

message("Reading beta values and phenotype data...")
betas <- readRDS(beta_file)
pheno <- read.csv(pheno_file , stringsAsFactors = F , row.names = 1)
if(!identical(rownames(pheno) , colnames(betas))){
  stop("Row names in Phenotype data are not matched with column names in beta matrix!")
}
message("Total number of probes:",dim(betas)[1],"\nTotal number of samples:",dim(betas)[2])

message("Reading CrossHydridising Probes data..")
cross<-read.table(cross_hyd_ref, stringsAsFactors = FALSE)
crosslist<-cross[,1]

message("Reading SNP Probes data..")
snpProbes<-read.table(SNPPprob_ref, stringsAsFactors = FALSE, header = TRUE)
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
SNPlist<-snpProbes[,1]

print("Reading EPIC manifest file...")
epicManifest<-fread(EPIC_ref, skip = 7,fill=TRUE)
ANNOT<-data.frame(epicManifest)

##########################################################

message("Step 1: Removing CrossHydridisingProbes ...")
index <- !(rownames(betas) %in% crosslist)
betas.1 <-betas[index,]
message("        ",sum(index),"probes removed")

message("Step 2: Removing probes with SNPs...")
index <- !(rownames(betas.1) %in% SNPlist)
betas.2 <- betas.1[ ,]
message("        ",sum(index),"probes removed")

print("Step 3: Removing probes in Chr X and Y...")
sexchrX<-ANNOT[which(ANNOT$CHR=="X"),]
sexchrY<- ANNOT[which(ANNOT$CHR== "Y"),]
index<- !((rownames(betas.2) %in% sexchrX$IlmnID)|(rownames(betas.3) %in% sexchrY$IlmnID))
betas.3 <- betas.2[index,]
message("        ",sum(index),"probes removed")

print("Step 4: Removing probes with rs name...")
snps<-substr(rownames(betas.3),1,2)
rsbn<-which(snps=="rs")
if(length(rsbn)!=0)
{
  betas.4<-betas.3[-rsbn,]
}else
{
  betas.4<-betas.3
}
message("        ",length(rsbn),"probes removed")

betas <- betas.4

message("Final number of probes in filtered beta matrix: ",nrow(betas))

print("Saving results...")

saveRDS(betas,file = paste0(out_file,'.prepared.rds'))

############################################################################################

variables_numeric <- str_split(variables_numeric,pattern = ',',simplify = T)[1,]
variables_factor <- str_split(variables_factor,pattern = ',',simplify = T)[1,]
pheno <- pheno[,c(variables_numeric,variables_factor)]

print("Calculating PCs...")
pca <- prcomp(t(betas),center = T,scale. = T)
pca1 <- pca$x[,1:10]
pca1 <- data.frame(pca1)

cor_ <- matrix(data = NA,nrow =10,ncol = length(c(variables_numeric,variables_factor)) )
colnames(cor_) <- c(variables_numeric,variables_factor)
rownames(cor_) <- colnames(pca1)[1:10]

cor_pval <- cor_
print("Calculatin correlations...")
for (i in 1:10) {
  for (j in 1:length(c(variables_numeric,variables_factor))) {
    var_ <- c(variables_numeric,variables_factor)[j]
    
    if(var_ %in% variables_factor){
      res1<-cor.test(as.numeric(pca1[,i]),as.numeric(as.factor(pheno[,var_])), method="spearman",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }else{
      res1<-cor.test(as.numeric(pca1[,i]),as.numeric(pheno[,var_]), method="pearson",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }
  }
}

textMatrix = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix) = dim(cor_)
print("Saving Plot...")
pdf(file = paste0(out_prefix,".PCA.pdf"))
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = rownames(cor_),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("PCA Analysis"))

dev.off()














