
print("Loading requiered libraries...")

suppressMessages(library(data.table))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))

args = commandArgs(T)

betas_file <- args[1]
manifest_file <- args[2]
out_prefix <- args[3]
ScriptDir <- args[4]
##########################################################################

message("Reading beta values and phenotype data...")
betas <- readRDS(betas_file)

message("Total number of probes:",dim(betas)[1],"\nTotal number of samples:",dim(betas)[2])

message("Reading CrossHydridising Probes data..")
crosslist<-read.table(paste0(ScriptDir , "/References/CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)[,1]

message("Reading SNP Probes data..")
snpProbes<-read.table(paste0(ScriptDir , "/References/SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
SNPlist<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),1]

message("Reading manifest file...")
manifest<-fread(EPIC_ref, skip = 7,fill=TRUE, data.table = F)

##########################################################

message("Step 1: Removing CrossHydridisingProbes ...")
index <- !(rownames(betas) %in% crosslist)
betas.1 <-betas[index,]
message("        ",sum(index),"/",nrow(betas)," probes removed")

message("Step 2: Removing probes with SNPs...")
index <- !(rownames(betas.1) %in% SNPlist)
betas.1 <- betas.1[ ,]
message("        ",sum(index),"/",nrow(betas)," probes removed")

print("Step 3: Removing probes in Chr X and Y...")
sexchrX<-manifest[which(manifest$CHR=="X"),]
sexchrY<- manifest[which(manifest$CHR== "Y"),]
index<- !((rownames(betas.2) %in% sexchrX$IlmnID)|(rownames(betas.3) %in% sexchrY$IlmnID))
betas.1 <- betas.1[index,]
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

# variables_numeric <- str_split(variables_numeric,pattern = ',',simplify = T)[1,]
# variables_factor <- str_split(variables_factor,pattern = ',',simplify = T)[1,]
# pheno <- pheno[,c(variables_numeric,variables_factor)]
# 
# print("Calculating PCs...")
# pca <- prcomp(t(betas),center = T,scale. = T)
# pca1 <- pca$x[,1:10]
# pca1 <- data.frame(pca1)
# 
# cor_ <- matrix(data = NA,nrow =10,ncol = length(c(variables_numeric,variables_factor)) )
# colnames(cor_) <- c(variables_numeric,variables_factor)
# rownames(cor_) <- colnames(pca1)[1:10]
# 
# cor_pval <- cor_
# print("Calculatin correlations...")
# for (i in 1:10) {
#   for (j in 1:length(c(variables_numeric,variables_factor))) {
#     var_ <- c(variables_numeric,variables_factor)[j]
#     
#     if(var_ %in% variables_factor){
#       res1<-cor.test(as.numeric(pca1[,i]),as.numeric(as.factor(pheno[,var_])), method="spearman",exact = FALSE)
#       cor_[i,var_]<-as.numeric(res1$estimate)
#       cor_pval[i,var_]<-as.numeric(res1$p.value)
#     }else{
#       res1<-cor.test(as.numeric(pca1[,i]),as.numeric(pheno[,var_]), method="pearson",exact = FALSE)
#       cor_[i,var_]<-as.numeric(res1$estimate)
#       cor_pval[i,var_]<-as.numeric(res1$p.value)
#     }
#   }
# }
# 
# textMatrix = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
# dim(textMatrix) = dim(cor_)
# print("Saving Plot...")
# pdf(file = paste0(out_prefix,".PCA.pdf"))
# par(mar = c(6, 8.5, 3, 3))
# labeledHeatmap(Matrix = cor_,
#                xLabels = colnames(cor_),
#                yLabels = rownames(cor_),
#                ySymbols = rownames(cor_),
#                colorLabels = FALSE,
#                colors = greenWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.5,
#                zlim = c(-1,1),
#                main = paste("PCA Analysis"))
# 
# dev.off()
# 

