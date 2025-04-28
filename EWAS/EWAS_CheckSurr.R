args<-commandArgs(TRUE)

beta_file <- args[1] 
pheno_file <- args[2] 
trait <- args [3] 
covars_fact <- args[4] 
covars_num <- args[5] 
batchs <- args[6] 
model_combat <- args[7] 
model_sva <- args[8] 
max_sor_var <- args[9] 
model_lm <- args[10] 
out_pref <- args[11] 
calc_sv <- args[12]
sv_file <- args[13]

print("Arguments: ")
print(paste0("  beta value file= ",beta_file))
print(paste0("  Phenotype file= ",pheno_file))
print(paste0("  Trait var= ",trait))
print(paste0("  Factor covariates= ",covars_fact))
print(paste0("  Numeric covariates= ",covars_num))
print(paste0("  Batchs= ",batchs))
print(paste0("  Combat model= ",model_combat))
print(paste0("  SVA model= ",model_sva))
print(paste0("  Max number of sur var= ",max_sor_var))
print(paste0("  lm model= ",model_lm))
print(paste0("  Out files prefix= ",out_pref))
print(paste0("  Calculate Surrogate variables= ",if(as.numeric(calc_sv)==1) "Yes" else "No"))
if(as.numeric(calc_sv)==0){
  print(paste0("  Name of SV file to load= ",sv_file))
}

cat('\n')
print("Loading Libraries...")
suppressMessages(library("Gviz"))
suppressMessages(library(GenomicRanges))
suppressMessages(library(parallel))
suppressMessages(library(Biobase))
suppressMessages(library(sva))
suppressMessages(library(limma))
suppressMessages(library(qqman))
suppressMessages(library(stringr))
suppressMessages(library(WGCNA))
suppressMessages(library(ggplot2))

print("Reading Inputs...")
calc_sv <- as.numeric(calc_sv)
max_sor_var <- as.numeric(max_sor_var)
betas <- readRDS(beta_file)
pheno <- read.csv(file = pheno_file,row.names = 1,stringsAsFactors = F)

print("Is pheno and beta files match?")
identical(pheno$Basename,colnames(betas))

covars_fact <- str_split(covars_fact,pattern = ",",simplify = T)[1,]
covars_num <- str_split(covars_num,pattern = ",",simplify = T)[1,]
pheno.2 <- as.factor(pheno[,trait])

for (i in 1:length(covars_fact)) {
  pheno.2 <- cbind.data.frame(pheno.2,as.factor(pheno[,covars_fact[i]]))
}
for (i in 1:length(covars_num)) {
  pheno.2 <- cbind.data.frame(pheno.2,as.numeric(pheno[,covars_num[i]]))
}
names(pheno.2) <- c(trait,covars_fact,covars_num)

if(("sentrixid" %in% tolower(covars_fact))&("basename" %in% tolower(covars_fact))){
  i= match("sentrixid",tolower(names(pheno.2)))
  j= match("basename",tolower(names(pheno.2)))
  pheno.2[,i] = as.factor(str_split(pheno.2[,j],pattern = '_',simplify = T)[,1])
}

if(batchs!="0"){
  
  print("Removing batches...")
  batchs <- str_split(batchs,pattern = ",",simplify = T)[1,]
  if(model_combat!="0"){
    model_combat <- model.matrix(as.formula(model_combat),data = pheno.2)
  }else{
    model_combat <- NULL
  }
  combat1 <- betas
  for(i in 1:length(batchs)){
    combat1 <- ComBat(dat = combat1,batch = pheno.2[,batchs[i]],mod =model_combat )
  }
  betas <- combat1
  pheno.3 <- as.data.frame(lapply(pheno.2, as.numeric))
  print("PCA analysis after removing batches...")
  pca <- prcomp(betas,center = T,scale. = T,retx = F)
  pca1 <- as.data.frame(pca$rotation)
  pca1 <- cbind.data.frame(pca1,pheno.2)
  pca1$sdev <- pca$sdev
  pca1$eigs <- pca1$sdev^2
  pca1$Proportion = round(pca1$eigs/sum(pca1$eigs),digits = 4)*100
  pca1$Cumulative = cumsum(pca1$eigs)/sum(pca1$eigs)
  pca1$lab <- c(1:length(pca1$eigs))
  pca1 <- pca1[order(pca1$Proportion,decreasing=T),]
  max1 <- pca1$lab[1]
  max2 <- pca1$lab[2]
  
  pca2 <- as.data.frame(pca1[,c(1:10)])
  cor_ <- matrix(data = NA, ncol = dim(pheno.2)[2], nrow = ncol(pca2))
  colnames(cor_) <- names(pheno.2)
  rownames(cor_) <- paste0("PC",1:10)
  cor_pval <- cor_
  
  for(i in 1:ncol(pca2)){
    for (j in 1:ncol(pheno.3)) {
      res1<-cor.test(pca2[,i],pheno.3[,j], method="spearman",exact = FALSE)
      cor_[i,j]<-res1$estimate
      cor_pval[i,j]<-res1$p.value
    }
  }
  
  textMatrix = paste(signif(cor_, 2), "\n(",
                     signif(cor_pval, 1), ")", sep = "")
  dim(textMatrix) = dim(cor_)
  
  pdf(file = paste0(out_pref,"_Corrected_Batch_PCA.pdf"))
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
}

##SVA
if(calc_sv==1){
  print("Calculating surrogate variables...")
  mod1 = model.matrix(as.formula(model_sva), data=pheno.2)
  mod0 = cbind(mod1[,1])
  set.seed(1234)
  svs= sva::sva(betas,mod1,mod0)$sv
  svs<-data.frame(svs)
  names(svs) <- paste0("sv",1:dim(svs)[2])
  pheno <- pheno.2
  save(betas, pheno,svs, file=paste0(out_pref,"_sva.Rdata"))
}else{
  print("Loading surrogate variables...")
  load(str_trim(sv_file,side = "both"))
  pheno <- pheno.2
}
#######################################################################
#######################################################################

print("Running linear regression and saving results...")

EWAS <- function(x,model_lm,k,surro,pheno){
  x <- as.numeric(x)
  lm_vars <- all.vars(as.formula(model_lm))
  for(j in 1:length(lm_vars)){
    
    assign(lm_vars[j], pheno[,lm_vars[j]])
  }
  model_lm <- paste0("x",model_lm)
  if(k>0){
    model_lm <- paste0(model_lm,"+",paste(names(surro)[1:k],collapse = "+"))
    for(j in 1:k){
      
      assign(names(surro)[j], surro[,j])
    }
  }
  fit<-lm(formula = as.formula(model_lm),family = 'binomial')
  return(coef(summary(fit))[2,])
}

cl<- makeCluster(21)
inf_indx <- vector(mode = "numeric",length = max_sor_var+1)
best_sur <- 0
best_dist <- 100
for(i in 1:(max_sor_var+1)){
  
  EWAS_RESULTS<-t(parApply(cl, betas , 1, EWAS, model_lm,(i-1),svs,pheno))
  
  chisq <- qchisq(1-(EWAS_RESULTS[,4]),1)
  inf_indx[i] <- median(chisq)/qchisq(0.5,1)
  print(paste0("number of Sor var: ",(i-1)," Inflation index: ",inf_indx[i]))
  colnames(EWAS_RESULTS) <- c( "Estimate", "Std.Error"   , "t.value"  , "p.value")
  pdf(file = paste0(out_pref,'_qqplot_',round(median(chisq)/qchisq(0.5,1),digits = 3),'_SV',(i-1),'.pdf'))
  qq(EWAS_RESULTS[,4])
  dev.off()
  save(EWAS_RESULTS, file=paste0(out_pref,"_EWAS_RESULTS_SV",(i-1),".Rdata"))
  if(abs(inf_indx[i]-1)<best_dist){
    best_dist <- abs(inf_indx[i]-1)
    best_sur <- (i-1)
  }
}

stopCluster(cl) 

print(paste0("Best number of SVs is ",best_sur))

print("All done!")
#######################################################################
#######################################################################
