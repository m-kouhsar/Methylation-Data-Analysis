args<-commandArgs(TRUE)

out_pref <- args[1]     # Output files prefix
run_lm <- args[2]       #1 means yes and 0 means No
model_lm <- args[3]     # Model for linear regression (If run_lm set to Yes)
num_sur <- args[4]      # Number of surrogate variables
sva_file <- args[5]     # Input files (an .Rdata file contains 3 dataframes: 
                                       # pheno, a dataframe contains clinical information
                                       # betas, a matrix or dataframe contains the beta values, the variables in lm model should exist in the column names of pheno
                                       # svs, a dataframe which its columns contain the surrugate variables)
ewas_file <- args[6]
EPIC_ref <- args[7]

print("Arguments: ")
print(paste0("  Output files prefix= ",out_pref))
print(paste0("  Number of Surrogate variable= ",num_sur))
print(paste0("  Run linear regression to generate EWAS results= ",if(as.numeric(run_lm)==1) "Yes" else "No"))
if(as.numeric(run_lm)==1){
  print(paste0("  Name of SV file to load= ",sva_file))
}
if(as.numeric(run_lm)==1){
  print(paste0("  lm model= ",model_lm))
}
if(as.numeric(run_lm)==0){
  print(paste0("  EWAS result file to load= ",ewas_file))
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
suppressMessages(library(data.table))

run_lm <- as.numeric(run_lm)

if(run_lm==1){
  print("Loading SVA file...")
  load(sva_file) 
  
  print("Running linear regression...")
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
  
  cl<- makeCluster(32)
  EWAS_RESULTS<-t(parApply(cl, betas , 1, EWAS, model_lm,num_sur,svs,pheno))
  stopCluster(cl) 
  colnames(EWAS_RESULTS) <- c( "Estimate", "Std.Error"   , "t.value"  , "p.value")
  table(EWAS_RESULTS[,4]<1e-4)
  
  chisq <- qchisq(1-(EWAS_RESULTS[,4]),1)
  print("Inflation Index:")
  median(chisq)/qchisq(0.5,1)
  
  print("Saving QQ plot...")
  pdf(file = paste0(out_pref,'_qqplot_',round(median(chisq)/qchisq(0.5,1),digits = 3),'_',num_sur,'.pdf'))
  qq(EWAS_RESULTS[,4])
  dev.off()
  
  print("Saving EWAS results...")
  save(EWAS_RESULTS, file=paste0(out_pref,"_EWAS_RESULTS_",num_sur,".Rdata"))
}else{
  print("Loading EWAS file...")
  load(ewas_file)
}
p.bonf<-p.adjust(EWAS_RESULTS[,4], "bonferroni")

EWAS_RESULTS<-cbind(EWAS_RESULTS,p.bonf)
EWAS_RESULTS<-data.frame(EWAS_RESULTS)
print("Summary of bonferroni pvalues:")
summary(EWAS_RESULTS$p.bonf)

print("Loading EPIC references...")
epicManifest<-fread(file=EPIC_ref, skip = 7,stringsAsFactors = F,fill = T)
ANNOT<-data.frame(epicManifest)
#head(ANNOT)
print("Adding annotations to EWAS results...")
x<-intersect(rownames(EWAS_RESULTS),ANNOT[,1])
ANNOT.sub<-ANNOT[which(ANNOT[,1] %in% x),]
EWAS_RESULTS<-EWAS_RESULTS[x,]

ANNOT.sub<-ANNOT.sub[order(ANNOT.sub[,1]),]
EWAS_RESULTS<-EWAS_RESULTS[order(rownames(EWAS_RESULTS)),]
#identical(rownames(EWAS_RESULTS), as.character(ANNOT.sub[,1]))
#[1] TRUE

############################################################
EWAS_RESULTS.ANNOT<-cbind(EWAS_RESULTS,ANNOT.sub)
#dim(EWAS_RESULTS.ANNOT)
#head(EWAS_RESULTS.ANNOT)
print("Saving annotated EWAS file...")
save(EWAS_RESULTS.ANNOT, file=paste0(out_pref,"_EWAS_RESULTS_Annot_",num_sur,".Rdata"))
#EWAS_RESULTS.ANNOT$p.bh <- p.adjust(EWAS_RESULTS.ANNOT[,4], "BH")

print("Number of DMPs with corrected pvalue < 0.05:")
length(which(EWAS_RESULTS.ANNOT$p.bonf<0.05))
#[1] 0
#bonf<-EWAS_RESULTS.ANNOT[which(EWAS_RESULTS.ANNOT$p.bonf<0.05),]
#table(as.character(bonf$UCSC_RefGene_Name))
#table(as.character(bonf$CHR))

############################################################
############################################################
print("Saving manhattan and QQ plots...")
manhattan <- data.frame("SNP"= as.character(EWAS_RESULTS.ANNOT$IlmnID),
                        "CHR"= 	as.numeric(as.character(EWAS_RESULTS.ANNOT$CHR)), 
                        "BP"= 	EWAS_RESULTS.ANNOT$MAPINFO,
                        "P"= 	EWAS_RESULTS.ANNOT[,4])
#manhattan <- manhattan[!is.na(manhattan$SNP)&!is.na(manhattan$CHR)&!is.na(manhattan$BP)&!is.na(manhattan$P),]
manhattan<-manhattan[order(manhattan$BP),]
manhattan<-manhattan[order(manhattan$CHR),]
#head(manhattan)

library(qqman)

tiff(filename = paste0(out_pref,"_EWAS.MANHATTAN_",num_sur,".tiff"), units="in", width=7, height=5, res=300)
manhattan(manhattan,annotatePval=3.162277e-06,annotateTop=F,suggestiveline = FALSE, genomewideline=(-log10(3.162277e-06)), cex.axis = 0.6 ,cex = 0.6, na.rm=na.omit,col = c("aquamarine4", "aquamarine1"))
dev.off()
tiff(filename = paste0(out_pref,"_EWAS.qqplot_",num_sur,".tiff"), units="in", width=3, height=3, res=300)
qq(manhattan$P, cex = 0.4, cex.axis = 0.6 , na.rm=na.omit)
dev.off()
############################################################
############################################################
print("Saving DMR file...")
EWAS_RESULTS.ANNOT$MAPINFO.1<-as.numeric(EWAS_RESULTS.ANNOT$MAPINFO)
EWAS_RESULTS.ANNOT$MAPINFO.1.1<-(EWAS_RESULTS.ANNOT$MAPINFO.1)+1

dmr <- data.frame("chrom" = 	paste0("chr", as.numeric(as.character(EWAS_RESULTS.ANNOT$CHR))), 
                  "start" = 	EWAS_RESULTS.ANNOT[rownames(EWAS_RESULTS.ANNOT), "MAPINFO.1"],
                  "end" 	= 	EWAS_RESULTS.ANNOT[rownames(EWAS_RESULTS.ANNOT), "MAPINFO.1.1"],
                  "pvalue" = 	EWAS_RESULTS.ANNOT[,4])

colnames(dmr) <- c("chrom", "start", "end", "pvalue")
dmr.10396<-dmr[order(dmr[,1],dmr[,2]),]

write.table(dmr.10396, file=paste0(out_pref,"_dmr_",num_sur,".txt"),sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)
#comb-p pipeline -c 4 --dist 500 --seed 1.0e-3 --anno hg19 -p  /mnt/data1/Ehsan/22q11DS/dmr.22q.cnv /mnt/data1/Ehsan/22q11DS/dmr.22q.del.txt
print("All done!")

