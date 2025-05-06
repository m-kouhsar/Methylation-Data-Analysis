
print("Loading libraries...")
suppressMessages(library(methylGSA))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library("reactome.db"))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19 ))
library("funr")
library(stringr)
source(paste0(str_remove(sys.script(),pattern="EWAS_methylGSA.Enrichment.R"),"methylGSA.methylRRA.R"))

args = commandArgs(T)

EWAS_Result_file <- args[1]
pval_col <- as.numeric(args[2])
sig.cpg_file <- args[3]
array.type <- args[4]
method <- args[5]
GS.type <- args[6]
minsize <- as.numeric(args[7])
maxsize <- as.numeric(args[8])
out_pref <- args[9]
cat('\n')
print("Arguments:")
print(paste("    EWAS file:",EWAS_Result_file))
print(paste("    Pvalue column:",pval_col))
print(paste("    Sig CpG file:",sig.cpg_file))
print(paste("    Array type:",array.type))
print(paste("    Enrichment method:",method))
print(paste("    GS type:",GS.type))
print(paste("    Min size:",minsize))
print(paste("    Max size:",maxsize))
print(paste("    Output prefix:",out_pref))
cat('\n')
print("Reading inputs...")
cpgs <- read.csv(file = EWAS_Result_file,header = T,stringsAsFactors = F,row.names = 1)
cpgs = cpgs[order(cpgs[,pval_col],decreasing=F),]
cpg.pval = cpgs[,pval_col]
names(cpg.pval) = rownames(cpgs)
sig.cpg= rownames(read.csv(file = sig.cpg_file,header = T,stringsAsFactors = F,row.names = 1))

print("Running enrichment analysis...")
results <- methylRRA(cpg.pval = cpg.pval,array.type = array.type,method = method,sig.cpg = sig.cpg,GS.type = GS.type,minsize = minsize,maxsize = maxsize)

print("Saving results...")

if(GS.type=="GO"){
  results.BP = results[results$ONTOLOGY=="BP",]
  results.MF = results[results$ONTOLOGY=="MF",]
  results.CC = results[results$ONTOLOGY=="CC",]
  
  write.csv(results.BP,file = paste0(out_pref,".methylGSA.GO.BP.min.",minsize,".max.",maxsize,".Enrichment.csv"))
  pdf(paste0(out_pref,".methylGSA.GO.BP.min",minsize,".max.",maxsize,".Enrichment.pdf"))
  print(barplot(results.BP, num = 20, colorby = "pvalue"))
  dev.off()
  write.csv(results.MF,file = paste0(out_pref,".methylGSA.GO.MF.min.",minsize,".max.",maxsize,".Enrichment.csv"))
  pdf(paste0(out_pref,".methylGSA.GO.MF.min.",minsize,".max.",maxsize,".Enrichment.pdf"))
  print(barplot(results.MF, num = 20, colorby = "pvalue"))
  dev.off()
  write.csv(results.CC,file = paste0(out_pref,".methylGSA.GO.CC.min.",minsize,".max.",maxsize,".Enrichment.csv"))
  pdf(paste0(out_pref,".methylGSA.GO.CC.min.",minsize,".max.",maxsize,".Enrichment.pdf"))
  print(barplot(results.CC, num = 20, colorby = "pvalue"))
  dev.off()
}else{
  write.csv(results,file = paste0(out_pref,".methylGSA.",GS.type,".min.",minsize,".max.",maxsize,".Enrichment.csv"))
  pdf(paste0(out_pref,".methylGSA.",GS.type,".min.",minsize,".max.",maxsize,".Enrichment.pdf"))
  print(barplot(results, num = 20, colorby = "pvalue"))
  dev.off()
}

