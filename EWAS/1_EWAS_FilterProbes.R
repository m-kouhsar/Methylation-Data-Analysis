
message("Loading requiered libraries...")

suppressMessages(library(data.table))

args = commandArgs(T)

betas_file <- args[1]
out_prefix <- args[2]
ScriptDir <- args[3]
##########################################################################

message("Reading beta matrix...")
betas <- readRDS(betas_file)

message("Total number of probes:",dim(betas)[1],"\nTotal number of samples:",dim(betas)[2])

message("Reading Cross Hydridising Probes data..")
crosslist<-read.table(paste0(ScriptDir , "/References/CrossHydridisingProbes_McCartney.txt"), stringsAsFactors = FALSE)[,1]

message("Reading SNP Probes data..")
snpProbes<-read.table(paste0(ScriptDir , "/References/SNPProbes_McCartney.txt"), stringsAsFactors = FALSE, header = TRUE)
SNPlist<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),1]

message("Reading manifest file...")
manifest<-fread(paste0(ScriptDir , "/References/EPIC.V1.Manifest.tsv"), data.table = F)

message("**********************************************************************\n")

message("Step 1: Removing Cross Hydridising Probes ...")
index <- rownames(betas) %in% crosslist
betas.1 <-betas[!index,]
message("        ",sum(index)," probes are Cross Hydridising.")

message("Step 2: Removing SNP probes...")
index <- rownames(betas.1) %in% SNPlist
betas.1 <- betas.1[!index,]
message("        ",sum(index)," probes are SNPs.")

print("Step 3: Removing probes in the sex chromosomes...")
sexchrX<-manifest[which(manifest$CHR=="X"),]
sexchrY<- manifest[which(manifest$CHR== "Y"),]
index<- (rownames(betas.1) %in% sexchrX$IlmnID)|(rownames(betas.1) %in% sexchrY$IlmnID)
betas.1 <- betas.1[!index,]
message("        ",sum(index)," are in the sex chromosomes.")

print("Step 4: Removing probes with rs name...")
index <- substr(rownames(betas.1),1,2)=="rs"
betas.1 <- betas.1[!index,]
message("        ",sum(index)," probes start with 'rs'")

message("Final number of probes after filtering: ",nrow(betas.1))

print("Saving results...")

saveRDS(betas.1,file = paste0(out_prefix,'.betas.Filtered.rds'))

