###=========================================================================================================================================###
# Modified gometh Function in missMethyl package to run enrichment analysis on a list of CpGs in KEGG, GO, Wikipathway and Reactome databases #
# By: Morteza P Kouhsar                                                                                                                       #
###=========================================================================================================================================###

suppressMessages(library(missMethyl))
suppressMessages(library(clusterProfiler))
suppressMessages(library(ReactomePA))

###############################################################################################################################################
new.getKEGG <- function () 
{
  GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = "hsa", 
                                           convert = TRUE)
  isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
  GeneID.PathID <- GeneID.PathID[!isna, ]
  ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2], sep = ".")
  d <- !duplicated(ID.ID)
  GeneID.PathID <- GeneID.PathID[d, ]
  PathID.PathName <- limma::getKEGGPathwayNames(species.KEGG = "hsa", 
                                                remove.qualifier = TRUE)
  GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by = "PathwayID")
  kegg <- tapply(GeneID.PathID$GeneID, GeneID.PathID$PathwayID, 
                 list)
  list(idList = kegg, idTable = PathID.PathName)
}

# environment(New.getKEGG) <- asNamespace('missMethyl')
# assignInNamespace(".getKEGG", New.getKEGG, ns = "missMethyl")
###############################################################################################################################################
new.gometh <- function (sig.cpg, all.cpg = NULL, collection = c("GO", "KEGG", "WP" , "RA"), 
          array.type = c("450K", "EPIC"), plot.bias = FALSE, prior.prob = TRUE, 
          anno = NULL, equiv.cpg = TRUE, fract.counts = TRUE, genomic.features = c("ALL", 
                                                                                   "TSS200", "TSS1500", "Body", "1stExon", "3'UTR", "5'UTR", 
                                                                                   "ExonBnd"), sig.genes = T,
          minSize = 10, maxSize = 500) #, Adj.method = "BH") 
{
  array.type <- match.arg(toupper(array.type), c("450K", "EPIC"))
  collection <- match.arg(toupper(collection), c("GO", "KEGG", "WP" , "RA"))
  genomic.features <- match.arg(genomic.features, c("ALL", 
                                                    "TSS200", "TSS1500", "Body", "1stExon", "3'UTR", "5'UTR", 
                                                    "ExonBnd"), several.ok = TRUE)
  if (length(genomic.features) > 1 & any(grepl("ALL", genomic.features))) {
    message("All input CpGs are used for testing.")
    genomic.features <- "ALL"
  }
  if (array.type == "450K" & any(grepl("ExonBnd", genomic.features))) {
    stop("'ExonBnd' is not an annotated feature on 450K arrays,\n\n           please remove it from your genomic.feature parameter\n\n           specification.")
  }
  if (collection == "GO") {
    go <- missMethyl:::.getGO()
    result <- missMethyl::gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = go$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(go$idTable, result, by.x = "GOID", by.y = "row.names")
    rownames(result) <- result$GOID
  }
  else if (collection == "KEGG") {
    kegg <- new.getKEGG()
    result <- missMethyl::gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = kegg$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(kegg$idTable, result, by.x = "PathwayID", 
                    by.y = "row.names")
    rownames(result) <- result$PathwayID
  } else if (collection == "WP"){
    
    wp.all <- clusterProfiler:::get_wp_data(organism = "Homo sapiens")
    
    wp = vector(mode = "list" , length = 2)
    names(wp) = c("idTable" , "idList")
    wp$idTable = data.frame(WPID = wp.all$wpid , Description = wp.all$name)
    wp$idTable = dplyr::distinct(wp$idTable)
    wp$idList = vector(mode = "list" , length = length(unique(wp.all$wpid)))
    names(wp$idList) = unique(wp.all$wpid)
    for (wpid in unique(wp.all$wpid)) {
      wp$idList[[wpid]] = wp.all$gene[wp.all$wpid == wpid]
    }
    rm(wp.all)
    
    result <- missMethyl::gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = wp$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(wp$idTable, result, by.x = "WPID", 
                    by.y = "row.names")
    rownames(result) <- result$WPID
  } else if (collection == "RA"){
    
    RA.all = ReactomePA:::get_Reactome_DATA(organism = "human")
    RA = vector(mode = "list" , length = 2)
    names(RA) = c("idTable" , "idList")
    RA$idTable = data.frame(ReactomeID = names(RA.all$PATHID2NAME) , Description = RA.all$PATHID2NAME)
    RA$idList = RA.all$PATHID2EXTID
    
    rm(RA.all)
    
    result <- missMethyl::gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = RA$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(RA$idTable, result, by.x = "ReactomeID", 
                    by.y = "row.names")
    rownames(result) <- result$ReactomeID
    
  }
  result = result[, -1]
  result <- result[result$N >= minSize & result$N <= maxSize, ]
  #result$FDR <- p.adjust(result$P.DE, method = Adj.method)
  ## Order by p-value
  result <- result[order(result$P.DE), ]
  result
}

################################################################
