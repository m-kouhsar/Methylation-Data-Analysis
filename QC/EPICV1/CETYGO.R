####################################################################
# Estimating cell types using CETYGO package
#
# Reference: https://github.com/ejh243/BrainFANS
#            https://github.com/ejh243/CETYGO
#
###################################################################

adultBrainCETYGO <- function(betas, cellType){
  
  predPropAll<-list() # store output in list
  CETYGOplots <- list() # store CETYGO score plots in list
  propPlots <- list() # store proportion plots in list
  counter = 1 # counter for storing list elements from for loop
  maxCETYGO <- 0 # vector to hold max CETYGO score from each model/ref
  minCETYGO <- 1 # vector to hold min CETYGO score from each model/ref
  
  # run CETYGO using both ANOVA and IDOL methods and all available cell ref panals 
  for(method in names(modelBrainCoef)){
    for(j in 1:length(modelBrainCoef[[method]])){
      if(!is.null(modelBrainCoef[[method]][[j]])){
        predPropAll[[method]][[j]]<-projectCellTypeWithError(betas, modelBrainCoef[[method]][[j]])
      }
    } 
  }
  
  # Find max and min CETYGO scores across all methods/ref panals to calibrate axis limits
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[i]]
    for(j in 1:length(type)){
      plotdf <- as.data.frame(type[[j]])
      
      if (ncol(plotdf) > 0) {
        
        if (max(plotdf$CETYGO) > maxCETYGO){
          maxCETYGO <- max(plotdf$CETYGO)}
        
        if (min(plotdf$CETYGO) < minCETYGO){
          minCETYGO <- min(plotdf$CETYGO)}
      }
    }
  }
  
  # create boxplots of CETYGO scores and proportions
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[i]]
    for(j in 1:length(type)){
      plotdf <- as.data.frame(type[[j]])
      
      if (ncol(plotdf) > 0) {
        
        # CETYGO score boxplot
        pCETYGO <- ggplot(plotdf, aes(factor(0), CETYGO))+
          geom_boxplot()+
          coord_cartesian(ylim=c(minCETYGO - 0.01, maxCETYGO + 0.01))+
          ggtitle(paste0(cellType, "_", names(predPropAll[i]), "-", j))+
          theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.title = element_text(size=6))+
          geom_hline(yintercept=0.07, linetype="dashed", color = "red")
        
        CETYGOplots[[counter]] <- pCETYGO
        
        # proportions boxplot
        plotdf <-as.data.frame(plotdf[,1:(ncol(plotdf)-2)])
        plotdf$Basename <- rownames(plotdf)
        plotdf <- reshape2::melt(plotdf, id.vars = "Basename")
        colnames(plotdf) <- c("Basename", "CellType", "Proportion")
        
        p <- ggplot(plotdf, aes(x=CellType, y = Proportion)) +
          geom_boxplot()+
          ggtitle(paste0(cellType, "_", names(predPropAll[i]), "-", j))+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        propPlots[[counter]] <- p
        
        counter = counter + 1
      }
    }
  }
  return(list(Plot1 = CETYGOplots , Plot2=propPlots , Data=predPropAll))
  
}


adultBloodCETYGO <- function(betas){
  
  # run CETYGO using blood ref panal
  rowIndex<-rownames(betas)[rownames(betas) %in% rownames(modelBloodCoef)]
  predProp<-as.data.frame(projectCellTypeWithError(betas, modelBloodCoef[rowIndex,]))


  # boxPlot CETYGO score
  pCETYGO <- ggplot(predProp, aes(factor(0), CETYGO))+
    geom_boxplot()+
  geom_hline(yintercept=0.07, linetype="dashed", color = "red")+
    xlab("")

  
  # box plot proportions
  plotdf <-as.data.frame(predProp[,1:(ncol(predProp)-2)])
  plotdf$Basename <- rownames(plotdf)
  plotdf <- reshape2::melt(plotdf, id.vars = "Basename")
  colnames(plotdf) <- c("Basename", "CellType", "Proportion")

  p <- ggplot(plotdf, aes(x=CellType, y=Proportion)) +
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  return(list(Plot1=pCETYGO , Plot2 = p , Data = predProp))

}
