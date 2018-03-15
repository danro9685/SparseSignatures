library(data.table)
library(ggplot2)
library(gridExtra)

load("../data/mutation_categories.RData")

plotSignatures=function(beta, useColNames=T, firstBackground=T, xlabels=T){
  
  #Name the signatures
  if(firstBackground){
    rownames(beta) = c("Background", paste0("Signature ", 1:(nrow(beta)-1)))
  }else{
    rownames(beta) = paste0("Signature ", 1:nrow(beta))
  }
  
  #Name the categories if not given
  if(!useColNames){
    colnames(beta) = sort(mutation_categories$cat)
  }
  
  #Separate context and alteration
  x = as.data.table(melt(beta, varnames=c("signature", "cat")))
  x[, Context := paste0(substr(cat, 1,1), ".", substr(cat, 7,7))]
  x[, alt := paste0(substr(cat, 3,3), ">", substr(cat, 5,5))]
  
  #Plot
  glist = list()
  
  for(i in 1:nrow(beta)){
    
    plt = ggplot(x[signature == rownames(beta)[i]]) + 
      geom_bar(aes(x = Context, y = value, fill = alt), stat = "identity", position = "identity") + 
      facet_wrap(~alt, nrow=1, scales = "free_x") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      ggtitle(rownames(beta)[i]) + theme(legend.position="none")+ylab("Frequency of mutations")
    
    if(!xlabels){
      plt = plt+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    }
    
    glist[[i]] = plt
  } 
  
  grid.arrange(grobs = glist, ncol=ceiling(nrow(beta)/3)) 
  
}
