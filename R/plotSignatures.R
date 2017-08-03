require("data.table")
require("ggplot2")
require("gridExtra")

# function to plot extracted signatures
"plotSignatures" <- function( beta, patients_ids = colnames(beta), genomeFreq = TRUE ) {
  
  # set rownames and colnames
  if(genomeFreq) {
    rownames(beta) = c("Background",paste0("S",1:(nrow(beta)-1)))
  }
  else {
    rownames(beta) = paste0("S",1:nrow(beta))
  }
  colnames(beta) = patients_ids
  
  # make the plot
  x = as.data.table(melt(beta))
  x[,context := paste0(substr(Var2,1,1),".",substr(Var2,7,7))]
  x[,alt := paste0(substr(Var2,3,3),">",substr(Var2,5,5))]
  
  glist = list()
  
  for(i in 1:nrow(beta)) {
    sig = rownames(beta)[i]
    glist[[i]] = ggplot(x[Var1 == sig]) + 
      geom_bar(aes(x = context, y = value, fill = alt), stat = "identity", position = "identity") + 
      facet_grid(.~alt) + 
      theme(legend.position="none", 
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_rect(fill = "white", colour = NA), 
            panel.grid.major = element_line(colour = "grey80"), 
            panel.grid.minor = element_blank()) + 
      ggtitle(sig)
  } 
  grid.arrange(grobs=glist, ncol=ceiling(nrow(beta)/4))
}
