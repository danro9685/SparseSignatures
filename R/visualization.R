# plot patients
"plotPatient" <- function(countMatrix, patientName, xlabels=TRUE, freq=FALSE) {
  
  # separate contexts and alterations
  x = as.data.table(suppressWarnings(melt(countMatrix[patientName,], variable.name = "cat", value.name = "N")))
  x[, Context := paste0(substr(cat, 1,1), ".", substr(cat, 7,7))]
  x[, alt := paste0(substr(cat, 3,3), ">", substr(cat, 5,5))]
  
  # convert to frequency if needed
  if(freq) {
    x[, freq:=N/x[, sum(N)]]
  }
  
  # make the plot
  if(freq) {
    plt = ggplot(x) + 
      geom_bar(aes(x = Context, y = freq, fill = alt), stat = "identity", position = "identity") + 
      facet_wrap(~alt, nrow=1, scales = "free_x") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(legend.position="none")+ggtitle(patientName)+ylab("Frequency of mutations")
  } else {
    plt = ggplot(x) + 
      geom_bar(aes(x = Context, y = N, fill = alt), stat = "identity", position = "identity") + 
      facet_wrap(~alt, nrow=1, scales = "free_x") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(legend.position="none")+ggtitle(patientName)+ylab("Number of mutations")
  }

  if(!xlabels) {
    plt = plt+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  return(plt)
  
}

# plot signatures
"plotSignatures" <- function(beta, useColNames=TRUE, firstBackground=TRUE, xlabels=TRUE) {
  
  # name the signatures
  if(firstBackground) {
    rownames(beta) = c("Background", paste0("Signature ", 1:(nrow(beta)-1)))
  }else {
    rownames(beta) = paste0("Signature ", 1:nrow(beta))
  }
  
  # name the categories if not given
  if(!useColNames) {
    colnames(beta) = sort(mutation_categories$cat)
  }
  
  # separate context and alteration
  x = as.data.table(melt(beta, varnames=c("signature", "cat")))
  x[, Context := paste0(substr(cat, 1,1), ".", substr(cat, 7,7))]
  x[, alt := paste0(substr(cat, 3,3), ">", substr(cat, 5,5))]
  
  # make the plot
  glist = list()
  
  for(i in 1:nrow(beta)) {
    
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
