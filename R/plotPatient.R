library(data.table)
library(ggplot2)

#Plot patient counts
plotPatient = function(countMatrix, patientName, xlabels=T){
  
  #Separate contexts and alterations
  x = as.data.table(melt(countMatrix[patientName,]))
  x[, context := paste0(substr(variable, 1,1), ".", substr(variable, 7,7))]
  x[, alt := paste0(substr(variable, 3,3), ">", substr(variable, 5,5))]
  
  #Plot
  plt = ggplot(x) + 
    geom_bar(aes(x = context, y = value, fill = alt), stat = "identity", position = "identity") + 
    facet_wrap(~alt, nrow=1, scales = "free_x") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position="none")+ggtitle(patientName)
  
  if(!xlabels){
    plt = plt+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  return(plt)
  
}

