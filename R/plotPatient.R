library(data.table)
library(ggplot2)

#Plot patient counts
plotPatient = function(countMatrix, patientName, xlabels=T, freq=F){
  
  #Separate contexts and alterations
  x = as.data.table(suppressWarnings(melt(countMatrix[patientName,], variable.name = "cat", value.name = "N")))
  x[, Context := paste0(substr(cat, 1,1), ".", substr(cat, 7,7))]
  x[, alt := paste0(substr(cat, 3,3), ">", substr(cat, 5,5))]
  
  #Convert to frequency if needed
  if(freq){
    x[, freq:=N/x[, sum(N)]]
  }
  
  #Plot
  if(freq){
    plt = ggplot(x) + 
      geom_bar(aes(x = Context, y = freq, fill = alt), stat = "identity", position = "identity") + 
      facet_wrap(~alt, nrow=1, scales = "free_x") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(legend.position="none")+ggtitle(patientName)+ylab("Frequency of mutations")
  }else{
    plt = ggplot(x) + 
      geom_bar(aes(x = Context, y = N, fill = alt), stat = "identity", position = "identity") + 
      facet_wrap(~alt, nrow=1, scales = "free_x") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(legend.position="none")+ggtitle(patientName)+ylab("Number of mutations")
  }

  if(!xlabels){
    plt = plt+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  

  
  return(plt)
  
}

