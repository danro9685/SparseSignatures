library(data.table)
library(deconstructSigs)
library(BSgenome.Hsapiens.1000genomes.hs37d5)

#Example data
ssm560 = fread("../data/ssm560.txt")
load("../data/patients.RData")

#List bases
bases = c("A", "C", "G", "T")
alts = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

#List categories
types = as.data.table(expand.grid(bases, alts, bases))
setnames(types, c("prevbase", "alt", "nextbase"))
types[, cat:=paste0(prevbase, "[", alt, "]", nextbase)]
types = types[order(cat)]

importData = function(input, sample = "sample", chr = "chrom", pos = "pos", ref = "ref", alt = "alt", bsg = NULL) 
{
  #Check that input is a data frame or data table
  if (class(input)=="data.frame") {
    inp = as.data.table(input)
    cat("found input data frame\n")
  }else if(class(input)[1]=="data.table"){
    inp = copy(inputname)
    cat("found input data table\n")
    }else if (file.exists(input)){
      inp = fread(input, header = TRUE, as.is = FALSE, check.names = FALSE)
    }else {
      stop("input is neither a file nor a loaded data frame")
    }
  
  #Filter input: recognized bases only
  inp = inp[, .(sample, chrom, pos, ref, alt)]
  inp = inp[ref %in% c("A", "T", "C", "G") & alt %in% c("A", "T", "C", "G")]
  
  #Check that bsg is a bsgenome object
  if (is.null(bsg)|class(bsg) != "BSgenome") {
      stop("The bsg parameter needs to be a BSgenome object.")
    }else{
      unknownChrom = inp[!chrom %in% GenomeInfoDb::seqnames(bsg)]
      if(nrow(unknownChrom)>0){
        warning(paste0("Check chromosome names -- not all match ", bsg, " object\n"))
        inp = inp[chrom %in% GenomeInfoDb::seqnames(bsg)]
      }
    }
  
  #Find context for each mutation
  inp$context = BSgenome::getSeq(bsg, chrom, pos-1, pos+1, as.character = T)
  
  #Check for mismatches with BSgenome context
  bad = inp[ref!=substr(context,2,2)]
  if(nrow(bad)>0){
    warning("Check ref bases -- not all match context:\n ")
  }
  
  #If ref is G or A, convert alteration
  inp[,alt := paste(ref, ">", alt, sep = "")]
  inp[,stdalt := alt]
  inp[ref %in% c("G", "A"), stdalt:=tolower(stdalt)]
  inp[ref %in% c("G", "A"), stdalt:=gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", stdalt))))]
  
  #If ref is G or A, take reverse complement of context
  
  inp[,stdcontext := context]
  mut$std.context[c(gind, tind)] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.context[c(gind, tind)]))))
  mut$std.context[c(gind, tind)] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.context[c(gind, tind)]))))
  mut$std.context[c(gind, tind)] <- sapply(strsplit(mut$std.context[c(gind, tind)], split = ""), function(str) {paste(rev(str), collapse = "")})
  mut$tricontext = paste(substr(mut$std.context, 1, 1), "[", mut$std.mutcat, "]", substr(mut$std.context, 3, 3), sep = "")
  final.matrix = matrix(0, ncol = 96, nrow = length(unique(mut[, sample.id])))
  colnames(final.matrix) = all.tri
  rownames(final.matrix) = unique(mut[, sample.id])
  for (i in unique(mut[, sample.id])) {
    tmp = subset(mut, mut[, sample.id] == i)
    beep = table(tmp$tricontext)
    for (l in 1:length(beep)) {
      trimer = names(beep[l])
      if (trimer %in% all.tri) {
        final.matrix[i, trimer] = beep[trimer]
      }
    }
  }
  final.df = data.frame(final.matrix, check.names = F)
  bad = names(which(rowSums(final.df) <= 50))
  if (length(bad) > 0) {
    bad = paste(bad, collapse = ", ")
    warning(paste("Some samples have fewer than 50 mutations:\n ", bad, sep = " "))
  }
  return(final.df)
}