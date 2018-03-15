library(data.table)
library(Biostrings)
library(GenomicRanges)

load("../data/mutation_categories.RData")

#Example data
#ssm560 = fread("../data/ssm560.txt")
#load("../data/patients.RData")

#Input requirements
#input: a data.frame/data.table object or file with 5 columns: sample name, chromosome, position, ref, alt
#bsg: a BSgenome object for the reference genome. Chromosome nmaes must match the input table.

importData = function(input, bsg) 
{
  #Check that input is a data frame or data table
  if (!("data.frame" %in% class(input))){
    if (file.exists(input)){
      input = fread(input, header = TRUE, as.is = FALSE, check.names = FALSE)
    }else{
      stop("input is neither a file nor a loaded data frame")
    }
  }
  
  #Set column names
  colnames(input = c("sample", "chrom", "pos", "ref", "alt"))
  
  #Convert input to GRanges
  inp = GRanges(input$chrom, IRanges(start=input$pos-1, width=3), ref=DNAStringSet(input$ref), alt=DNAStringSet(input$alt), sample=input$sample)

  #Filter input: recognized bases only
  bad = which(!(inp$ref %in% c("A", "T", "C", "G") & inp$alt %in% c("A", "T", "C", "G")))
  if(length(bad)>0){
    warning(paste0("Removing ", length(bad), " entries containing nonstandard bases. \n"))
    inp = inp[inp$ref %in% c("A", "T", "C", "G") & inp$alt %in% c("A", "T", "C", "G")]
  }
  
  #Check that bsg is a bsgenome object
  if(is.null(bsg)|class(bsg) != "BSgenome") {
      stop("The bsg parameter needs to be a BSgenome object.")
    }
  
  #Check that all chromosomes match bsg
  if(length(setdiff(seqnames(inp), GenomeInfoDb::seqnames(bsg)))>0){
        warning(paste0("Check chromosome names -- not all match ", bsg, " object\n"))
        inp[seqnames(inp) %in% GenomeInfoDb::seqnames(bsg)]
    }

  #Find context for each mutation
  inp$context = getSeq(bsg, inp)

  #Check for mismatches with BSgenome context
  if(any(subseq(inp$context,2,2)!=inp$ref)){
    warning("Check ref bases -- not all match context\n ")
  }
  
  #Get complements and reverse complements
  inp$cref = complement(inp$ref)
  inp$calt = complement(inp$alt)
  inp$rccontext=reverseComplement(inp$context)
  
  #Identify motif
  inp$cat = ifelse(inp$ref %in% c("C", "T"), 
                     paste0(subseq(inp$context,1,1), "[", inp$ref, ">", inp$alt, "]", subseq(inp$context, 3, 3)),
                     paste0(subseq(inp$rccontext,1,1), "[", inp$cref, ">", inp$calt, "]", subseq(inp$rccontext, 3, 3)))
  
  #Count number of mutations per sample, category
  counts = merge(mutation_categories[, .(cat)], data.table(sample=inp$sample, cat=inp$cat)[, .N, by=.(sample, cat)], by="cat", all=T)
  counts = dcast(counts, sample~cat, value.var = "N")
  counts = counts[!is.na(sample)]
  counts[is.na(counts)]=0
  
  #Make count marix
  countMatrix = as.matrix(counts[, 2:ncol(counts)])
  rownames(countMatrix) = counts$sample
  
  #Order matrix columns alphabetically
  countMatrix = countMatrix[,order(colnames(countMatrix))]
  
  #Identify samples with <100 mutations
  if(!is.null(nrow(countMatrix))){
    bad = names(rowSums(countMatrix) <= 100)
  }else{
    bad = (sum(countMatrix) <= 100)
  }
  if (length(bad) > 0) {
    warning(paste("Some samples have fewer than 100 mutations:\n ", paste(bad, collapse = ", "), sep = " "))
  }
  
  #Return matrix
  return(countMatrix)
}