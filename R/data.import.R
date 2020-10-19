#' Make trinucleotides counts matrix from input data for a given reference genome.
#'
#' @examples
#' \dontrun{
#' data(ssm560_reduced)
#' library("BSgenome.Hsapiens.1000genomes.hs37d5")
#' trinucleotides_counts = import.trinucleotides.counts(data=ssm560_reduced, 
#'      reference=BSgenome.Hsapiens.1000genomes.hs37d5)
#' }
#'
#' @title import.trinucleotides.counts
#' @param data a data.frame with variants having 6 columns: sample name, chromosome, start position, end position, ref, alt.
#' @param reference a BSgenome object with the reference genome to be used to retrieve flanking bases.
#' @return A matrix with trinucleotides counts per patient.
#' @export import.trinucleotides.counts
#' @importFrom data.table data.table dcast .N
#' @importFrom Biostrings DNAStringSet complement reverseComplement subseq
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom BSgenome getSeq
#'
"import.trinucleotides.counts" <- function( data, reference = NULL ) {

    # check that reference is a BSgenome object
    if(is.null(reference)|class(reference)!="BSgenome") {
        stop("The reference genome provided as input needs to be a BSgenome object.")
    }

    # preprocessing input data
    data <- as.data.frame(data)
    colnames(data) <- c("sample","chrom","start","end","ref","alt")

    # consider only single nucleotide variants involving (A,C,G,T) bases
    data <- data[which(data[,"start"]==data[,"end"]),,drop=FALSE]
    data <- data[which(as.matrix(data[,"ref"])%in%c("A","C","G","T")),,drop=FALSE]
    data <- data[which(as.matrix(data[,"alt"])%in%c("A","C","G","T")),,drop=FALSE]
    data <- data[,c("sample","chrom","start","ref","alt"),drop=FALSE]
    colnames(data) <- c("sample","chrom","pos","ref","alt")
    data <- unique(data)
    data <- data[order(data[,"sample"],data[,"chrom"],data[,"pos"]),,drop=FALSE]

    # convert data to GRanges
    data <- GRanges(data$chrom,IRanges(start=(data$pos-1),width=3),ref=DNAStringSet(data$ref),alt=DNAStringSet(data$alt),sample=data$sample)

    # check that all chromosomes match reference
    if(length(setdiff(seqnames(data),GenomeInfoDb::seqnames(reference)))>0) {
        warning("Check chromosome names, not all match reference genome.")
    }

    # find context for each mutation
    data$context <- getSeq(reference,data)

    # check for any mismatch with BSgenome context
    if(any(subseq(data$context,2,2)!=data$ref)) {
        warning("Check reference bases, not all match context.")
    }

    # get complements and reverse complements
    data$cref <- complement(data$ref)
    data$calt <- complement(data$alt)
    data$rccontext <- reverseComplement(data$context)

    # identify trinucleotides motif
    data$cat <- ifelse(data$ref%in%c("C","T"),paste0(subseq(data$context,1,1),"[",data$ref,">",data$alt,"]",subseq(data$context,3,3)), 
                       paste0(subseq(data$rccontext,1,1),"[",data$cref,">",data$calt,"]",subseq(data$rccontext,3,3)))

    # create 96 trinucleotides mutation categories
    categories_context <- NULL
    categories_alt <- rep(c(rep("C>A",4),rep("C>G",4),rep("C>T",4),rep("T>A",4),rep("T>C",4),rep("T>G",4)),4)
    categories_cat <- NULL
    cont <- 0
    for(i in c("A","C","G","T")) {
        for(j in 1:6) {
            for(k in c("A","C","G","T")) {
                cont <- cont + 1
                categories_context <- c(categories_context,paste0(k,":",i))
                categories_cat <- c(categories_cat,paste0(k,"[",categories_alt[cont],"]",i))
            }
        }
    }
    mutation_categories <- data.table::data.table(context=categories_context,alt=categories_alt,cat=categories_cat)
    
    # count number of mutations per sample for each category
    input1 <- data.table::data.table(mutation_categories[,"cat"])
    input2 <- data.table::data.table(sample=data$sample,cat=data$cat)
    input2 <- input2[,.N,by=.(sample,cat)]
    data <- merge(input1,input2,by="cat",all=TRUE)
    data <- data.table::dcast(data,sample~cat,value.var="N")
    data <- data[!is.na(sample),drop=FALSE]
    data[is.na(data)] <- 0

    # make trinucleotides counts matrix
    samples_names <- data$sample
    data <- as.matrix(data[,2:ncol(data),drop=FALSE])
    rownames(data) <- samples_names
    data <- data[sort(rownames(data)),,drop=FALSE]
    data <- data[,sort(colnames(data)),drop=FALSE]
    trinucleotides_counts <- array(0,c(nrow(data),96))
    rownames(trinucleotides_counts) <- rownames(data)
    colnames(trinucleotides_counts) <- sort(as.character(mutation_categories$cat))
    rows_contexts <- rownames(data)
    cols_contexts <- colnames(trinucleotides_counts)[which(colnames(trinucleotides_counts)%in%colnames(data))]
    trinucleotides_counts[rows_contexts,cols_contexts] <- data[rows_contexts,cols_contexts]

    # return trinucleotides counts matrix
    return(trinucleotides_counts)

}
