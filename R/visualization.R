#' Plot trinucleotides counts for a set of given patients.
#'
#' @examples
#' data(patients)
#' patients.plot(trinucleotides_counts=patients,samples=c("PD10010a","PD10011a","PD10014a"))
#'
#' @title patients.plot
#' @param trinucleotides_counts trinucleotides counts matrix.
#' @param samples name of the samples. This should match a rownames in trinucleotides_counts.
#' @param freq boolean value; shall I display rates instead of counts?
#' @param xlabels boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export patients.plot
#' @import ggplot2
#' @import gridExtra
#' @importFrom reshape2 melt
#' @importFrom data.table as.data.table :=
#'
"patients.plot" <- function( trinucleotides_counts, samples = rownames(trinucleotides_counts), freq = FALSE, xlabels = FALSE ) {

    # make samples data
    trinucleotides_counts <- trinucleotides_counts[samples,,drop=FALSE]
    if(freq) {
        trinucleotides_counts <- trinucleotides_counts / rowSums(trinucleotides_counts)
    }

    # separate context and alteration
    x <- as.data.table(reshape2::melt(as.matrix(trinucleotides_counts),varnames=c("patient","cat")))
    x[,Context:=paste0(substr(cat,1,1),".",substr(cat,7,7))]
    x[,alt:=paste0(substr(cat,3,3),">",substr(cat,5,5))]

    # make the ggplot2 object
    glist <- list()
    for(i in 1:nrow(trinucleotides_counts)) {

        plt <- ggplot(x[patient==rownames(trinucleotides_counts)[i]]) + 
            geom_bar(aes(x=Context,y=value,fill=alt),stat="identity",position="identity") + 
            facet_wrap(~alt,nrow=1,scales="free_x") + 
            theme(axis.text.x=element_text(angle=90,hjust=1),panel.background=element_blank(),axis.line=element_line(colour="black")) + 
            ggtitle(rownames(trinucleotides_counts)[i]) + theme(legend.position="none") + ylab("Number of mutations")

        if(freq) {
            plt <- plt + ylab("Frequency of mutations")
        }

        if(!xlabels) {
            plt <- plt + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs=glist,ncol=ceiling(nrow(trinucleotides_counts)/3))

}


#' Plot the inferred mutational signatures.
#'
#' @examples
#' data(nmf_LassoK_example)
#' signatures.plot(beta=nmf_LassoK_example$beta)
#'
#' @title signatures.plot
#' @param beta matrix with the inferred mutational signatures.
#' @param useRowNames boolean value; shall I use the rownames from beta as names for the signatures?
#' @param xlabels boolean value; shall I display x labels?
#' @return A ggplot2 object.
#' @export signatures.plot
#' @import ggplot2
#' @import gridExtra
#' @importFrom data.table as.data.table :=
#'
"signatures.plot" <- function( beta, useRowNames = FALSE, xlabels = FALSE ) {

    # set names of the signatures
    if(!useRowNames) {
        rownames(beta) <- paste0("Signature ",1:nrow(beta))
    }

    # separate context and alteration
    x <- as.data.table(reshape2::melt(as.matrix(beta),varnames=c("signature","cat")))
    x[,Context:=paste0(substr(cat,1,1),".",substr(cat,7,7))]
    x[,alt:=paste0(substr(cat,3,3),">",substr(cat,5,5))]

    # make the ggplot2 object
    glist <- list()
    for(i in 1:nrow(beta)) {

        plt <- ggplot(x[signature==rownames(beta)[i]]) + 
            geom_bar(aes(x=Context,y=value,fill=alt),stat="identity",position="identity") + 
            facet_wrap(~alt,nrow=1,scales="free_x") + 
            theme(axis.text.x=element_text(angle=90,hjust=1),panel.background=element_blank(),axis.line=element_line(colour="black")) + 
            ggtitle(rownames(beta)[i]) + theme(legend.position="none") + ylab("Frequency of mutations")

        if(!xlabels) {
            plt <- plt + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
        }

        glist[[i]] <- plt

    }

    # make the final plot
    grid.arrange(grobs=glist,ncol=ceiling(nrow(beta)/3))

}
