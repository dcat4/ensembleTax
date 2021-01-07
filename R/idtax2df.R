#' Converts outputs of DECIPHER's idtaxa algorithm into a dataframe compatible
#' with the algorithms used in ensembleTax.

#' @details For consistency with DADA2's assignTaxonomy function, when used with
#' Silva, RDP, or GreenGenes it subsamples the ranks c("domain", "phylum",
#' "class", "order", "family", "genus"). Set db = NULL and supply ranks for
#' databases that aren't directly supported. The output taxonomy table is sorted
#' by the ASV-identifying data supplied in the rubric.
#'
#' CAUTION: the idtaxa algorithm does not return any ASV-identifying data in its
#' output "taxon" object. The elements of tt should thus be supplied in the same
#' order as the elements in rubric. This will typically be the case so long as
#' there is no tampering with the rubric or taxon object in between implementing
#' idtaxa and their use here.
#'
#' @author Dylan Catlett
#' @author Connie Liang
#'
#' @param tt The taxonomy table output by DECIPHER's idtaxa algorithm
#' @param db The database you ran idtaxa against. Either "pr2", "silva", "rdp",
#' or "gg" are supported.
#' @param ranks NULL, or a character vector of column names if db is set to NULL
#' @param boot The bootstrap threshold below which taxonomic assignments should
#' be set to NA. This can also be done with DECIPHER's idtaxa but is included
#' here for convenience.
#' @param rubric a DNAStringSet (see Biostrings package) with ASV sequences
#' named by your preferred ASV identifier. Both the ASV sequence and identifier
#' will be merged with the output dataframe. If NULL, ASV-identifying data is
#' not included in the output dataframe.
#' @param return.conf If TRUE, returns a list where the first element is your
#' formatted taxonomy table and the second element is a dataframe of bootstrap
#' confidence values. If FALSE, your formatted taxonomy table is returned as a
#' dataframe.
#'
#' @return a dataframe formatted for use with taxmapper and/or ensembleTax
#'
#' @seealso bayestax2df, ensembleTax, taxmapper
#'
#' @examples
#' data("idtax.pr2.sample")
#' data("rubric.sample")
#' head(idtax.pr2.sample)
#' head(rubric.sample)
#' df <- idtax2df(tt = idtax.pr2.sample, db = "pr2", ranks = NULL, boot = 0,
#' rubric = NULL, return.conf = FALSE)
#' head(df)
#' df <- idtax2df(tt = idtax.pr2.sample, db = "pr2", ranks = NULL, boot = 0,
#' rubric = rubric.sample, return.conf = FALSE)
#' head(df)
#' df <- idtax2df(tt = idtax.pr2.sample, db = "pr2", ranks = NULL, boot = 60,
#' rubric = rubric.sample, return.conf = FALSE)
#' head(df)
#' df <- idtax2df(tt = idtax.pr2.sample, db = "pr2", ranks = NULL, boot = 60,
#' rubric = rubric.sample, return.conf = TRUE)
#' head(df)
#'
#' @export
idtax2df <- function(tt, db = "pr2", ranks = NULL, boot = 0, rubric = NULL,
                     return.conf = FALSE){
  if (db == "pr2" || is.null(db)) {
    taxonomy<-c()
    conf <- c()
    notu <- length(tt)
    for(j in 1:notu){
      taxonomy<-append(taxonomy,tt[[j]]$taxon)
      conf <- append(conf,tt[[j]]$confidence)
    }
    yydf <- data.frame(matrix(unlist(taxonomy), nrow=notu, byrow=TRUE), stringsAsFactors = FALSE)
    yydf <- yydf[,-1]
    confdf <- data.frame(matrix(unlist(conf), nrow=notu, byrow=TRUE))
    confdf <- confdf[,-1]
    yydf[confdf < boot] <- NA

    if (db == "pr2") {
      ranks <- c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
      colnames(yydf) <- ranks
      colnames(confdf) <- ranks
    } else if (is.null(db)) {
      colnames(yydf) <- ranks
      colnames(confdf) <- ranks
    }

    if (!(is.null(rubric))) {
      rubdf <- data.frame(svN = names(rubric), ASV = as.character(rubric, use.names = FALSE), stringsAsFactors = FALSE)
      yydf <- cbind(rubdf, yydf)
      confdf <- cbind(rubdf, confdf)

      confdf <- ensembleTax::sort_my_taxtab(confdf, ranknames = ranks)
      yydf <- ensembleTax::sort_my_taxtab(yydf, ranknames = ranks)
    }

    rownames(yydf) <- NULL
    rownames(confdf) <- NULL

    # align formatting with taxmapper
    yydf <- base::apply(yydf, MARGIN = 2, FUN = as.character)
    confdf <- base::apply(confdf, MARGIN = 2, FUN = as.character)
    yydf <- base::as.data.frame(yydf, stringsAsFactors = FALSE)
    confdf <- base::as.data.frame(confdf, stringsAsFactors = FALSE)

    if (return.conf) {
      return(list(yydf, confdf))
    } else if (!return.conf) {
      return(yydf)
    }
  } else if (db == "silva" || db == "rdp" || db == "gg") {

    ranksub <- c("domain", "phylum", "class", "order", "family", "genus")

    # subset ranks from tax tab - adapted from dada2 tutorial
    tax <- t(sapply(tt, function(x) {
      m <- match(ranksub, x$rank)
      taxa <- x$taxon[m]
      taxa
    }))

    # subset ranks from bootstrap matrix - adapted from dada2 tutorial
    conf <- t(sapply(tt, function(x) {
      m <- match(ranksub, x$rank)
      conf <- x$confidence[m]
      conf
    }))
    # this loop fills in NAs in between names (introduced by rank subsampling):
    # does so by mirroring dada2's filler names in bayesian-silva taxonomy:
    for (i in 1:nrow(tax)) {
      nana <- which(is.na(tax[i,]))
      bubu <- which(!is.na(tax[i,]))
      if (any(nana > bubu)) {
        for (j in 1:length(nana)) {
          ah <- nana[j] - bubu
          if (any(ah < 0)) {
            # you need to replace it with unclassified_ closest upstream name:
            eh <- which(ah == min(ah[ah > 0])) # this is the index in bubu of the closest name...
            suff <- substr(ranksub[nana[j]], start = 1, stop = 2) # suffix based on rank you're filling in
            tax[i, nana[j]] <- paste0(tax[i, bubu[eh]], "_", suff)
            conf[i, nana[j]] <- conf[i, bubu[eh]] # propagate the confidence to these ranks too since you're not adding any real info...
          }
        }
      }
    }
    colnames(tax) <- ranksub
    colnames(conf) <- ranksub

    tax <- data.frame(tax, stringsAsFactors = FALSE)
    conf <- data.frame(conf, stringsAsFactors = FALSE)
    tax[conf < boot] <- NA # NA out below the boot threshold supplied

    # add rubric data
    if (!(is.null(rubric))) {
      rubdf <- data.frame(svN = names(rubric), ASV = as.character(rubric, use.names = FALSE), stringsAsFactors = FALSE)
      tax <- cbind(rubdf,tax)
      conf <- cbind(rubdf,conf)

      conf <- sort_my_taxtab(conf, ranknames = ranksub)
      tax <- sort_my_taxtab(tax, ranknames = ranksub)
    }

    rownames(tax) <- NULL
    rownames(conf) <- NULL

    # align formatting with taxmapper
    tax <- base::apply(tax, MARGIN = 2, FUN = as.character)
    conf <- base::apply(conf, MARGIN = 2, FUN = as.character)
    tax <- base::as.data.frame(tax, stringsAsFactors = FALSE)
    conf <- base::as.data.frame(conf, stringsAsFactors = FALSE)

    if (return.conf) {
      return(list(tax, conf))
    } else if (!return.conf) {
      return(tax)
    }
  }
}
