#' Converts the output of DADA2's assignTaxonomy, which implements a naive
#' bayesian classifier, into a dataframe compatible with the algorithms used in
#' ensembleTax
#'
#' @author Dylan Catlett
#' @author Connie Liang
#'
#' @details For consistency with dada2's assignTaxonomy function, when used with
#' Silva, RDP, or GreenGenes it subsamples the ranks c("domain", "phylum",
#' "class", "order", "family", "genus"). Set db = NULL and supply ranks for
#' databases that aren't directly supported. If a rubric is supplied with
#' ASV-identifying meta-data (this is highly recommended), the output taxonomy
#' table is sorted by the (first returned column of) ASV-identifying data.
#'
#' @param tt The taxonomy table output by DADA2's assignTaxonomy function.
#' @param db The database you ran assignTaxonomy against. Either "pr2", "silva",
#' "rdp", or "gg" are supported. You may set to NULL and include a character
#' vector of rank (column) names for other databases.
#' @param ranks NULL, or a character vector of column names if db is set to NULL
#' @param boot The bootstrap threshold below which taxonomic assignments should
#' be set to NA. This can also be done with DADA2's assignTaxonomy but is
#' included here for convenience.
#' @param rubric NULL, or a DNAStringSet (see Biostrings package) with ASV
#' sequences named by your preferred ASV identifier. Both the ASV sequence and
#' identifier will be merged with the output dataframe. If NULL, ASV-identifying
#' data are excluded in the output dataframe.
#' @param return.conf If TRUE, returns a list where the first element is your
#' formatted taxonomy table and the second element is a dataframe of bootstrap
#' confidence values. If FALSE, your formatted taxonomy table is returned as a
#' dataframe.
#'
#' @return a dataframe formatted for use with taxmapper and/or ensembleTax
#'
#' @seealso idtax2df, ensembleTax, taxmapper
#'
#' @examples
#' data("bayes.sample")
#' data("rubric.sample")
#' head(bayes.sample)
#' head(rubric.sample)
#' df <- bayestax2df(tt = bayes.sample, db = "pr2", boot = 0, rubric = NULL,
#' return.conf = FALSE)
#' head(df)
#' df <- bayestax2df(tt = bayes.sample, db = "pr2", boot = 0,
#' rubric = rubric.sample, return.conf = FALSE)
#' head(df)
#' df <- bayestax2df(tt = bayes.sample, db = "pr2", boot = 60,
#' rubric = rubric.sample, return.conf = FALSE)
#' head(df)
#' df <- bayestax2df(tt = bayes.sample, db = "pr2", boot = 60,
#' rubric = rubric.sample, return.conf = TRUE)
#' head(df)
#'
#' @export
bayestax2df <- function(tt, db = "pr2", ranks = NULL, boot = 0,
                        rubric = NULL, return.conf = FALSE){
  taxonomy <- tt[[1]]
  conf <- tt[[2]]
  taxdf <- base::data.frame(taxonomy, stringsAsFactors = FALSE)

  confdf <- base::data.frame(conf, stringsAsFactors = FALSE)
  taxdf[confdf < boot] <- NA

  if (db == "pr2") {
    ranks <- c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
    colnames(taxdf) <- ranks
    colnames(confdf) <- ranks
  } else if (db == "silva" || db == "rdp" || db == "gg") {
    ranks <- c("domain", "phylum", "class", "order", "family", "genus")
    colnames(taxdf) <- ranks
    colnames(confdf) <- ranks
  } else if (is.null(db)) {
    colnames(taxdf) <- ranks
    colnames(confdf) <- ranks
  }

  if (!(is.null(rubric))) {
    rubdf <- base::data.frame(svN = names(rubric), ASV = as.character(rubric, use.names = FALSE), stringsAsFactors = FALSE)
    taxdf$ASV <- rownames(taxdf)
    taxdf <- base::merge(taxdf, rubdf, by.x = "ASV", by.y = "ASV")
    confdf$ASV <- rownames(confdf)
    confdf <- base::merge(confdf, rubdf, by.x = "ASV", by.y = "ASV")

    colorder <- c(colnames(rubdf), ranks)
    taxdf <- taxdf[ , colorder]
    confdf <- confdf[ , colorder]

    taxdf <- ensembleTax::sort_my_taxtab(taxdf, ranknames = ranks)
    confdf <- ensembleTax::sort_my_taxtab(confdf, ranknames = ranks)
  }

  rownames(taxdf) <- NULL
  rownames(confdf) <- NULL

  # align formatting w taxmapper
  taxdf <- base::apply(taxdf, MARGIN = 2, FUN = as.character)
  confdf <- base::apply(confdf, MARGIN = 2, FUN = as.character)
  taxdf <- base::as.data.frame(taxdf, stringsAsFactors = FALSE)
  confdf <- base::as.data.frame(confdf, stringsAsFactors = FALSE)

  if (return.conf) {
    return(list(taxdf, confdf))
  } else if (!return.conf) {
    return(taxdf)
  }
}
