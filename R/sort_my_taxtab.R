#' Sorts taxonomy table by ASV-identifying columns.
#'
#' @details A helper function for the ...2df family of pre-processing functions.
#' If multiple columns are available to sort, it uses the left-most column.
#'
#' @author Dylan Catlett
#'
#' @param tt A taxonomy table supplied as a dataframe of type char (no factors)
#' @param ranknames A character vector of the names of columns of tt that
#' contain taxonomic assignments. tt is sorted by columns not included in
#' ranknames.
#'
#' @seealso LCA2df, bayestax2df, idtax2df
#'
#' @examples
#' data("bayes.sample")
#' data("rubric.sample")
#' bayes.pretty <- bayestax2df(bayes.sample, rubric = rubric.sample)
#' sort_my_taxtab(bayes.pretty,
#' ranknames = c("kingdom", "supergroup", "division", "class", "order",
#' "family", "genus", "species"))
#'
#' @export
sort_my_taxtab <- function(tt, ranknames) {

  sort.by <- base::setdiff(colnames(tt), ranknames)
  if (length(sort.by)) {
    sort.by <- sort.by[1]
  }

  ii <- base::sort(tt[ , sort.by], index.return = TRUE)
  tt.out <- tt[ii$ix , ]
  return(tt.out)
}
