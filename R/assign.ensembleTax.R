#' Computes ensemble taxonomic assignments for each ASV in an amplicon data set
#'
#' @details
#' The algorithm takes as input a list of taxonomy tables (dataframes of type
#' character or list; no factors) and assumes rows correspond to ASVs/OTUs and
#' columns correspond to taxonomic assignments at ranks listed in descending
#' order in the input ranknames. All taxonomy tables should follow the same
#' taxonomic nomenclature (naming and ranking conventions), should include
#' ASV/OTU-identifying columns (e.g. ASV sequences or a column of asv numbers,
#' etc), and each row of each taxonomy table should represent the same ASV/OTU.
#' Use of the functions bayestax2df, idtax2df, and/or taxmapper will ensure your
#' taxonomy tables meet these requirements. Be advised that rownames of each
#' taxonomy table are set to NULL by assign.ensembleTax.
#'
#' Ensemble taxonomic assignments are computed by finding the highest-frequency
#' taxonomic assignment for each ASV across all input taxonomy tables. Several
#' parameters can be controlled by the user to weight the assignments of
#' specific taxonomy tables more highly than others (weights), to favor
#' assignments by a specific table in the event that multiple assignments are
#' found at the same (weighted) highest frequency (tiebreakz), to set a
#' (weighted) frequency threshold above which a taxonomic assignment must be
#' found to be assigned in the ensemble (assign.threshold), and finally to
#' ignore non-assignments signalled by NA in the frequency and assignment
#' computations (count.na).
#'
#' The output is a dataframe of ASVs and corresponding ensemble taxonomic
#' assignments.
#'
#' @author Dylan Catlett
#' @author Kevin Son
#'
#' @param x A list of dataframes of type character or list (no factors) that
#' contain an arbitrary number of meta-data columns (e.g. ASV sequences or
#' numbers), and other columns named according to ranknames that include
#' taxonomic assignments for each ASV in the data set
#' @param tablenames A character vector of the names of each taxonomy table
#' provided in x. Default is names(x)
#' @param ranknames The names of ranks (columns) of the taxonomy tables included
#' in x. These are used to track ASV-identifying data through the ensemble
#' calculations.
#' @param weights A numeric vector with length = length(x) that specifies
#' relative weights to the taxonomic assignments in the corresponding element of
#' x. Default is a vector with all elements =1 to specify equal weighting of
#' all taxonomy tables assignments. All values must be integers.
#' @param tiebreakz NULL is the default. Alternatively, a character vector
#' containing the tablenames in order of priority to be used as a tie-breaker
#' in the event that multiple taxonomic names are found at equal (weighted)
#' highest frequencies (above assign.threshold).
#' @param count.na TRUE or FALSE indicating whether you would like NA
#' assignments considered in the ensemble calculation. TRUE considers NA
#' assignments, FALSE does not consider NA assignments. assign.threshold is
#' implemented differently depending on whether this is TRUE or FALSE.
#' @param assign.threshold A number between 0 and 1 that indicates the
#' (weighted) proportion at which a particular taxonomic name must be assigned
#' in the input taxonomy tables in order to be assigned to the ensemble
#' taxonomic assignment. When count.na=FALSE, proportions are calculated only
#' relative to the number of tables with no NA assignments. When count.na=TRUE,
#' proportions are calculated relative to the sum of the weights argument.
#'
#' @return a dataframe containing ensemble taxonomic assignments
#'
#' @seealso idtax2df, bayestax2df, taxmapper
#'
#' @examples
#' fake1.pr2 <- data.frame(ASV = c("AAAA", "ATCG", "GCGC", "TATA", "TCGA"),
#'          kingdom = c("Eukaryota", "Eukaryota", "Eukaryota", "Eukaryota",
#'          "Eukaryota"),
#'          supergroup = c(NA, "Stramenopiles", "Rhizaria", "Stramenopiles",
#'          "Alveolata"),
#'          division = c(NA, "Ochrophyta", "Radiolaria", "Opalozoa",
#'          "Dinoflagellata"),
#'          class = c(NA, "Bacillariophyta", "Polycystinea", "MAST-12",
#'          "Syndiniales"),
#'          order = c(NA, "Bacillariophyta_X", "Collodaria", "MAST-12A", NA),
#'          family = c(NA, "Polar-centric-Mediophyceae", "Collophidiidae", NA,
#'          NA),
#'          genus = c(NA, NA, "Collophidium", NA, NA),
#'          species = as.character(c(NA, NA, NA, NA, NA)),
#'          stringsAsFactors = FALSE)
#' fake2.pr2 <- data.frame(ASV = c("AAAA", "ATCG", "GCGC", "TATA", "TCGA"),
#'          kingdom = c("Eukaryota", "Eukaryota", "Eukaryota", "Eukaryota",
#'          "Eukaryota"),
#'          supergroup = c(NA, "Stramenopiles", "Rhizaria", "Stramenopiles",
#'          "Alveolata"),
#'          division = c(NA, "Opalozoa", "Radiolaria", "Opalozoa",
#'          "Dinoflagellata"),
#'          class = c(NA, NA, "Polycystinea", NA, "Dinophycese"),
#'          order = c(NA, NA, "Collodaria", NA, NA),
#'          family = c(NA, NA, "Collophidiidae", NA, NA),
#'          genus = c(NA, NA, "Collophidium", NA, NA),
#'          species = as.character(c(NA, NA, NA, NA, NA)),
#'          stringsAsFactors = FALSE)
#' head(fake1.pr2)
#' head(fake2.pr2)
#' xx <- list(fake1.pr2, fake2.pr2)
#' names(xx) <- c("fake1", "fake2")
#' xx
#' eTax <- assign.ensembleTax(xx,
#'            tablenames = names(xx),
#'            ranknames = c("kingdom", "supergroup", "division","class","order",
#'            "family","genus","species"),
#'            tiebreakz = NULL,
#'            count.na=TRUE,
#'            assign.threshold = 0.5,
#'            weights=rep(1,length(xx)))
#' head(eTax)
#' eTax <- assign.ensembleTax(xx,
#'                     tablenames = names(xx),
#'                     ranknames = c("kingdom", "supergroup", "division",
#'                     "class","order","family","genus","species"),
#'                     tiebreakz = NULL,
#'                     count.na=FALSE,
#'                     assign.threshold = 0.5,
#'                     weights=c(2,1))
#' head(eTax)
#'
#' @export
assign.ensembleTax <- function(x, tablenames = names(x), ranknames = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                        weights=rep(1,length(x)), tiebreakz = NULL, count.na=TRUE, assign.threshold = 0) {
  names(x) <- tablenames

  # split meta-data from each table by ranknames:
  msplitz1 <- function(y) {
    eh <- colnames(y) %in% ranknames
    z <- y[ , eh, drop = FALSE]
    if (!is.null(row.names(z))) {
      row.names(z) <- NULL
    }
    z
  }
  msplitz2 <- function(y) {
    eh <- colnames(y) %in% ranknames
    zz <- y[ , !eh, drop = FALSE]
    if (!is.null(row.names(zz))) {
      row.names(zz) <- NULL
    }
    zz
  }
  metameta <- lapply(x, msplitz2)
  x <- lapply(x, msplitz1)

  # set up meta-data to merge back w/ ensemble later:
  metameta <- unique(metameta)
  if (length(metameta) > 1) {
    stop("The ASV identifiers supplied with each taxonomy table are not
         consistent. Please check that the ASV identifier columns of these
         taxonomy tables are identical and try again")
  }
  metameta <- metameta[[1]]

  ## set up
  n.rows <- nrow(x[[1]]) # get the number of ASV's aka number of rows
  n.cols <- ncol(x[[1]]) # get the number of columns to recreate ensemble taxonomy table
  n.dfs <- length(x) # getting the number of taxonomy tables inputted

  # set up empty consensus taxonomy table to add by each row
  consensus.tax <- base::data.frame(matrix(ncol=n.cols, nrow=0, dimnames=list(NULL, names(x[[1]]))))

  # determine if tiebreaks were specified
  if (is.null(tiebreakz)) {
    tiebreaker <- NA
  } else {
    # convert list of tiebreakers to a dataframe
    tiebreaker <- base::data.frame(table = tiebreakz, stringsAsFactors=FALSE)
    # set priority numbers
    tiebreaker$priority <- base::as.numeric(rownames(tiebreaker))
  }

  # iterate through each row and find the consensus of each row
  for (row in 1:n.rows) {
    # initialize the consensus row
    c.row <- base::data.frame(matrix(rep(NA, n.cols), ncol=n.cols, nrow = 1, dimnames=list(NULL, names(consensus.tax))))
    # create an alltax dataframe to track heirarchical assignments:
    alltax <- base::data.frame(matrix(NA, nrow = n.dfs, ncol = n.cols, dimnames=list(names(x), names(consensus.tax))), stringsAsFactors = FALSE)
    for (j in 1:n.dfs) {
      df <- x[[j]]
      alltax[j , ] <- df[row , ]
    }
    # create a collection of arrays to track the input tables used for assigning each rank
    # if a name that's about to be assigned does not come from a previous-used table, you need to break it:
    tmp.tblname <- tablenames
    for (col in ranknames) {
      # collects the taxs from each data frame
      taxs <- vector()
      # corresponding vector to know which df the tax came from
      df.idx <- vector()
      for (i in 1:length(tmp.tblname)) {
        df <- x[[tmp.tblname[i]]]
        # weights are represented by the amount of repeated taxs are
        taxs <- c(taxs, rep(df[row, col], weights[i]))
        df.idx <- c(df.idx, rep(tmp.tblname[i], weights[i]))
      }
      if (count.na) {
        freq.df <- base::as.data.frame(base::table(taxs, exclude=NULL), stringsAsFactors=FALSE)
      }
      else if (!count.na) {
        freq.df <- base::as.data.frame(base::table(taxs), stringsAsFactors=FALSE)
      }
      # if entires exist in the frequency table, determine the majority
      # else that means the NA's wasn't counted but will be set to NA as default
      if (nrow(freq.df) > 0) {
        # determine the proportion to compare to threshold by majority
        # if count.na you'll always consider all input tables for the proportion
        # if not you'll only consider those tables that were retained in the
        # frequency table
        if (count.na) {
          freq.df$prop <- freq.df$Freq / sum(weights)
        } else if (!count.na) {
          freq.df$prop <- freq.df$Freq / sum(freq.df$Freq)
        }
        # get the one with the greatest proportion aka the majority
        # see if the proportion of the majority is at least the threshold
        max.prop <- max(freq.df$prop)
        if (max.prop >= assign.threshold) {
          c.tax <- base::as.character(freq.df[which(freq.df$prop == max.prop), "taxs"])
          # if there are multiple taxs as the majority, we need to tie break it
          if (length(c.tax) > 1) {
            # see which data frame it is coming from
            # check if the data frame chosen is an option
            if (base::is.data.frame(tiebreaker)) {
              # create a data frame with taxs and tablename it came from
              pairs <- data.frame()
              for (i in 1:length(c.tax)) {
                # find the corresponding tablename of tax
                tax <- c.tax[i]
                idx <- taxs %in% tax
                tbl <- df.idx[idx]
                pairs <- base::rbind(pairs, cbind(tbl, rep(tax, times = length(tbl))))
              }
              # data frame of ties
              colnames(pairs) <- c("table","tax")
              # first assign exact matches
              exact <- base::merge(pairs, tiebreaker, by=c("table"), all.x=TRUE)

              # go back and consider taxnames only with NA
              alltt <- base::merge(exact, tiebreaker[ , c("table","priority")], by="table", all.x=TRUE)

              # resolve priorities
              alltt$priority <- dplyr::coalesce(alltt$priority.x, alltt$priority.y)
              alltt$priority.x <- NULL
              alltt$priority.y <- NULL

              # sort by priority to get tiebreaker at the top of the data frame
              sorted <- alltt[base::order(alltt$priority), ]
              # assign top row as consensus
              if (is.na(sorted[1, "priority"])) {
                c.row[, col] <- NA
              } else {
                c.row[, col] <- sorted[1, "tax"]
                bob <- base::data.frame(tblnam = tmp.tblname, tt = alltax[tmp.tblname , col], stringsAsFactors = FALSE)
                tina <- base::as.character(bob[which(bob$tt == c.row[, col]) , "tblnam"])
                tmp.tblname <- tina # this is now trimmed to only include tables who's names have been used at higher-order ranks
              }
            } else {
              # at this point just set it as NA
              c.row[, col] <- NA
            }
          } else {
            c.row[, col] <- c.tax
            bob <- base::data.frame(tblnam = tmp.tblname, tt = alltax[tmp.tblname , col])
            tina <- base::as.character(bob[which(bob$tt == c.row[, col]) , "tblnam"])
            tmp.tblname <- tina # this is now trimmed to only include tables who's names have been used at higher-order ranks
          }
        }
      }
      # if the assignment was NA, break (all downstream ranks are already NA)
      if (is.na(c.row[, col])) {
        break
      }
    }
    c.row <- base::data.frame(base::lapply(c.row, as.character), stringsAsFactors=FALSE)

    # after iterating through the columns add the consensus row to the data frame
    consensus.tax <- base::rbind(consensus.tax, c.row)
  }

  df <- base::cbind(metameta,consensus.tax)
  df <- base::data.frame(base::lapply(df, base::as.character), stringsAsFactors=FALSE)
  return(df)
}
