#' Maps an input taxonomy table onto a different taxonomic nomenclature.
#'
#' @author Dylan Catlett
#' @author Kevin Son
#'
#' @param tt The input taxonomy table you would like to map onto a new
#' taxonomic nomenclature. Should be a dataframe of type char or list (no
#' factors).
#' @param tt.ranks A character vector of the column names where taxonomic
#' names are found in tt. Supply them heirarchically (e.g. kingdom --> species)
#' @param tax2map2 The taxonomic nomenclature you would like to map onto. pr2
#' v4.12.0, Silva SSU v138 nr, GreenGenes v13.8 clustered at 97% similarity, and
#' the RDP train set 16 are included in the ensembleTax package. You can map to
#' these by specifying "pr2", "Silva", "gg", or "rdp". Otherwise should be a
#' dataframe of type character or list (no factors) with each column
#' corresponding to a taxonomic rank.
#' @param exceptions A character vector of taxonomic names at the basal/root
#' rank of tt that will be propagated onto the mapped taxonomy. ASVs assigned
#' to these names will retain these names at their basal/root rank in the mapped
#' taxonomy. All other ranks are assigned NA.
#' @param ignore.format If TRUE, the algorithm modifies taxonomic names in tt to
#' account for common variations in taxonomic name syntax and/or formatting
#' commonly encountered in reference databases (e.g. Pseudo-nitzschia will map
#' to Pseudonitzschia). If FALSE, formatting issues may preclude mapping of
#' synonymous taxonomic names (e.g. Pseudonitzschia will NOT map to
#' Pseudo-nitzschia). An exhaustive list of formatting details is included in
#' Details.
#' @param synonym.file If "default", taxmapper uses taxonomic synonyms included
#' with the ensembleTax package. If a custom taxonomic synonym file is
#' preferred, a string corresponding to the name of the csv file should be
#' supplied. Taxonomic synonyms are searched when exact name matches are not
#' found in tax2map2. ignore.format applies to synonyms if TRUE. Specify NULL if
#' you wish to forego synonym searches.
#' @param streamline If TRUE, only the mapped version of tt is returned as a
#' dataframe. If FALSE, a 3-element list is returned where element 1 is the
#' mapping key returned as a dataframe, element 2 is a character vector of all
#' names that could not be mapped (no exact matches found in tax2map2), and
#' element 3 is the mapped version of tt (a dataframe).
#' @param outfilez If NULL, mapping files are not saved to the current working
#' directory. Otherwise should be a 3-element character vector including, in
#' this order, the name of the file to store the taxonomic mapping key, the name
#' of the file to store the names that could not be mapped, and the name of the
#' file to store the ASVs supplied with tt with their mapped taxonomic
#' assignments. Each element of the vector should end in csv (only csv files
#' may be saved)
#'
#' @details Exceptions should be used when the user knows a particular taxonomic
#' group is not found in tax2map2. The user is responsible for supplying valid
#' taxonomic names as these must be found in tt and will be propagated as
#' given to all ASVs that are assigned this name in tt. This should only be
#' used for high-level taxonomic groups that are not found in a database (e.g.
#' for retaining Eukaryota when mapping onto a prokaryote-only taxonomic
#' nomenclature).
#'
#' When ignore.format = TRUE, names for which taxmapper cannot find exact
#' matches in tax2map2 are altered in case an exact match was not found due to
#' formatting issues. To do this taxmapper first checks for hyphens "-",
#' underscores "_", and single spaces " ". If these are found, variants of the
#' name with the hyphen/underscore/spaces replaced by each of the other two, as
#' well as all subnames spearated by these characters, and all subnames pasted
#' together with none of these special characters, are searched against tax2map2
#' for exact matches. It also creates all-lower and all-upper case versions of
#' these elements and again searches for exact name matches for these names.
#' To prevent matching of arbitrary names often used in reference databases like
#' "Clade X", after creating all of the above alternative names, those names
#' that begin with any variant of the words "clade" or "group" and those names
#' that are 2 characters or less are removed prior to re-searching tax2map2. All
#' alternative names created when ignore.format = TRUE are also searched for
#' synonyms in synonym.file. Be advised that setting ignore.format = TRUE does
#' not guarantee a more finely resolved mapped taxonomy table, and can actually
#' result in a less-resolved mapped taxonomy table in some circumstances.
#'
#' For high-throughput implementation of taxmapper, it's recommended to set
#' streamline = TRUE.
#'
#' @return If streamline = TRUE, a dataframe formatted for use with ensembleTax
#' that contains mapped taxonomic assignments for each ASV/OTU in the data set.
#'
#' If streamline = FALSE, a 3-element list where the first element is a
#' dataframe that contains all unique input taxonomic assignments and their
#' corresponding mapped outputs, the second element is a character vector that
#' contains all taxonomic names that could not be mapped, and the third element
#' contains mapped taxonomic assignments for each ASV in the data set.
#'
#' If is.null(outfilez) = FALSE, three csv files are saved in the current
#' working directory containing each of the three list elements above.
#'
#' @seealso idtax2df, bayestax2df, ensembleTax
#'
#' @examples
#' fake.silva <- data.frame(ASV = c("AAAA", "ATCG", "GCGC", "TATA", "TCGA"),
#' domain = c("Bacteria", "Eukaryota", "Eukaryota", "Eukaryota", "Eukaryota"),
#' phylum = c("Firmicutes", "Diatomea", "Retaria", "MAST-12", "Diatomea"),
#' class = c(NA, "Coscinodiscophytina_cl", "Polycystinea", "MAST-12A",
#' "Mediophyceae"),
#' order = c(NA, "Fragilariales", "Collodaria", NA, NA),
#' family = c(NA, "Fragilariales_fa", "Collodaria_fa", NA, NA),
#' genus = c(NA, "Podocystis", "Collophidium", NA, NA),
#' stringsAsFactors = FALSE)
#' head(fake.silva)
#' mapped.silva <- taxmapper(fake.silva,
#'                           tt.ranks = colnames(fake.silva)[2:ncol(fake.silva)],
#'                           tax2map2 = "pr2",
#'                           exceptions = c("Archaea", "Bacteria"),
#'                           ignore.format = FALSE,
#'                           synonym.file = "default",
#'                           streamline = TRUE,
#'                           outfilez = NULL)
#'
#' @export
taxmapper <- function(tt,
                      tt.ranks = colnames(tt),
                      tax2map2 = "pr2",
                      exceptions = c("Archaea", "Bacteria"),
                      ignore.format = FALSE,
                      synonym.file = "default",
                      streamline = TRUE,
                      outfilez = NULL) {

  if (length(tt.ranks) == ncol(tt)) {
    stop("You have not included any ASV-identifying data in your input
         taxonomy table. Please do this and try again.")
  }

  if (is.data.frame(tax2map2)){
    # do nothing
  } else if (tax2map2 == "pr2") {
    tax2map2 <- ensembleTax::pr2v4.12.0
  } else if (tax2map2 == "Silva") {
    tax2map2 <- ensembleTax::silva.nr.v138
  } else if (tax2map2 == "rdp") {
    tax2map2 <- ensembleTax::rdp_train_set_16
  } else if (tax2map2 == "gg") {
    tax2map2 <- ensembleTax::gg_13_8_train_set_97
  } else {
    stop("No valid tax2map2 object supplied.")
  }
  tax2map2.ranks <- colnames(tax2map2)

  # function to remove hyphens, underscores, upper case of name
  preprocessTax <- function(taxonomy) {
    alt.full <- c(stringr::str_replace_all(taxonomy, "-"," "),
                  stringr::str_replace_all(taxonomy, "-","_"),
                  stringr::str_replace_all(taxonomy, " ","-"),
                  stringr::str_replace_all(taxonomy, " ","_"),
                  stringr::str_replace_all(taxonomy, "_","-"),
                  stringr::str_replace_all(taxonomy, "_"," "))
    # split terms by hyphens
    no.hyphen <- base::strsplit(taxonomy, "-")
    # split terms by underscores
    no.underscore <- base::strsplit(taxonomy, "_")
    # split by space:
    no.spc <- base::strsplit(taxonomy, " ")
    # split terms by first instance of underscores and combine previous splits
    taxs <- c(taxonomy,
              alt.full,
              no.hyphen[[1]], no.underscore[[1]], no.spc[[1]],
              paste(no.hyphen[[1]], sep = '', collapse = ''),
              paste(no.underscore[[1]], sep = '', collapse = ''),
              paste(no.spc[[1]], sep = '', collapse = ''))
    # remove duplicates
    taxs <- base::unique(taxs)
    # convert all to lower and uppercase
    no.upper <- base::tolower(taxs)
    no.lower <- base::toupper(taxs)
    # create alternative suffixes for certain taxonomies
    final.taxs <- unique(c(taxs, no.upper, no.lower))
    # remove names that start with clade/group (and format variants), and/or
    # short names as these are likely to be numbers
    cl <- stringr::str_locate(final.taxs, "clade")
    cl2 <- stringr::str_locate(final.taxs, "Clade")
    cl3 <- stringr::str_locate(final.taxs, "CLADE")
    gr <- stringr::str_locate(final.taxs, "group")
    gr2 <- stringr::str_locate(final.taxs, "Group")
    gr3 <- stringr::str_locate(final.taxs, "GROUP")
    ll <- stringr::str_length(final.taxs)
    rm.rows <- c(which(cl[, "start"] == 1),
                 which(cl2[, "start"] == 1),
                 which(cl3[, "start"] == 1),
                 which(gr[, "start"] == 1),
                 which(gr2[, "start"] == 1),
                 which(gr3[, "start"] == 1),
                 which(ll <= 2))
    final.taxs <- final.taxs[-rm.rows]
    # this ensures that the original name is always the first one mapped:
    final.taxs <- c(taxonomy, final.taxs)
    return(final.taxs)
  }

  # function to search through the tax2map2 to find a match for the taxonomy name inputted
  findMapping <- function(taxonomy, tax2map2) {
    # iterate through the most specific ranking to the most generic ranking
    cols <- base::rev(names(tax2map2))
    for (i in 1:length(cols)) {
      matchings <- tax2map2[which(tax2map2[, cols[i]] == taxonomy), ] # find rows that match at that rank
      if (nrow(matchings) != 0) {
        # create respective row for the match
        # make everything downstream of the rank found to be NA's
        matched.row <- base::data.frame(matrix(rep(NA, length(cols)), ncol = length(cols), nrow = 1))
        colnames(matched.row) <- base::rev(cols)
        # grab only the first match found
        matched.row[1:(length(cols)-i+1)] <- matchings[1, ][1:(length(cols)-i+1)]
        return(matched.row)
      }
    }
    return(NA)
  }

  # function to search through the synonyms data frame to find synonyms for given taxonomy name
  getSynonyms <- function(taxonomy, syn.df) {
    if (is.null(syn.df)) {
      return(c(taxonomy))
    }
    found.rows <- syn.df[which(syn.df == taxonomy, arr.ind=TRUE)[,'row'],] # find rows for synonym
    if (length(found.rows) > 0) {
      # populate the taxonomy with its synonyms
      v <- as.character(as.matrix(found.rows))
      return (unique(c(taxonomy, v[!is.na(v)])))
    }
    else {
      # if no synonyms found, just return the taxonomy
      return(c(taxonomy))
    }
  }

  # rename tax2map2 columns for uniqueness
  colnames(tax2map2) <- base::paste("tax2map2", colnames(tax2map2), sep="_")

  # grab only the taxonomies part of the dataframes
  taxin.u <- base::unique(tt[, (names(tt) %in% tt.ranks)])
  tax2map2.u <- base::unique(tax2map2)

  # read in the synonyms file
  if (is.null(synonym.file)) {
    synonyms <- NULL
  } else if (synonym.file != "default") {
    synonyms <- utils::read.csv(synonym.file, stringsAsFactors = FALSE)
    synonyms <- synonyms[, colnames(synonyms)[startsWith(colnames(synonyms), "Name")]]
  } else if (synonym.file == "default") {
    synonyms <- ensembleTax::synonyms_20200816
    synonyms <- synonyms[, colnames(synonyms)[startsWith(colnames(synonyms), "Name")]]
  }

  taxin.cols <- base::rev(names(taxin.u))

  # keep track of the taxonomy names that are not mapped
  not.mapped <- vector()

  # finialized mapping table from taxin to tax2map2 with only the taxonomy names
  mapped <- base::data.frame(matrix(ncol=(ncol(taxin.u) + ncol(tax2map2.u)),nrow=0, dimnames=list(NULL, c(names(taxin.u), names(tax2map2.u)))))

  # iterate through each row and column of taxin data frame
  for (row in 1:nrow(taxin.u)) {
    # keep track of the most generic taxonomy name
    highest.tax <- taxin.u[row, taxin.cols[ncol(taxin.u)]]
    # see if it is in the exceptions to skip the row
    if (base::is.element(highest.tax, exceptions)) {
      # create a NA row assignment since part of exceptions
      null.row <- base::data.frame(matrix(rep(NA, ncol(tax2map2.u)), ncol = ncol(tax2map2.u), nrow = 1, dimnames=list(NULL, names(tax2map2.u))))
      null.row[1] <- highest.tax
      combined <- base::cbind(taxin.u[row, ], null.row)
      mapped <- base::rbind(mapped, combined)
    }
    else {
      for (col in 1:ncol(taxin.u)) {
        # keep track of the original taxonomy name
        orig.tax <- taxin.u[row, taxin.cols[col]]
        # process the name to get alternatives by igorning its format
        if (ignore.format) {
          pos.taxs <- preprocessTax(orig.tax)
          for(tax in pos.taxs) {
            pos.taxs <- c(pos.taxs, getSynonyms(tax, synonyms))
          }
          pos.taxs <- base::unique(c(orig.tax, getSynonyms(orig.tax, synonyms), pos.taxs))
        }
        else {
          pos.taxs <- c(orig.tax, getSynonyms(orig.tax, synonyms))
        }
        # flag to keep track of when the row is already matched
        matched <- FALSE
        # counter to keep track what column number we are on
        counter <- 1
        # iterate through all alternatives of the taxonomy name with original one first
        for (taxonomy in pos.taxs) {
          last <- FALSE
          if (counter == length(pos.taxs)) {
            last <- TRUE
          }
          if (!is.na(taxonomy)) {
            # find matching
            match <- findMapping(taxonomy, tax2map2.u)
            if (is.data.frame(match)) {
              combined <- cbind(taxin.u[row, ], match)
              mapped <- rbind(mapped, combined)
              matched <- TRUE
              break
            }
            else { # if no matching is found, add to not.mapped
              if (last) {
                not.mapped <- c(not.mapped, orig.tax)
              }
            }
          }
          counter <- counter + 1
        }
        # if a match is found, we can move onto the next row
        if (matched) {
          break
        }
      }
    }
  }

  # use the mapped table created to left join the original data frame with metadata
  asv.mapped <- base::merge(x=tt, y=mapped, by=colnames(taxin.u), all.x=TRUE)
  asv.mapped <- asv.mapped[ , !(colnames(asv.mapped) %in% colnames(taxin.u))]
  # remove the unique addition of tax2map2 column names
  colnames(asv.mapped) <- base::gsub("tax2map2_", "", colnames(asv.mapped))

  # filter out duplicates for not mapped taxonomy names
  not.mapped <- base::unique(not.mapped)

  zz <- base::apply(asv.mapped, MARGIN = 2, FUN = as.character)
  df <- base::as.data.frame(zz, stringsAsFactors = FALSE)
  asv.mapped <- df

  asv.mapped <- sort_my_taxtab(asv.mapped, ranknames = tax2map2.ranks)
  rownames(asv.mapped) <- NULL

  if (!(base::is.null(outfilez))) {
    utils::write.csv(mapped, outfilez[1], row.names=FALSE)
    not.mapped.df <- base::as.data.frame(not.mapped)
    utils::write.table(not.mapped.df, outfilez[2], row.names=FALSE, col.names=FALSE)
    utils::write.csv(asv.mapped, outfilez[3], row.names=FALSE)
  }

  if (streamline) {
    return(asv.mapped)
  } else {
    return(list(mapped, not.mapped, asv.mapped))
  }
}

