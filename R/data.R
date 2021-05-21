#' All unique taxonomic assignments from the pr2 reference database v4.12.0
#'
#' @format ^ dataframe with 45352 rows and 8 columns
#' \describe{
#'   \item{kingdom}{kingdom assignments}
#'   \item{supergroup}{supergroup assignments}
#'   \item{division}{division assignments}
#'   \item{class}{class assignments}
#'   \item{order}{order assignments}
#'   \item{family}{family assignments}
#'   \item{genus}{genus assignments}
#'   \item{species}{species assignments}
#' }
"pr2v4.12.0"

#' All unique taxonomic assignments from the Silva SSU nr database v138
#'
#' @format ^ dataframe with 6011 rows and 6 columns
#' \describe{
#'   \item{domain}{domain assignments}
#'   \item{phylum}{phylum assignments}
#'   \item{class}{class assignments}
#'   \item{order}{order assignments}
#'   \item{family}{family assignments}
#'   \item{genus}{genus assignments}
#' }
"silva.nr.v138"

#' All unique taxonomic assignments from the RDP Train Set 16
#'
#' @format ^ dataframe with 2472 rows and 6 columns
#' \describe{
#'   \item{domain}{domain assignments}
#'   \item{phylum}{phylum assignments}
#'   \item{class}{class assignments}
#'   \item{order}{order assignments}
#'   \item{family}{family assignments}
#'   \item{genus}{genus assignments}
#' }
"rdp_train_set_16"

#' All unique taxonomic assignments from the GreenGenes v13.8 clusted at 97%
#'
#' @format ^ dataframe with 4163 rows and 7 columns
#' \describe{
#'   \item{domain}{domain assignments}
#'   \item{phylum}{phylum assignments}
#'   \item{class}{class assignments}
#'   \item{order}{order assignments}
#'   \item{family}{family assignments}
#'   \item{genus}{genus assignments}
#'   \item{species}{genus assignments}
#' }
"gg_13_8_train_set_97"

#' Taxonomic synonyms searched by the taxmapper algorithm
#'
#' @format ^ dataframe with 174 rows and 11 columns
#' \describe{
#'   \item{Name_1}{first synonym}
#'   \item{Name_2}{second synonym}
#'   \item{Name_3}{third synonym}
#'   \item{Name_4}{fourth synonym}
#'   \item{Name_5}{fifth synonym}
#'   \item{Name_6}{sixth synonym}
#'   \item{Name_7}{seventh synonym}
#'   \item{References}{Reference for some synonyms}
#'   \item{Notes.References}{Notes from references}
#'   \item{X}{Additional references for some synonyms}
#'   \item{X.1}{Additional references for some synonyms}
#' }
"synonyms_v2"

#' Example rubric with ASV-identifying data
#'
#' @format ^ DNAStringSet with 5 elements
#' \describe{
#' \item{sv1}{sample ASV 1}
#' \item{sv2}{sample ASV 2}
#' \item{sv3}{sample ASV 3}
#' \item{sv4}{sample ASV 4}
#' \item{sv5}{sample ASV 5}
#' }
"rubric.sample"

#' Example output of dada2 assignTaxonomy function
#'
#' @format ^ list with 2 elements
#' \describe{
#' \item{tax}{taxonomic assignments}
#' \item{boot}{bootstrap confidence estimates}
#' }
"bayes.sample"

#' Example output of DECIPHER idtaxa function with pr2 taxonomy
#'
#' @format ^ list with 5 elements
#' \describe{
#' \item{1}{tax data for ASV 1}
#' \item{2}{tax data for ASV 2}
#' \item{3}{tax data for ASV 3}
#' \item{4}{tax data for ASV 4}
#' \item{5}{tax data for ASV 5}
#' }
"idtax.pr2.sample"

#' Example output of DECIPHER idtaxa function with silva taxonomy
#'
#' @format ^ list with 5 elements
#' \describe{
#' \item{1}{tax data for ASV 1}
#' \item{2}{tax data for ASV 2}
#' \item{3}{tax data for ASV 3}
#' \item{4}{tax data for ASV 4}
#' \item{5}{tax data for ASV 5}
#' }
"idtax.silva.sample"

