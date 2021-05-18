## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  # dada2 assignTaxonomy:
#  taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#  # this is optional and was not used in the example data below:
#  taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz")
#  
#  # DECIPHER idtaxa:
#  library(DECIPHER); packageVersion("DECIPHER")
#  dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
#  load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
#  ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors

## ---- eval=FALSE--------------------------------------------------------------
#  rubric <- DNAStringSet(getSequences(seqtab.nochim))
#  # this creates names (sv1, sv2, ..., svX) for each ASV
#  snam <- vector(mode = "character", length = length(rubric))
#  for (i in 1:length(rubric)) {
#      snam[i] <- paste0("sv", as.character(i))
#  }
#  names(rubric) <- snam

## -----------------------------------------------------------------------------
library("ensembleTax")
library("Biostrings")

data("idtax.pr2.sample")
data("idtax.silva.sample")
data("bayes.sample")
data("rubric.sample")

## ---- eval = FALSE------------------------------------------------------------
#  head(idtax.pr2.sample)
#  head(idtax.silva.sample)
#  head(bayes.sample)
#  head(rubric.sample)

## -----------------------------------------------------------------------------
idtax.pr2.pretty <- idtax2df(idtax.pr2.sample, 
                             db = "pr2", 
                             ranks = NULL,
                             boot = 50,
                             rubric = rubric.sample,
                             return.conf = FALSE)
idtax.silva.pretty <- idtax2df(idtax.silva.sample, 
                             db = "silva", 
                             ranks = NULL,
                             boot = 50,
                             rubric = rubric.sample,
                             return.conf = FALSE)
bayes.pr2.pretty <- bayestax2df(bayes.sample, 
                             db = "pr2", 
                             ranks = NULL,
                             boot = 50,
                             rubric = rubric.sample,
                             return.conf = FALSE)

## -----------------------------------------------------------------------------
# remove ASV columns from tax tables for easier viewing:
idtax.pr2.pretty <- idtax.pr2.pretty[ , -which(names(idtax.pr2.pretty) %in% c("ASV"))]
idtax.silva.pretty <- idtax.silva.pretty[ , -which(names(idtax.silva.pretty) %in% c("ASV"))]
bayes.pr2.pretty <- bayes.pr2.pretty[ , -which(names(bayes.pr2.pretty) %in% c("ASV"))]

head(idtax.pr2.pretty)
head(idtax.silva.pretty)
head(bayes.pr2.pretty)

## -----------------------------------------------------------------------------

idtax.silva.mapped2pr2 <- taxmapper(idtax.silva.pretty,
                      tt.ranks = colnames(idtax.silva.pretty)[2:ncol(idtax.silva.pretty)],
                      tax2map2 = "pr2",
                      exceptions = c("Archaea", "Bacteria"),
                      ignore.format = TRUE,
                      synonym.file = "default",
                      streamline = TRUE,
                      outfilez = NULL)
head(idtax.silva.mapped2pr2)


## -----------------------------------------------------------------------------

xx <- list(idtax.pr2.pretty, idtax.silva.mapped2pr2, bayes.pr2.pretty)
names(xx) <- c("idtax-pr2", "idtax-silva", "bayes-pr2")
eTax1 <- assign.ensembleTax(xx, 
                     tablenames = names(xx), 
                     ranknames = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                     tiebreakz = NULL, 
                     count.na=TRUE, 
                     assign.threshold = 0, 
                     weights=rep(1,length(xx)))
# show the 3 individuals again for easy viewing:
lapply(xx, FUN = head)
head(eTax1)


## ---- eval=FALSE--------------------------------------------------------------
#  install.packages(c("ggplot2", "reshape2"))

## -----------------------------------------------------------------------------
library("ggplot2")
library("reshape2")

nasum <- function(taxdf){
    notuz <- nrow(taxdf)
    x <- is.na(taxdf)
    ii <- colSums(x) / notuz
    return(ii)
}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
# create a list with individuals and ensemble(s):
xx <- list(idtax.pr2.pretty, idtax.silva.mapped2pr2, bayes.pr2.pretty, eTax1)
names(xx) <- c("idtax-pr2", "idtax-silva", "bayes-pr2", "ensemble")
# remove meta-data columns (NOTE the hard-coded column removal!!!)
xx <- lapply(xx, function(x) x[, -c(1)])
tina <- lapply(xx, nasum)
yaboi <- matrix(unlist(tina), nrow=length(xx), byrow=TRUE)
rownames(yaboi) <- names(xx)
colnames(yaboi) <- colnames(xx[[1]])
yaboi <- as.data.frame(t(yaboi))
yaboi$rankz <- rownames(yaboi)
yaboi <- melt(yaboi, id.vars = "rankz")
plt <- ggplot(yaboi, aes(fill = variable, x = rankz, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Taxonomic Rank", y = "Proportion of ASVs Unassigned") + 
    scale_x_discrete(limits = colnames(xx[[1]])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy table",
                      breaks = names(xx),
                      values = cbPalette)
print(plt)

## ---- warning=FALSE, message = FALSE------------------------------------------
library("dplyr")
# fcn that compares two taxonomy tables:
tblcomper <- function(y,x) {
  perf <- dplyr::intersect(x, y) # perfectly-matching rows (assignments)

  tmp.x <- dplyr::setdiff(x, perf)
  tmp.y <- dplyr::setdiff(y, perf)
  if (!identical(tmp.x[, 1], tmp.y[, 1])) {
    tmp.x <- ensembleTax::sort_my_taxtab(tmp.x, ranknames = colnames(tmp.x)[3:ncol(tmp.x)])
    tmp.y <- ensembleTax::sort_my_taxtab(tmp.y, ranknames = colnames(tmp.y)[3:ncol(tmp.x)])
  }
  xna <- is.na(tmp.x)
  yna <- is.na(tmp.y)
  # This warning happens 1 million times when there's no NA's. It doesn't matter b/c inf is always bigger so suppressing.
  ## Warning in min(which(z)): no non-missing arguments to min; returning Inf
  x.minna <- apply(xna, MARGIN = 1, function(z) suppressWarnings(min(which(z))))
  y.minna <- apply(yna, MARGIN = 1, function(z) suppressWarnings(min(which(z))))
  x.mo <- which(x.minna > y.minna)
  y.mo <- which(y.minna > x.minna)
  
  yunder.i <- c()
  yover.i <- c()
  mis.i <- c()
  
  # subset where x has more resolved assignments and only to cols where 
  # both have assignments, then see if names match
  yunder <- 0
  ymis <- 0
  if (length(x.mo) > 0){
    for (i in 1:length(x.mo)){
      if ((y.minna[x.mo[i]]-1) < 3){ # if kingdom is unassigned just add to under-classification
        yunder <- yunder+1
      } else {
        tmp.tmpx <- tmp.x[x.mo[i] , 3:(y.minna[x.mo[i]]-1)]
        tmp.tmpy <- tmp.y[x.mo[i] , 3:(y.minna[x.mo[i]]-1)]
    
        if (all(tmp.tmpx == tmp.tmpy)) {
          yunder <- yunder+1
        } else {
          ymis <- ymis+1
        }
      }
    }
  }
  
  # repeat above where y is more resolved than x:
  yover <- 0
  xmis <- 0
  if (length(y.mo) > 0){
    for (i in 1:length(y.mo)){
      if ((x.minna[y.mo[i]]-1) < 3){ # if kingdom is unassigned just add to under-classification
        yover <- yover+1
      } else {
        tmp.tmpx <- tmp.x[y.mo[i] , 3:(x.minna[y.mo[i]]-1)]
        tmp.tmpy <- tmp.y[y.mo[i] , 3:(x.minna[y.mo[i]]-1)]
    
        if (all(tmp.tmpx == tmp.tmpy)) {
          yover <- yover+1
        } else {
          xmis <- xmis+1
        }
      }
    }
  }
  
  if (yover+xmis != length(y.mo) || yunder+ymis != length(x.mo)) {
    stop("somethings wrong i think")
  }
  
  perf <- nrow(perf)
  # whatever's left is where both tables have the same ranks named, but different names. this is a misclassification:
  moremis <- nrow(x) - (xmis+ymis+yover+yunder+perf)
  result <- data.frame(all.match = c(perf), mis = c(ymis+xmis+moremis), over = c(yover), under = c(yunder))
  if (sum(result) == nrow(x)) {
    return(result)
  } else {
    stop("noooooooooooooo")
  }
}

# make a list of all individual taxonomy tables (not the ensemble):
tbl.list <- list(idtax.pr2.pretty, idtax.silva.mapped2pr2, bayes.pr2.pretty)
all.comp <- lapply(tbl.list, FUN = tblcomper, x = eTax1) # apply above comparison fcn to all idnividual tables
all.comp <- base::do.call(base::rbind.data.frame, all.comp) # clean results
# name the rows (order matches the order you put them into the list)
row.names(all.comp) <- c("idtax-pr2", "idtax-silva", "bayes-pr2")
all.comp <- all.comp / rowSums(all.comp) # convert to proportions

# prep data for plotting:
all.comp$tbl <- row.names(all.comp)
plt.all.comp <- melt(all.comp)

# plot the results:
plt2 <- ggplot(plt.all.comp, aes(fill = tbl, x = variable, y = value)) + 
    geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) + 
    labs(x = "Relative to ensemble", y = "Proportion of ASVs") + 
    scale_x_discrete(breaks = c("all.match","mis","over","under"), 
                   labels = c("Agree (all ranks)", "Disagree (any rank)", "More ranks assigned","Less ranks assigned")) +
    scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                   limits = c(0, 1), expand = c(0,0)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"),
          axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_fill_manual(name = "Taxonomy table",
                      breaks = c("idtax-pr2", "bayes-pr2", "idtax-silva"),
                      values = cbPalette[c(1:3)])
print(plt2)

