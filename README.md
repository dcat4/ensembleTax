The ensembleTax R package
================
D Catlett
5/15/2021

# ensembleTax overview

ensembleTax is an R package that allows incorporation of information from multiple taxonomic assignment algorithms and/or reference databases to compute ensemble taxonomic assignments for ASVs/OTUs generated by common marker gene sequence analyses. The package was built to conveniently compliment the dada2 R package and pipeline, but can be used with any combination of taxonomic assignment algorithms and reference databases if the data can be read into R and converted to a dataframe.

Here you can find instructions for downloading and installing the ensembleTax package, and a brief demonstration of a sample ensembleTax "pipeline".

Please ask questions and/or report bugs on the issue-tracker associated with this repository. Contributions of R-able wrappers or implementations of additional taxonomic assignment algorithms and/or reference database taxonomic nomenclatures are welcome.

## The problem(s) addressed by ensembleTax

Taxonomic assignment of marker gene sequences is a critical step of marker gene workflows as it imparts ecological significance and understanding to genetic data.

Many taxonomic assignment algorithms have been proposed to assign taxonomy to marker gene sequences (or OTUs/ASVs). Similarly, analysts are often forced to choose from one of several reference databases containing representative marker gene sequences with known taxonomic identities. The "best" assignment algorithm and/or reference database for a particular scientific question is often not obvious. To complicate things further, different reference databases generally do not share consistent taxonomic naming or ranking conventions.

ensembleTax solves this problem by providing flexible algorithms that synthesize information from multiple taxonomic assignment algorithm/reference database combinations and compute a single ensemble taxonomic assignment for each ASV/OTU in a marker gene data set.

## Quick reference:

This vignette gives a broad overview of the ensembleTax workflow.

Vignettes for more details on the core ensembleTax algorithms are linked below: taxmapper: <https://github.com/dcat4/ensembleTax/blob/master/taxmapper_tutorial.md> assign.ensembleTax: <https://github.com/dcat4/ensembleTax/blob/master/assign.ensembleTax_tutorial.md>

## Download instructions

The ensembleTax package is available on Github and CRAN. To install from Github, first install ensembleTax's dependencies from CRAN and Bioconductor (and install Bioconductor if you don't have it), then use devtools to install from Github as follows:

``` r
install.packages(c("dplyr", "stringr", "usethis", "devtools"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DECIPHER", "Biostrings", "dada2"))

library(devtools)
devtools::install_github("dcat4/ensembleTax", build_manual = TRUE, build_vignettes = TRUE)
packageVersion("ensembleTax")
```

Alternatively, you can install from CRAN by doing:

``` r
install.packages(c("ensembleTax"))
packageVersion("ensembleTax")
```

## ensembleTax algorithms

The core algorithms employed by ensembleTax are *taxmapper* and *assign.ensembleTax*. *taxmapper* maps, or 'translates', one taxonomic nomenclature onto another by exact name matching. *taxmapper* is rank-agnostic, meaning it does not consider the hierarchical structure of a taxonomy and assumes that a taxonomic name means the same thing regardless of which reference database employs it (there are times when this assumption is not valid, for example when "Clade\_X" or similar is used as a stand-alone taxonomic name, but *taxmapper* can handle these cases without introducing errors in taxonomic assignments).

*assign.ensembleTax* computes ensemble taxonomic assignments based on assignments determined by any number of individual taxonomic assignment algorithm/reference database combinations. *assign.ensembleTax* includes several user-tuneable parameters that balance obtaining assignments at lower taxonomic ranks for a larger number of ASVs (at the expense of a likely increase in false-positive annotations) with obtaining more robust assignments for fewer ASVs that are supported by multiple methods.

Additional functions are included for pre-processing taxonomic assignments generated by specific taxonomic assignment algorithms and reference databases. These functions are designed to conveniently plug in downstream of the dada2 pipeline, but other pipelines may be used if the data is formatted properly for use with *taxmapper* and/or *assign.ensembleTax* (see the documentation for the objects these algorithms can work with).

The outputs of the following taxonomic assignment algorithms are explicitly supported by ensembleTax:

1.  RDP bayesian classifier as implemented in dada2's assignTaxonomy.
2.  idtaxa algorithm as implemented in DECIPHER.

Supported reference databases include:

1.  Silva SSU NR reference database v138 (silva).
2.  Protistan Ribosomal Reference database v4.12.0 (pr2).
3.  RDP train set v16
4.  GreenGenes v13.8 clustered at 97% similarity

Note that other databases may still be used with ensembleTax, but they must be mapped onto the taxonomic nomenclatures employed by one of the above supported databases using *taxmapper*, or the user must extract all unique taxonomic assignments from the database and format them properly for use with *taxmapper*. A vignette that shows how to extract taxonomic assignments from a database not supported by the ensembleTax package can be found here: <https://github.com/dcat4/ensembleTax/blob/master/add_tax_dbs.md>

A vignette that shows how to incorporate taxonomic assignments produced outside of R into an ensembleTax workflow can be found here: <https://github.com/dcat4/ensembleTax/blob/master/add_tax_tab_from_csv.md>

### ensembleTax 'pipeline' demonstration

Here we step through a simple example of an ensembleTax workflow to compute ensemble taxonomic assignments for a small set of 18S-V9 protist ASVs.

Since ensembleTax is built to plug in directly downstream of the dada2 pipeline, to start the ensembleTax pipeline demonstration we include here code that can be used to generate initial sets of taxonomic assignments using the RDP Bayesian classifier implemented in dada2's assignTaxonomy and DECIPHER's idtaxa algorithm. The below snippet comes directly from the dada2 tutorial. Please see the dada2 pipeline tutorial for further information: <https://benjjneb.github.io/dada2/tutorial.html>

``` r
# dada2 assignTaxonomy:
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# this is optional and was not used in the example data below:
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz") 

# DECIPHER idtaxa:
library(DECIPHER); packageVersion("DECIPHER")
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
```

We encourage creating an ASV "rubric" to track ASVs through the ensembleTax pipeline. The rubric is a DNAStringSet (see the Biostrings package) object produced by extracting ASV sequences from the seqtab output by dada2, and giving them arbitrary names (like sv1, sv2, etc). You can create a "rubric" from an "ASV table" output by dada2 by doing the following:

``` r
rubric <- DNAStringSet(getSequences(seqtab.nochim))
# this creates names (sv1, sv2, ..., svX) for each ASV
snam <- vector(mode = "character", length = length(rubric))
for (i in 1:length(rubric)) {
    snam[i] <- paste0("sv", as.character(i))
}
names(rubric) <- snam
```

For this vignette, instead of running the assignments we'll just load some data included with the ensembleTax package. These are outputs of dada2's assignTaxonomy implemented against pr2, and of DECIPHER's idtaxa implemented against both pr2 and silva. The rubric.sample is an example of a "rubric", which ensembleTax uses to track ASV-identifying information.

``` r
library("ensembleTax")
library("Biostrings")
```

    ## Loading required package: BiocGenerics

    ## Warning: package 'BiocGenerics' was built under R version 4.0.5

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
data("idtax.pr2.sample")
data("idtax.silva.sample")
data("bayes.sample")
data("rubric.sample")
```

If you want to see what they look like, do this:

``` r
head(idtax.pr2.sample)
head(idtax.silva.sample)
head(bayes.sample)
head(rubric.sample)
```

#### ensembleTax pre-processing

We see from the above that the data structures returned by our two taxonomic assignment algorithms are different. It is critically important that the order of sequences in the rubric and in the idtaxa-returned Taxon object are the same. idtaxa does not return sequence names, but if you did not alter the order of your sequences as you provided them to idtaxa and/or DNAStringSet when creating your rubric, the ordering should be preserved and you should be good to go.

Here we'll run these tables through ensembleTax's pre-processing functions. Supplying a rubric allows ensembleTax to give each taxonomy table the same ASV-identifying information and to better track and organize your data.

``` r
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
```

Again to view the data, do this:

``` r
# remove ASV columns from tax tables for easier viewing:
idtax.pr2.pretty <- idtax.pr2.pretty[ , -which(names(idtax.pr2.pretty) %in% c("ASV"))]
idtax.silva.pretty <- idtax.silva.pretty[ , -which(names(idtax.silva.pretty) %in% c("ASV"))]
bayes.pr2.pretty <- bayes.pr2.pretty[ , -which(names(bayes.pr2.pretty) %in% c("ASV"))]

head(idtax.pr2.pretty)
```

    ##       svN   kingdom    supergroup   division           class             order
    ## 1 sv14136 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 2 sv17278 Eukaryota Stramenopiles Ochrophyta            <NA>              <NA>
    ## 3 sv20747 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 4  sv3579      <NA>          <NA>       <NA>            <NA>              <NA>
    ## 5  sv4298 Eukaryota Stramenopiles       <NA>            <NA>              <NA>
    ##                       family       genus species
    ## 1 Polar-centric-Mediophyceae        <NA>    <NA>
    ## 2                       <NA>        <NA>    <NA>
    ## 3 Polar-centric-Mediophyceae Chaetoceros    <NA>
    ## 4                       <NA>        <NA>    <NA>
    ## 5                       <NA>        <NA>    <NA>

``` r
head(idtax.silva.pretty)
```

    ##       svN    domain       phylum   class order family genus
    ## 1 sv14136 Eukaryota         <NA>    <NA>  <NA>   <NA>  <NA>
    ## 2 sv17278 Eukaryota Eukaryota_ph  MOCH-2  <NA>   <NA>  <NA>
    ## 3 sv20747 Eukaryota         <NA>    <NA>  <NA>   <NA>  <NA>
    ## 4  sv3579 Eukaryota         <NA>    <NA>  <NA>   <NA>  <NA>
    ## 5  sv4298 Eukaryota       MAST-3 MAST-3B  <NA>   <NA>  <NA>

``` r
head(bayes.pr2.pretty)
```

    ##       svN   kingdom     supergroup     division           class
    ## 1 sv14136 Eukaryota  Stramenopiles   Ochrophyta Bacillariophyta
    ## 2 sv17278 Eukaryota  Stramenopiles   Ochrophyta          MOCH-2
    ## 3 sv20747 Eukaryota  Stramenopiles   Ochrophyta Bacillariophyta
    ## 4  sv3579 Eukaryota Archaeplastida Streptophyta   Embryophyceae
    ## 5  sv4298 Eukaryota  Stramenopiles     Opalozoa          MAST-3
    ##               order                     family       genus        species
    ## 1 Bacillariophyta_X Polar-centric-Mediophyceae  Minidiscus Minidiscus_sp.
    ## 2          MOCH-2_X                  MOCH-2_XX  MOCH-2_XXX MOCH-2_XXX_sp.
    ## 3 Bacillariophyta_X Polar-centric-Mediophyceae Chaetoceros           <NA>
    ## 4   Embryophyceae_X           Embryophyceae_XX        <NA>           <NA>
    ## 5           MAST-3B                  MAST-3B_X  MAST-3B_XX MAST-3B_XX_sp.

We see that each taxonomy table is now a dataframe sorted by the column "svN".

#### The taxmapper algorithm

After pre-processing our taxonomic assignment data sets above, we see we still can't make apples-to-apples comparisons between the "idtax-silva" table and the other two because they employ different ranking and (though this may not be as obvious) naming conventions. *taxmapper* was created to solve this problem.

Here we'll use *taxmapper* to 'translate' the idtax-silva taxonomic assignments onto the same taxonomic nomenclature as the other two tables.

For more detailed examples illustrating the behavior of *taxmapper* with different input arguments, there's a vignette here: <https://github.com/dcat4/ensembleTax/blob/master/taxmapper_tutorial.md>

Note that if you want to use your own custom synonym.file, there's a vignette showing how to do this here: <https://github.com/dcat4/ensembleTax/blob/master/how_to_add_synonyms.md>

``` r
idtax.silva.mapped2pr2 <- taxmapper(idtax.silva.pretty,
                      tt.ranks = colnames(idtax.silva.pretty)[2:ncol(idtax.silva.pretty)],
                      tax2map2 = "pr2",
                      exceptions = c("Archaea", "Bacteria"),
                      ignore.format = TRUE,
                      synonym.file = "default",
                      streamline = TRUE,
                      outfilez = NULL)
head(idtax.silva.mapped2pr2)
```

    ##       svN   kingdom    supergroup   division  class   order family genus
    ## 1 sv14136 Eukaryota          <NA>       <NA>   <NA>    <NA>   <NA>  <NA>
    ## 2 sv17278 Eukaryota Stramenopiles Ochrophyta MOCH-2    <NA>   <NA>  <NA>
    ## 3 sv20747 Eukaryota          <NA>       <NA>   <NA>    <NA>   <NA>  <NA>
    ## 4  sv3579 Eukaryota          <NA>       <NA>   <NA>    <NA>   <NA>  <NA>
    ## 5  sv4298 Eukaryota Stramenopiles   Opalozoa MAST-3 MAST-3B   <NA>  <NA>
    ##   species
    ## 1    <NA>
    ## 2    <NA>
    ## 3    <NA>
    ## 4    <NA>
    ## 5    <NA>

Inspection of the mapped taxonomy table shows that it now mirrors the naming and ranking conventions of the other two taxonomy tables.

#### The assign.ensembleTax algorithm

Now we have three different taxonomy tables with independent taxonomic assignments for each ASV in our example data set. From these we can compute ensemble taxonomic assignments with the *assign.ensembleTax* algorithm.

Here's a run with the default parameters:

``` r
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
```

    ## $`idtax-pr2`
    ##       svN   kingdom    supergroup   division           class             order
    ## 1 sv14136 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 2 sv17278 Eukaryota Stramenopiles Ochrophyta            <NA>              <NA>
    ## 3 sv20747 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 4  sv3579      <NA>          <NA>       <NA>            <NA>              <NA>
    ## 5  sv4298 Eukaryota Stramenopiles       <NA>            <NA>              <NA>
    ##                       family       genus species
    ## 1 Polar-centric-Mediophyceae        <NA>    <NA>
    ## 2                       <NA>        <NA>    <NA>
    ## 3 Polar-centric-Mediophyceae Chaetoceros    <NA>
    ## 4                       <NA>        <NA>    <NA>
    ## 5                       <NA>        <NA>    <NA>
    ## 
    ## $`idtax-silva`
    ##       svN   kingdom    supergroup   division  class   order family genus
    ## 1 sv14136 Eukaryota          <NA>       <NA>   <NA>    <NA>   <NA>  <NA>
    ## 2 sv17278 Eukaryota Stramenopiles Ochrophyta MOCH-2    <NA>   <NA>  <NA>
    ## 3 sv20747 Eukaryota          <NA>       <NA>   <NA>    <NA>   <NA>  <NA>
    ## 4  sv3579 Eukaryota          <NA>       <NA>   <NA>    <NA>   <NA>  <NA>
    ## 5  sv4298 Eukaryota Stramenopiles   Opalozoa MAST-3 MAST-3B   <NA>  <NA>
    ##   species
    ## 1    <NA>
    ## 2    <NA>
    ## 3    <NA>
    ## 4    <NA>
    ## 5    <NA>
    ## 
    ## $`bayes-pr2`
    ##       svN   kingdom     supergroup     division           class
    ## 1 sv14136 Eukaryota  Stramenopiles   Ochrophyta Bacillariophyta
    ## 2 sv17278 Eukaryota  Stramenopiles   Ochrophyta          MOCH-2
    ## 3 sv20747 Eukaryota  Stramenopiles   Ochrophyta Bacillariophyta
    ## 4  sv3579 Eukaryota Archaeplastida Streptophyta   Embryophyceae
    ## 5  sv4298 Eukaryota  Stramenopiles     Opalozoa          MAST-3
    ##               order                     family       genus        species
    ## 1 Bacillariophyta_X Polar-centric-Mediophyceae  Minidiscus Minidiscus_sp.
    ## 2          MOCH-2_X                  MOCH-2_XX  MOCH-2_XXX MOCH-2_XXX_sp.
    ## 3 Bacillariophyta_X Polar-centric-Mediophyceae Chaetoceros           <NA>
    ## 4   Embryophyceae_X           Embryophyceae_XX        <NA>           <NA>
    ## 5           MAST-3B                  MAST-3B_X  MAST-3B_XX MAST-3B_XX_sp.

``` r
head(eTax1)
```

    ##       svN   kingdom    supergroup   division           class             order
    ## 1 sv14136 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 2 sv17278 Eukaryota Stramenopiles Ochrophyta          MOCH-2              <NA>
    ## 3 sv20747 Eukaryota Stramenopiles Ochrophyta Bacillariophyta Bacillariophyta_X
    ## 4  sv3579 Eukaryota          <NA>       <NA>            <NA>              <NA>
    ## 5  sv4298 Eukaryota Stramenopiles   Opalozoa          MAST-3           MAST-3B
    ##                       family       genus species
    ## 1 Polar-centric-Mediophyceae        <NA>    <NA>
    ## 2                       <NA>        <NA>    <NA>
    ## 3 Polar-centric-Mediophyceae Chaetoceros    <NA>
    ## 4                       <NA>        <NA>    <NA>
    ## 5                       <NA>        <NA>    <NA>

We see that the assignments made at the highest frequency across the 3 assignment algorithms are assigned as the ensemble taxonomic assignment.

For more detailed examples illustrating the behavior of *assign.ensembleTax* with different input arguments, there's a vignette here: <https://github.com/dcat4/ensembleTax/blob/master/assign.ensembleTax_tutorial.md>

#### Comparisons of ensemble assignments with individual methods

As you can see above and in the *assign.ensembleTax* vignette, ensemble taxonomic assignments can be computed in different ways depending on your scientific objectives. We encourage you to try a few different settings that you think are appropriate for your science, and compare your ensemble assignments with those predicted by individual methods.

We demonstrate some possibilities for doing so below. We'll use the reshape2 and ggplot2 packages to help with this, so install those if you don't have them.

First we compare the proportion of ASVs that are unassigned (assigned NA) at each rank:

``` r
install.packages(c("ggplot2", "reshape2"))
```

``` r
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
```

![](README_files/figure-markdown_github/unnamed-chunk-12-1.png)

For our simple 5-ASV example in this vignette, you could more or less see these patterns just by manually inspecting the data. For thousands of ASVs though, this type of analysis can be useful.

A couple things to note here: in most cases, if you've mapped assignments from one database onto another, some information will likely be lost because reference databases occasionally do not contain the same taxonomic names and/or organisms. Therefore you shouldn't read too much into it if assignments made with one reference database consistently feature a higher proportion of ASVs that remain unassigned. In this case, the idtax-silva was mapped onto the pr2 nomenclature, so we don't really know whether the lower proportion of assigned ASVs is due to lower coverage of protists in silva, or due to a loss of information during mapping. We could of course map the pr2 tables onto silva and compare the information lost in both directions to give a better idea if there are systematic differences in assignments made with silva vs. pr2, but this is for another day.

We can however see that the when applied to the same database, the RDP classifier classified more of our ASVs than the idtaxa classifier. This was noted in the idtaxa paper and may be due to high overclassification error rates associated with the RDP classifier. Again though, we don't want to get too much into the weeds here.

Finally, we see that the ensemble is somewhere in the middle of the individual methods. This is because it has only assigned taxonomy where an assignment was supported by two assignment methods; if those were overclassification errors in the bayes-pr2 table, we've reduced the impact of those but it looks like we also had enough support to make taxonomy predictions in some places where the idtax-pr2 table was unassigned (and perhaps was too conservative?). This hints at some of the expected benefits of determining ensemble assignments, though further work is needed to validate these expectations.

We'll demonstrate one more comparison you might try. Here, we're comparing taxonomic annotations for each ASV in our data set according to each independent method with each ASV's ensemble taxonomic assignment. In these comparisons, we'll designate four possible outcomes for each ASV's taxonomic assignments: the assignments predicted by the individual method can be in perfect agreement with the ensemble at all ranks; the individual method can assign taxonomy for an ASV at more or less ranks than the ensemble; or the taxonomy predicted by the individual method can disagree with the ensemble assignment at any rank. Let's do it:

``` r
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
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

And we have a few things to point out here. First, note that for each taxonomy table (bar color), the bars sum to 1. You could thus employ a stacked bar or similar, but I prefer this format; see the ggplot2 help for tweaking the plot. For interpretation's sake, this just means that each ASV's assignment comparison must fall into one of the four categories we designated.

We also see that in this very small data set, we had no ASVs for which different taxonomic assignments were predicted by the individual methods. There were of course instances where one method assigned taxonomy to more or fewer ranks, but this lack of disagreement is a good thing. Rates of conflicting assignments tend to be low (&lt; 5-10%) in the larger data sets I've worked with too. This of course probably depends on the classifier and reference database employed, but is a good sign nonetheless.

One last thing to note: under certain implementations of *assign.ensembleTax*, assignments predicted by one or more methods can be intentionally prioritized by the user. Therefore, it is important to keep in mind that increased disagreements between an individual method and an ensemble does not necessarily indicate that the individual method is error prone.

While this small data set doesn't lend itself to many more interesting interpretations, hopefully the above comparisons help you get started in assessing the utility and optimal parameters for your ensemble taxonomic assignments.

This brings us to the end of our vignette. Please see the vignettes specifically geared toward demonstrating uses for the *taxmapper* and *assign.ensembleTax* algorithms for more information, and let us know on the issues page if/when you find issues!
