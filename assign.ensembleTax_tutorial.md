assign.ensembleTax\_tutorial
================
D Catlett
5/17/2021

# Introduction

This vignette provides detailed examples to demonstrate the functionality of the *assign.ensembleTax* algorithm included with the ensembleTax package. For a more general demonstration of the ensembleTax package functionality/workflow, please go here: <https://github.com/dcat4/ensembleTax>

## The assign.ensembleTax algorithm

*assign.ensembleTax*'s purpose is to synthesize taxonomic assignments made by any number of unique taxonomic assignment methods to determine an "ensemble" assignment for each ASV in a data set. *assign.ensembleTax* requires that each method's taxonomic assignments follow the same taxonomic nomenclature (naming and ranking conventions). If yours don't, use the *taxmapper* algorithm included with the ensembleTax package (vignette here: ADD LINK DYLAN).

By default, *assign.ensembleTax* determines the ensemble taxonomic assignment for each ASV by finding the highest-frequency taxonomic assignment across the input taxonomic assignments (presumably) determined with different methods.

*assign.ensembleTax* includes several user-tuneable parameters that balance obtaining assignments for a larger number of ASVs at lower taxonomic ranks (at the expense of a likely increase in false-positive annotations) with obtaining more robust assignments for fewer ASVs that are supported by multiple methods. Here we'll step through some examples to demonstrate how the algorithm works and how it's behavior can be modified by users.

### Examples

To demonstrate the functionality of *assign.ensembleTax*, we'll first create an artificial set of ASVs and corresponding taxonomic assignments obtained with 2 different artificial "methods". Note these "methods" follow the same naming and ranking conventions.

So first, load the ensembleTax package, and create the artificial data:

``` r
library("ensembleTax")
packageVersion("ensembleTax")
```

    ## [1] '1.1.1'

``` r
# create a fake taxonomy table of ASVs and taxonomic assignments
taxtab1 <- data.frame(ASV = c("sv1", "sv2", "sv3", "sv4"),
                          kingdom = c("Eukaryota", "Eukaryota", "Eukaryota", "Eukaryota"), 
                          supergroup = c("Stramenopile", "Stramenopile", "Alveolata", "Rhizaria"),
                          division = c("Ochrophyta", NA, "Dinoflagellata", NA),
                          class = c("Bacillariophyta", NA, NA, NA),
                          genus = c("Pseudo-nitzschia", NA, NA, NA))
taxtab2 <- data.frame(ASV = c("sv1", "sv2", "sv3", "sv4"),
                          kingdom = c("Eukaryota", "Eukaryota", "Eukaryota", "Eukaryota"), 
                          supergroup = c("Stramenopile", "Alveolata", "Alveolata", "Stramenopile"),
                          division = c("Ochrophyta", "Dinoflagellata", "Dinoflagellata", NA),
                          class = c("Bacillariophyta", NA, "Syndiniales", NA),
                          genus = c("Pseudo-nitzschia", NA, NA, NA))
# look at your artificial data:
taxtab1
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

We see across our 4 ASVs, we have certain ASVs for which the assigned taxonomy is identical across the two tables, or that vary in the number of ranks with assigned names, in the names that were assigned, or both. In what follows we'll see how one might obtain different ensemble assignments based on these data for various scientific questions and/or based on assumptions about the underlying methods employed to obtain each collection of taxonomic assignments.

Now would be a good time to review the *assign.ensembleTax* documentation to get a sense of the different parameter spaces available. Here we'll try to demonstrate what these different parameters are doing.

#### Example 1: Simply obtain the highest frequency assignments

In this example we'll run *assign.ensembleTax* with our two taxonomy tables. *assign.ensembleTax* expects a named list of dataframes with each element corresponding to a uniquely-named taxonomy table, so we'll create that first and then compute ensemble taxonomic assignments with the default parameters.

``` r
xx <- list(taxtab1, taxtab2)
names(xx) <- c("tab1","tab2")
eTax.def <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx)), 
                              tiebreakz = NULL, 
                              count.na=TRUE, 
                              assign.threshold = 0)
# show the initials and ensemble for ease-of-interpretation:
taxtab1
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.def
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota         <NA>           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota         <NA>           <NA>            <NA>             <NA>

Let's break this down. We saw that for sv1, the assignments were in perfect agreement across the two tables and this assignment has been retained in the ensemble. sv2 and sv4 disagreed at the supergroup rank, and so have been left unassigned at the supergroup rank in the ensemble (there is no "highest-frequency" assignment where only two tables are provided and they disagree). For sv3, the assignments were identical down to the division rank, but one table was unassigned (assigned NA) at the class rank while the other was assigned to *Syndiniales*. When *count.NA = TRUE*, NA values are counted as assignments and so again there was no highest frequency assignment at the class rank, resulting in the ensemble being assigned only to division where "Dinoflagellata" was assigned in both input tables.

More on that in a sec. First, we'll add in a third taxonomy table that is identical to taxtab1 and compute ensembles with all 3:

``` r
# create a 3rd fake taxonomy table of ASVs and taxonomic assignments
taxtab3 <- taxtab1
xx.with3 <- list(taxtab1, taxtab2, taxtab3)
names(xx.with3) <- c("tab1", "tab2", "tab3")
eTax.def <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx.with3)), 
                              tiebreakz = NULL, 
                              count.na=TRUE, 
                              assign.threshold = 0)
# show the initials and ensemble for ease-of-interpretation:
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.def
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

With 3 taxonomy tables, we see sv2 and sv4 are now assigned at the supergroup rank in the ensemble. For each of these, both taxtab1 and taxtab3 agreed in their assignments, meaning these assignments were found at a higher frequency than the conflicting supergroup assignments in taxtab2.

#### Example 2: The "count.na" argument

Here we'll create two ensembles again with the same combination of taxonomy tables, but we'll set *count.na = FALSE*. This adjustment is meant for users who want to increase the number of annotated ASVs but likely comes at the expense of an increase in false positive annotations.

``` r
eTax.nona2 <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx)), 
                              tiebreakz = NULL, 
                              count.na=FALSE, 
                              assign.threshold = 0)
# show the initials and ensemble for ease-of-interpretation:
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.nona2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota         <NA>           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota         <NA>           <NA>            <NA>             <NA>

``` r
eTax.nona3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx.with3)), 
                              tiebreakz = NULL, 
                              count.na=FALSE, 
                              assign.threshold = 0)
# ensemble with 3 tables:
eTax.nona3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

Here using the 2 taxonomy tables, we see the ensemble assignments for sv1, 2, and 4 are identical to the first example. Where there are conflicting assignments that are not NA, the *count.na* argument has no impact on ensemble determinations, as you can see in the first ensemble we computed. However, for sv3, we see that we now have an ensemble assignment at the class rank ("Syndiniales") because the NA assignment wasn't counted (in other words, there was 1 Syndiniales and 0 other assignments that were not NA, so the highest frequency assignment was Syndiniales).

You might be somewhat surprised to see that sv3's ensemble class assignment is still Syndiniales when we considered all 3 taxonomy tables. In this case NA was the highest frequency assignment, but by setting *count.na = FALSE* we ignored the NA assignments and so the highest frequency assignment in the absence of NA's was assigned as the ensemble.

#### Example 3: Breaking ties with the "tiebreakz" argument

In this example we'll count NA's again but we'll specify that we'd like to prioritize particular taxonomy tables in the event that multiple disagreeing assignments are found at the highest frequency. Again, using our same 3 tables as above.

``` r
eTax.tb2 <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx)), 
                              tiebreakz = c("tab2"), 
                              count.na=TRUE, 
                              assign.threshold = 0)
# show the initials and ensemble for ease-of-interpretation:
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.tb2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.tb3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx.with3)), 
                              tiebreakz = c("tab2"), 
                              count.na=TRUE, 
                              assign.threshold = 0)
# ensemble with 3 tables:
eTax.tb3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

In both ensembles we prioritized assignments in taxtab2 in the event of ties. When only 2 taxonomy tables were used to compute the ensemble, the ensemble assignments are identical to taxtab2 (the table we chose to break ties with). Makes sense. The ensemble computed with 3 tables was not impacted by tie-breaking because taxtab1 and 3 are identical, and so all assignments in these two tables will be found at the highest frequency and there will be no ties.

One more example computing an ensemble with taxtabs 1 and 2, but prioritizing taxtab1 to break ties and NOT counting NA's:

``` r
eTax.tb2 <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx)), 
                              tiebreakz = c("tab1"), 
                              count.na=FALSE, 
                              assign.threshold = 0)
# show the initials and ensemble for ease-of-interpretation:
taxtab1 
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.tb2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

And now we see where tie-breaking can actually be useful. If you'd like to primarily rely on assignments from one assignment method, but use a second method to fill in annotations where your favorite method is not assigned (assigned NA), you can specify a tiebreaker and ignore NA's in your ensemble assignment determinations.

#### Example 4: Prioritizing methods with the "weights" argument

Another way you can prioritize assignments from one (or more) particular assignment method is by using the weights argument. Weights are specified as integers in the order corresponding to the order of taxonomy tables in the list you supply to *assign.ensembleTax*. Weighting one table more highly than the other in ensembles determined from only two taxonomy tables will result in identical behavior as tie-breaking. Below we've weighted assignments in taxtab1 double those found in taxtab2:

``` r
# counting NA's:
eTax.wt2 <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=c(2,1), 
                              tiebreakz = NULL, 
                              count.na=TRUE, 
                              assign.threshold = 0)
# show the initials and ensemble for ease-of-interpretation:
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.wt2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
# NOT counting NA's:
eTax.wt2 <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=c(2,1), 
                              tiebreakz = NULL, 
                              count.na=FALSE, 
                              assign.threshold = 0)
# show the initials and ensemble for ease-of-interpretation:
eTax.wt2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

Just as expected. The ensemble assignments are identical to taxtab1 when we count NA's, but when we ignore NA's the Syndiniales assignment is filled in by taxtab2.

What happens when we compute a 3-table ensemble but weight the table that disagrees with the other two (taxtab2) 2x? Let's see:

``` r
eTax.wt3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=c(1,2,1), 
                              tiebreakz = NULL, 
                              count.na=TRUE, 
                              assign.threshold = 0)
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.wt3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota         <NA>           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota         <NA>           <NA>            <NA>             <NA>

We see that anywhere where taxonomic assignments were in disagreement, they are not assigned in the ensemble. This is because taxtab 1 and 3 are identical to one another, and we've weighted taxtab2 double. So where assignments disagree there are multiple assignments with the highest frequency. That means we need to specify a tiebreaker if we want to avoid the above scenario. Let's try:

``` r
eTax.wttb3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=c(1,2,1), 
                              tiebreakz = c("tab1"), 
                              count.na=TRUE, 
                              assign.threshold = 0)
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.wttb3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

And now our ensemble determinations match taxtab1 (and 3) where there were disagreements with taxtab2. We can prioritize taxtab2 and see the opposite:

``` r
eTax.wttb3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=c(1,2,1), 
                              tiebreakz = c("tab2"), 
                              count.na=TRUE, 
                              assign.threshold = 0)
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.wttb3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

Just as we expected.

#### Example 5: the "assign.threshold" argument

We have one last argument to address: *assign.threshold*. This argument interacts with some of the other arguments we've looked at above in ways that may be counter-intuitive to some, so make sure you understand what each of those arguments are doing before spending too much time with *assign.threshold*.

Let's take a look at some different thresholds alongside our default parameters first. We'll see how changing *assign.threshold* can impact tiebreaking and weighting first:

``` r
# tie-breaking to prioritize table 1, but with assign.threshold = 60%
eTax.at <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx)), 
                              tiebreakz = c("tab1"), 
                              count.na=TRUE, 
                              assign.threshold = 0.6)
# show the initials and ensemble for ease-of-interpretation:
taxtab1 # (remember taxtab3 is identical to this, so count 2x)
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
taxtab2
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota Stramenopile           <NA>            <NA>             <NA>

``` r
eTax.at
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota         <NA>           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota         <NA>           <NA>            <NA>             <NA>

``` r
# take away the tiebreaker and weight table 1 2x:
eTax.at <- assign.ensembleTax(xx, 
                              tablenames = names(xx), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=c(2,1), 
                              tiebreakz = NULL, 
                              count.na=TRUE, 
                              assign.threshold = 0.6)
eTax.at
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

Here we see that the *assign.threshold* argument "over-rules" tie-breaking if the threshold is not satisfied. Although we specified taxtab1 as the tie-breaker, by designating *assign.threshold = 0.6*, we've required ensemble assignments to be found in at least 60% of the (weighted) assignments. This means disagreements were left unassigned here.

If we take a look at the 2nd ensemble we computed (where tie-breaking was omitted and we instead weighted taxtab1 2x), we see that our ensemble mirrors taxtab1. This is because the *assign.threshold* argument operates on weighted assignment frequencies. In this case, where taxtab 1 and 2 disagreed, the assignments in taxtab1 comprised 66% of the weighted assignments, which was larger than our *assign.threshold*.

Let's try a couple thresholds applied to 3-table ensemble determinations:

``` r
# a low threshold:
eTax.at3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx.with3)), 
                              tiebreakz = NULL, 
                              count.na=TRUE, 
                              assign.threshold = 0.5)
eTax.at3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
# a high threshold (need all 3 to agree here for ensemble assignment):
eTax.at3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx.with3)), 
                              tiebreakz = NULL, 
                              count.na=TRUE, 
                              assign.threshold = 0.9)
eTax.at3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota         <NA>           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata            <NA>             <NA>
    ## 4 sv4 Eukaryota         <NA>           <NA>            <NA>             <NA>

In our first example, the threshold we applied (0.5) has no impact on the ensemble assignments because taxtab1 and 3 are identical, and therefore their assignments will always comprise more than 50% of the assignments for any ASV.

In our second example, the threshold of 0.9 (and in this case, any threshold &gt; 0.67) means that the ensemble is only assigned where all three input taxonomy tables are in agreement.

Finally, we'll take a look at how *assign.threshold* behaves when *count.na = FALSE*:

``` r
# a low threshold with count.na = FALSE:
eTax.at3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx.with3)), 
                              tiebreakz = NULL, 
                              count.na=FALSE, 
                              assign.threshold = 0.5)
eTax.at3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota Stramenopile           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota     Rhizaria           <NA>            <NA>             <NA>

``` r
# a high threshold with count.na = FALSE (need all 3 to agree here for ensemble assignment):
eTax.at3 <- assign.ensembleTax(xx.with3, 
                              tablenames = names(xx.with3), 
                              ranknames = colnames(taxtab1)[2:ncol(taxtab1)],
                              weights=rep(1,length(xx.with3)), 
                              tiebreakz = NULL, 
                              count.na=FALSE, 
                              assign.threshold = 0.9)
eTax.at3
```

    ##   ASV   kingdom   supergroup       division           class            genus
    ## 1 sv1 Eukaryota Stramenopile     Ochrophyta Bacillariophyta Pseudo-nitzschia
    ## 2 sv2 Eukaryota         <NA>           <NA>            <NA>             <NA>
    ## 3 sv3 Eukaryota    Alveolata Dinoflagellata     Syndiniales             <NA>
    ## 4 sv4 Eukaryota         <NA>           <NA>            <NA>             <NA>

We again see that when the threshold is low (0.5), it has little impact on our ensemble assignments.

However, implementing a high threshold (0.9) in conjunction with *count.na = FALSE* notably alters the ensemble assignments. For sv2 and sv4, where taxonomic names (and not NA's) were assigned and disagreed across the 3 tables, the ensemble remains unassigned because no single assignment was found at a frequency greater than 90%. However, for sv3, where 66% of the input class assignments were NA and 33% were Syndiniales, the ensemble is assigned to the class Syndiniales. This is because the *assign.threshold* does not consider NA assignments when *count.na = FALSE*. In other words, because NA's were ignored, in this example the Syndiniales assignment comprised 100% of the input assignments and thus surpassed the *assign.threshold*.

### Parameter summary and recommendations for different scientific objectives

As you can see in the above, the *assign.ensembleTax* algorithm allows for flexible computations of ensemble taxonomic assignments to suit different scientific questions and applications. There are trade-offs in adjusting each of the above parameters and here we'll discuss what those are and how you might implement *assign.ensembleTax* for your own objectives.

The big trade-off one should consider when implementing *assign.ensembleTax* is how to balance annotating more ASVs at lower ranks vs. only annotating ASVs where annotations are supported by multiple methods (and thus are likely to be quite robust). The former likely comes at the expense of increased false positive annotations (ASVs assigned to a lineage where there really is not enough information to make this determination), while the latter likely comes at the expense of increased false negative annotations (assigning NA where there IS enough information to assign the ASV to a lineage). Adjusting the parameters in
*assign.ensembleTax* allows you to decide where you fall on this spectrum.

First, we'll consider parameters that will promote robust annotations where ASVs are assigned to a taxonomic group at the expense of assigning taxonomy to a smaller number of ASVs/ranks. Arguments that would favor this strategy are: *count.na = TRUE* *assign.threshold = \[a high value, say &gt; 0.5\]* *tiebreakz = NULL* One could use all of the above settings, or select a few. Broadly, all of these parameters require a taxonomic assignment for a particular ASV to be found in multiple input taxonomy tables in order for it to be assigned to the ensemble. Implementing all of these simultaneously will result in extremely conservative, but extremely well-supported taxonomic assignments. These settings are probably most appropriate for studies that require precise identification of particular ASVs. While lower-rank assignments should in general be interpreted with caution, those supported by multiple methods are likely less error-prone than those determined with a single method.

Let's say you think (or know), that assignments made by one particular method are more accurate than the other methods you've considered. Perhaps the reference database uses a more robust annotation procedure, or benchmarking exercises show a particular classifier has very low error rates. You can prioritize assignments made by this method with the *tiebreakz* or *weights* argument. If you wanted to use secondary methods to "fill in" assignments where your preferred method determined that taxonomy could not be assigned to an ASV, you could subsequently change the *count.na* argument. These changes bring you closer to the other end of the spectrum...

At the other end of the spectrum, we'll consider parameters that promote obtaining annotations at lower ranks for a greater number of ASVs at the expense of potentially increasing the number of spurious classifications (either assigning ASVs to an incorrect lineage, or assigning an ASV to a lineage when there is not sufficient phylogenetic resolution to do so). Arguments that would favor this strategy are: *count.na = FALSE* *assign.threshold = \[a low value, say 0\]* *tiebreakz = \[specify the names of ALL input tables in order of priority\]* Broadly, all of these parameters favor an increase in the number of ASVs assigned at lower ranks, but these assignments may only be supported by a single method and so may be less robust. These settings are probably most appropriate for studies focused on very broad taxonomic groupings; lower-rank assignments should in general be interpreted with caution, and this is particularly true when ensemble assignments are computed with these parameters.

That brings us to the end of this tutorial. We include some code for comparing ensemble taxonomic assignments in our ensembleTax package overview vignette, and encourage you to perform such comparisons to test out different settings and optimize ensemble assignments for your particular scientific objectives.
