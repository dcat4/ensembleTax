Create a custom taxonomy to map onto
================
D Catlett
8/26/2020

Overview
========

This vignette demonstrates how to create a custom *tax2map2* object for use with *taxmapper*. A *tax2map2* object is just a collection of all unique taxonomic assignments that can be made using a particular reference database.

Database files used here
------------------------

The files this particular script works with are pulled from the dada2 webpage (<https://benjjneb.github.io/dada2/training.html>), and are standard fasta files with taxonomic names separated by semi-colons. The input files were too big for me to push onto Github, but if you want to follow along go grab them from the dada2 page and just change the file paths below.

RDP example
-----------

We'll do the RDP training set first, starting with reading it in as a DNAStringSet object using Biostrings.

``` r
rm(list = ls())
library("Biostrings")
```

    ## Loading required package: BiocGenerics

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
    ##     union, unique, unsplit, which, which.max, which.min

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
ff <- "~/Documents/R/rdp_train_set_16.fa" # CHANGE to your fasta file name
fastaFile <- Biostrings::readDNAStringSet(ff)
head(fastaFile)
```

    ##   A DNAStringSet instance of length 6
    ##     width seq                                               names               
    ## [1]  1445 GAACGCTGGCGGCGTGCTTAACA...ATCGGCGATTGGGACGAAGTCGT Bacteria;Actinoba...
    ## [2]  1536 GTTTGATCCTGGCTCAGATTGAA...CTAGGGGAACCTGCGGCTGGATC Bacteria;Proteoba...
    ## [3]  1482 TAGAGTTTGATCCTGGCTCAGGA...GAAGTCGTAACAAGGTAGCCGTA Bacteria;Actinoba...
    ## [4]  1386 AAGAGTTTGATCCTGGCTCAGAG...TTGACCTTAANCCGGNGAGCGAA Bacteria;Proteoba...
    ## [5]  1451 ATTGAACGCTGGCGGCATGCCTT...CGGCAGGGTTCGTGACTGGGGTG Bacteria;Proteoba...
    ## [6]  1470 GACGAACGCTGGCGGCGTGCTTA...CAAGGTAGCGTACCGGAAGGTGC Bacteria;Actinoba...

Looks good, now let's work up the sequence names and split them by the ";". We'll get help from the stringr package.

``` r
library("stringr")
seq_name = names(fastaFile)
eh <- stringr::str_split(seq_name, pattern = ";", simplify = TRUE)
eh <- base::as.data.frame(eh[, -ncol(eh)], stringsAsFactors = FALSE)
head(eh)
```

    ##         V1             V2                  V3               V4
    ## 1 Bacteria Actinobacteria      Actinobacteria  Actinomycetales
    ## 2 Bacteria Proteobacteria Gammaproteobacteria      Vibrionales
    ## 3 Bacteria Actinobacteria      Actinobacteria  Actinomycetales
    ## 4 Bacteria Proteobacteria Alphaproteobacteria Rhodospirillales
    ## 5 Bacteria Proteobacteria  Betaproteobacteria  Burkholderiales
    ## 6 Bacteria Actinobacteria      Actinobacteria  Actinomycetales
    ##                               V5            V6
    ## 1               Mycobacteriaceae Mycobacterium
    ## 2                   Vibrionaceae        Vibrio
    ## 3               Mycobacteriaceae Mycobacterium
    ## 4               Acetobacteraceae   Acetobacter
    ## 5 Burkholderiales_incertae_sedis     Aquincola
    ## 6             Micromonosporaceae Couchioplanes

We've successfully created a dataframe with ranks along the columns and rows corresponding to taxonomic assignments. Now we'll just assign rank names, make sure there are no NA's (these will create bugs in *taxmapper*) and make sure the taxonomic assignments are unique (to speed up *taxmapper*).

``` r
r <- c("domain", "phylum", "class", "order", "family", "genus")
colnames(eh) <- r
eh <- base::unique(eh, MARGIN = 1)
any(is.na(eh))
```

    ## [1] FALSE

And one last check of the data frame:

``` r
head(eh)
```

    ##     domain         phylum               class            order
    ## 1 Bacteria Actinobacteria      Actinobacteria  Actinomycetales
    ## 2 Bacteria Proteobacteria Gammaproteobacteria      Vibrionales
    ## 4 Bacteria Proteobacteria Alphaproteobacteria Rhodospirillales
    ## 5 Bacteria Proteobacteria  Betaproteobacteria  Burkholderiales
    ## 6 Bacteria Actinobacteria      Actinobacteria  Actinomycetales
    ## 7 Bacteria Proteobacteria Alphaproteobacteria  Rhodobacterales
    ##                           family          genus
    ## 1               Mycobacteriaceae  Mycobacterium
    ## 2                   Vibrionaceae         Vibrio
    ## 4               Acetobacteraceae    Acetobacter
    ## 5 Burkholderiales_incertae_sedis      Aquincola
    ## 6             Micromonosporaceae  Couchioplanes
    ## 7               Rhodobacteraceae Maritimibacter

``` r
class(eh)
```

    ## [1] "data.frame"

``` r
typeof(eh)
```

    ## [1] "list"

Alright. This data set is now ready to plug into *taxmapper*.

GreenGenes example
------------------

We'll repeat for greengenes:

``` r
ff <- "~/Documents/R/gg_13_8_train_set_97.fa"
fastaFile2 <- Biostrings::readDNAStringSet(ff)
seq_name2 = names(fastaFile2)
eh2 <- stringr::str_split(seq_name2, pattern = ";", simplify = TRUE)
eh2 <- base::as.data.frame(eh2[, -ncol(eh2)], stringsAsFactors = FALSE)
r2 <- c("domain", "phylum", "class", "order", "family", "genus","species")
colnames(eh2) <- r2
eh2 <- base::unique(eh2, MARGIN = 1)

head(eh2)
```

    ##        domain              phylum                  class                order
    ## 1 k__Bacteria p__Gemmatimonadetes              c__Gemm-1                  o__
    ## 2 k__Bacteria    p__Bacteroidetes      c__Flavobacteriia  o__Flavobacteriales
    ## 3 k__Bacteria      p__Chloroflexi         c__Chloroflexi   o__[Roseiflexales]
    ## 4 k__Bacteria   p__Proteobacteria c__Gammaproteobacteria o__Enterobacteriales
    ## 5 k__Bacteria              p__TM6               c__SJA-4                  o__
    ## 6 k__Bacteria             p__OP11              c__OP11-4                  o__
    ##                  family             genus species
    ## 1                   f__               g__     s__
    ## 2  f__Flavobacteriaceae g__Flavobacterium     s__
    ## 3   f__[Roseiflexaceae]    g__Roseiflexus     s__
    ## 4 f__Enterobacteriaceae               g__     s__
    ## 5                   f__               g__     s__
    ## 6                   f__               g__     s__

Ohh this one's different. But the stringr package is so dope\* this is no biggie:

\*I have no affiliation with the creator of the stringr package.

``` r
replace.me.plz <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
for (i in 1:ncol(eh2)) {
  eh2[ , i] <- str_replace_all(eh2[ , i], replace.me.plz[i], "")
}
any(is.na(eh2))
```

    ## [1] FALSE

``` r
head(eh2)
```

    ##     domain           phylum               class             order
    ## 1 Bacteria Gemmatimonadetes              Gemm-1                  
    ## 2 Bacteria    Bacteroidetes      Flavobacteriia  Flavobacteriales
    ## 3 Bacteria      Chloroflexi         Chloroflexi   [Roseiflexales]
    ## 4 Bacteria   Proteobacteria Gammaproteobacteria Enterobacteriales
    ## 5 Bacteria              TM6               SJA-4                  
    ## 6 Bacteria             OP11              OP11-4                  
    ##               family          genus species
    ## 1                                          
    ## 2  Flavobacteriaceae Flavobacterium        
    ## 3   [Roseiflexaceae]    Roseiflexus        
    ## 4 Enterobacteriaceae                       
    ## 5                                          
    ## 6

``` r
class(eh2)
```

    ## [1] "data.frame"

``` r
typeof(eh2)
```

    ## [1] "list"

Now it looks just like the other RDP *tax2map2* object we made above. Yeah!
