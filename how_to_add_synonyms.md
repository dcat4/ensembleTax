Adding a custom synonym file to the ensembleTax workflow
================
D Catlett
8/26/2020

Overview
========

This vignette demonstrates how to use a custom collection of taxonomic synonyms with *taxmapper*. The file in this repo called "tax\_synonyms\_FINAL.csv" contains taxonomic synonyms for eukaryotes compiled by the ensembleTax authors, and these are the default synonyms searched by *taxmapper* in the event that no matches are found when mapping a particular taxonomic name. This list is not comprehensive and contributed additions to this list are welcome.

Example synonym data
--------------------

Here we'll read in the file and how it looks so you can use a custom file with *taxmapper* if you'd like.

``` r
syn <- utils::read.csv("~/Documents/R/amplicon_bioinformatics/taxonomy_pipeline/tax_table_mapping/tax_synonyms_FINAL.csv")
class(syn)
```

    ## [1] "data.frame"

``` r
typeof(syn)
```

    ## [1] "list"

``` r
head(syn)
```

    ##              Name_1          Name_2 Name_3 Name_4 Name_5 Name_6 Name_7
    ## 1       Neoceratium       Ceratium    <NA>   <NA>   <NA>   <NA>   <NA>
    ## 2     Siphonophorae    Siphonophora   <NA>   <NA>   <NA>   <NA>   <NA>
    ## 3    Mamiellales_fa     Mamiellales   <NA>   <NA>   <NA>   <NA>   <NA>
    ## 4        Eurotiales Elaphomycetales   <NA>   <NA>   <NA>   <NA>   <NA>
    ## 5 Gymnodinium_clade     Gymnodinium   <NA>   <NA>   <NA>   <NA>   <NA>
    ## 6 Gymnodiniphycidae   Gymnodinoidia   <NA>   <NA>   <NA>   <NA>   <NA>
    ##                                                                                                                                                  References
    ## 1                                                                          https://img.algaebase.org/pdf/C13170E01dfb92CEDCGWI1A1AD34/gomez_neoceratium.pdf
    ## 2 https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=42759&lvl=3&p=has_linkout&p=blast_url&p=genome_blast&lin=f&keep=1&srchmode=1&unlock
    ## 3                                                                                                                                       dylans common sense
    ## 4  https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5042&lvl=3&p=has_linkout&p=blast_url&p=genome_blast&lin=f&keep=1&srchmode=1&unlock
    ## 5                                                                                                                                       dylans common sense
    ## 6                                                                                                       https://mmbr.asm.org/content/mmbr/57/4/953.full.pdf
    ##   Notes.References X                                                    X.1
    ## 1                                                                          
    ## 2                                                                          
    ## 3                                                                          
    ## 4                                                                          
    ## 5                                                                          
    ## 6                    http://faculty.une.edu/cas/lfritz/old_pages/notes8.htm

...that's it. That's what it is. You can have as many Name\_x and reference/note columns as you want in there. The important part is that the taxonomic synonyms are found in columns with names that start with "Name". When *taxmapper* can't find an exact match for a particular taxonomic name, it finds the row in the *synonym.file* containing that name, and compiles all the other values in that row (excluding NA) that are found in columns with "Name" in their name.
