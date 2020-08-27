Add non-supported taxonomic assignments to ensembleTax workflow
================
D Catlett
8/26/2020

Overview
========

This vignette demonstrates how to import a collection of taxonomic assignments from a csv file and use these in the ensembleTax package. This might be required if you want to use a taxonomic assignment method that's not explicitly supported by the ensembleTax package.

Data pre-processing
-------------------

The assignments we'll use come from an 18S-V9 rDNA data set and were assigned taxonomy by the lowest common ancestor algorithm in MEGAN 6. The csv files are in this repo if you want to follow along. There's one from running MEGAN-LCA against pr2 and another from running against Silva.

Start by reading in and looking at the data:

``` r
rm(list = ls())
lca.pr2 <- read.csv("tmp_MEGAN_LCA_pr2.csv", stringsAsFactors = FALSE, header = FALSE)
lca.silva <- read.csv("tmp_MEGAN_LCA_Silva.csv", stringsAsFactors = FALSE, header = FALSE)

head(lca.pr2)
```

    ##        V1
    ## 1 sv22340
    ## 2 sv21398
    ## 3 sv22912
    ## 4 sv23338
    ## 5 sv16360
    ## 6 sv12155
    ##                                                                                                                                                       V2
    ## 1                                        NCBI;cellular organisms;Eukaryota;Alveolata;Apicomplexa;Conoidasida;Gregarinasina;Eugregarinorida;Porosporidae;
    ## 2                             NCBI;cellular organisms;Eukaryota;Haptista;Haptophyta;Prymnesiophyceae;Prymnesiales;Chrysochromulinaceae;Chrysochromulina;
    ## 3                                                                              NCBI;cellular organisms;Eukaryota;Katablepharidophyta;Katablepharidaceae;
    ## 4                                                                            NCBI;cellular organisms;Eukaryota;Opisthokonta;Choanoflagellata;Craspedida;
    ## 5                                                                              NCBI;cellular organisms;Eukaryota;Rhizaria;Cercozoa;Chlorarachniophyceae;
    ## 6 NCBI;cellular organisms;Eukaryota;Stramenopiles;Bacillariophyta;Coscinodiscophyceae;Chaetocerotophycidae;Chaetocerotales;Chaetocerotaceae;Chaetoceros;

``` r
head(lca.silva)
```

    ##        V1
    ## 1 sv22340
    ## 2 sv21398
    ## 3 sv22912
    ## 4 sv23338
    ## 5 sv16360
    ## 6 sv12155
    ##                                                                                     V2
    ## 1                                                                        NCBI;No hits;
    ## 2 NCBI;cellular organisms;Eukaryota;Haptista;Haptophyta;Prymnesiophyceae;Prymnesiales;
    ## 3                 NCBI;cellular organisms;Eukaryota;Alveolata;Dinophyceae;Syndiniales;
    ## 4           NCBI;cellular organisms;Eukaryota;Opisthokonta;Fungi;Fungi incertae sedis;
    ## 5          NCBI;cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;
    ## 6                                                                        NCBI;No hits;

We've got some arbitrary names for our ASVs, and taxonomic assignments separated by semi-colons. We'll use the stringr package to separate the taxonomic names into columns of a matrix, with each column corresponding to a rank:

``` r
library("stringr")
new.pr2 <- cbind(lca.pr2$V1, stringr::str_split(lca.pr2$V2, pattern = ";", simplify = TRUE))
new.pr2
```

    ##      [,1]      [,2]   [,3]                 [,4]        [,5]                 
    ## [1,] "sv22340" "NCBI" "cellular organisms" "Eukaryota" "Alveolata"          
    ## [2,] "sv21398" "NCBI" "cellular organisms" "Eukaryota" "Haptista"           
    ## [3,] "sv22912" "NCBI" "cellular organisms" "Eukaryota" "Katablepharidophyta"
    ## [4,] "sv23338" "NCBI" "cellular organisms" "Eukaryota" "Opisthokonta"       
    ## [5,] "sv16360" "NCBI" "cellular organisms" "Eukaryota" "Rhizaria"           
    ## [6,] "sv12155" "NCBI" "cellular organisms" "Eukaryota" "Stramenopiles"      
    ## [7,] "sv22777" "NCBI" "cellular organisms" "Eukaryota" "Viridiplantae"      
    ##      [,6]                 [,7]                   [,8]                  
    ## [1,] "Apicomplexa"        "Conoidasida"          "Gregarinasina"       
    ## [2,] "Haptophyta"         "Prymnesiophyceae"     "Prymnesiales"        
    ## [3,] "Katablepharidaceae" ""                     ""                    
    ## [4,] "Choanoflagellata"   "Craspedida"           ""                    
    ## [5,] "Cercozoa"           "Chlorarachniophyceae" ""                    
    ## [6,] "Bacillariophyta"    "Coscinodiscophyceae"  "Chaetocerotophycidae"
    ## [7,] "Chlorophyta"        "Mamiellophyceae"      "Dolichomastigales"   
    ##      [,9]                   [,10]              [,11]         [,12]
    ## [1,] "Eugregarinorida"      "Porosporidae"     ""            ""   
    ## [2,] "Chrysochromulinaceae" "Chrysochromulina" ""            ""   
    ## [3,] ""                     ""                 ""            ""   
    ## [4,] ""                     ""                 ""            ""   
    ## [5,] ""                     ""                 ""            ""   
    ## [6,] "Chaetocerotales"      "Chaetocerotaceae" "Chaetoceros" ""   
    ## [7,] "Dolichomastigaceae"   ""                 ""            ""

``` r
new.silva <- cbind(lca.silva$V1, stringr::str_split(lca.silva$V2, pattern = ";", simplify = TRUE))
new.silva
```

    ##      [,1]      [,2]   [,3]                 [,4]        [,5]          
    ## [1,] "sv22340" "NCBI" "No hits"            ""          ""            
    ## [2,] "sv21398" "NCBI" "cellular organisms" "Eukaryota" "Haptista"    
    ## [3,] "sv22912" "NCBI" "cellular organisms" "Eukaryota" "Alveolata"   
    ## [4,] "sv23338" "NCBI" "cellular organisms" "Eukaryota" "Opisthokonta"
    ## [5,] "sv16360" "NCBI" "cellular organisms" "Eukaryota" "Opisthokonta"
    ## [6,] "sv12155" "NCBI" "No hits"            ""          ""            
    ## [7,] "sv22777" "NCBI" "cellular organisms" "Eukaryota" ""            
    ##      [,6]          [,7]                   [,8]           [,9]
    ## [1,] ""            ""                     ""             ""  
    ## [2,] "Haptophyta"  "Prymnesiophyceae"     "Prymnesiales" ""  
    ## [3,] "Dinophyceae" "Syndiniales"          ""             ""  
    ## [4,] "Fungi"       "Fungi incertae sedis" ""             ""  
    ## [5,] "Metazoa"     "Eumetazoa"            "Bilateria"    ""  
    ## [6,] ""            ""                     ""             ""  
    ## [7,] ""            ""                     ""             ""

Looks like the 2nd, 3rd, and last columns are not useful. If we had real data we would need to confirm that by inspecting the unique values, but since we have 7 ASVs we can see them all here.

Here we'll reformat to get these ready to plug into *taxmapper*:

``` r
lca.pr2 <- as.data.frame(new.pr2[ , -c(2,3, ncol(new.pr2))], stringsAsFactors = FALSE)
lca.silva <- as.data.frame(new.silva[ , -c(2,3, ncol(new.silva))], stringsAsFactors = FALSE)

colnames(lca.pr2) <- c("svN", "rank1", "rank2", "rank3", "rank4", "rank5", "rank6", 
                       "rank7", "rank8")
colnames(lca.silva) <- c("svN", "rank1", "rank2", "rank3", "rank4", "rank5")

for (i in colnames(lca.pr2)){
  empty.name <- which(str_length(lca.pr2[ , i]) == 0)
  lca.pr2[empty.name , i] <- NA
}
for (i in colnames(lca.silva)){
  empty.name <- which(str_length(lca.silva[ , i]) == 0)
  lca.silva[empty.name , i] <- NA
}

lca.pr2
```

    ##       svN     rank1               rank2              rank3                rank4
    ## 1 sv22340 Eukaryota           Alveolata        Apicomplexa          Conoidasida
    ## 2 sv21398 Eukaryota            Haptista         Haptophyta     Prymnesiophyceae
    ## 3 sv22912 Eukaryota Katablepharidophyta Katablepharidaceae                 <NA>
    ## 4 sv23338 Eukaryota        Opisthokonta   Choanoflagellata           Craspedida
    ## 5 sv16360 Eukaryota            Rhizaria           Cercozoa Chlorarachniophyceae
    ## 6 sv12155 Eukaryota       Stramenopiles    Bacillariophyta  Coscinodiscophyceae
    ## 7 sv22777 Eukaryota       Viridiplantae        Chlorophyta      Mamiellophyceae
    ##                  rank5                rank6            rank7       rank8
    ## 1        Gregarinasina      Eugregarinorida     Porosporidae        <NA>
    ## 2         Prymnesiales Chrysochromulinaceae Chrysochromulina        <NA>
    ## 3                 <NA>                 <NA>             <NA>        <NA>
    ## 4                 <NA>                 <NA>             <NA>        <NA>
    ## 5                 <NA>                 <NA>             <NA>        <NA>
    ## 6 Chaetocerotophycidae      Chaetocerotales Chaetocerotaceae Chaetoceros
    ## 7    Dolichomastigales   Dolichomastigaceae             <NA>        <NA>

``` r
lca.silva
```

    ##       svN     rank1        rank2       rank3                rank4        rank5
    ## 1 sv22340      <NA>         <NA>        <NA>                 <NA>         <NA>
    ## 2 sv21398 Eukaryota     Haptista  Haptophyta     Prymnesiophyceae Prymnesiales
    ## 3 sv22912 Eukaryota    Alveolata Dinophyceae          Syndiniales         <NA>
    ## 4 sv23338 Eukaryota Opisthokonta       Fungi Fungi incertae sedis         <NA>
    ## 5 sv16360 Eukaryota Opisthokonta     Metazoa            Eumetazoa    Bilateria
    ## 6 sv12155      <NA>         <NA>        <NA>                 <NA>         <NA>
    ## 7 sv22777 Eukaryota         <NA>        <NA>                 <NA>         <NA>

``` r
class(lca.pr2)
```

    ## [1] "data.frame"

``` r
typeof(lca.pr2)
```

    ## [1] "list"

``` r
class(lca.silva)
```

    ## [1] "data.frame"

``` r
typeof(lca.silva)
```

    ## [1] "list"

Looks good, except for the disagreements between the 2 tables. That's for another day though. This is how our data should look before plugging into *taxmapper*.

Running taxmapper
-----------------

MEGAN 6 actually maps all it's assignments onto the NCBI taxonomy, so we could go straight to ensembleTax by adding columns of NA's to lca.silva to give it the same number of ranks as lca.pr2.

But this is a demo so we'll map to the pr2 taxonomy first just for fun:

``` r
library("ensembleTax")
lca.pr2.m <- taxmapper(lca.pr2, tt.ranks = colnames(lca.pr2[2:ncol(lca.pr2)]),
                       tax2map2 = "pr2", 
                       exceptions = NULL, 
                       ignore.format = TRUE,
                       synonym.file = "default",
                       streamline = TRUE)
lca.silva.m <- taxmapper(lca.silva, tt.ranks = colnames(lca.pr2[2:ncol(lca.pr2)]),
                       tax2map2 = "pr2", 
                       exceptions = NULL, 
                       ignore.format = TRUE,
                       synonym.file = "default",
                       streamline = TRUE)

lca.pr2.m
```

    ##       svN   kingdom     supergroup            division                class
    ## 1 sv12155 Eukaryota  Stramenopiles          Ochrophyta      Bacillariophyta
    ## 2 sv16360 Eukaryota       Rhizaria            Cercozoa Chlorarachniophyceae
    ## 3 sv21398 Eukaryota       Hacrobia          Haptophyta     Prymnesiophyceae
    ## 4 sv22340 Eukaryota      Alveolata         Apicomplexa     Gregarinomorphea
    ## 5 sv22777 Eukaryota Archaeplastida         Chlorophyta      Mamiellophyceae
    ## 6 sv22912 Eukaryota       Hacrobia Katablepharidophyta   Katablepharidaceae
    ## 7 sv23338 Eukaryota   Opisthokonta    Choanoflagellida    Choanoflagellatea
    ##               order                     family            genus species
    ## 1 Bacillariophyta_X Polar-centric-Mediophyceae      Chaetoceros    <NA>
    ## 2              <NA>                       <NA>             <NA>    <NA>
    ## 3      Prymnesiales       Chrysochromulinaceae Chrysochromulina    <NA>
    ## 4   Eugregarinorida               Porosporidae             <NA>    <NA>
    ## 5 Dolichomastigales         Dolichomastigaceae             <NA>    <NA>
    ## 6              <NA>                       <NA>             <NA>    <NA>
    ## 7        Craspedida                       <NA>             <NA>    <NA>

``` r
lca.silva.m
```

    ##       svN   kingdom   supergroup       division            class        order
    ## 1 sv12155      <NA>         <NA>           <NA>             <NA>         <NA>
    ## 2 sv16360 Eukaryota Opisthokonta        Metazoa             <NA>         <NA>
    ## 3 sv21398 Eukaryota     Hacrobia     Haptophyta Prymnesiophyceae Prymnesiales
    ## 4 sv22340      <NA>         <NA>           <NA>             <NA>         <NA>
    ## 5 sv22777 Eukaryota         <NA>           <NA>             <NA>         <NA>
    ## 6 sv22912 Eukaryota    Alveolata Dinoflagellata      Syndiniales         <NA>
    ## 7 sv23338 Eukaryota Opisthokonta          Fungi             <NA>         <NA>
    ##   family genus species
    ## 1   <NA>  <NA>    <NA>
    ## 2   <NA>  <NA>    <NA>
    ## 3   <NA>  <NA>    <NA>
    ## 4   <NA>  <NA>    <NA>
    ## 5   <NA>  <NA>    <NA>
    ## 6   <NA>  <NA>    <NA>
    ## 7   <NA>  <NA>    <NA>

Looks good, most of the names were retained through the mapping and just popped onto the pr2 ranking conventions.

Calculating ensemble taxonomic assignments
------------------------------------------

Now let's get the ensemble assignments. We'll break ties by favoring lca.pr2.m (which means our ensemble will be identical to lca.pr2.m in this case):

``` r
xx <- list(lca.pr2.m, lca.silva.m)
names(xx) <- c("lca-pr2", "lca-silva")
eTax <- ensembleTax(xx, tiebreakz = c("lca-pr2"))
eTax
```

    ##       svN   kingdom     supergroup            division                class
    ## 1 sv12155 Eukaryota  Stramenopiles          Ochrophyta      Bacillariophyta
    ## 2 sv16360 Eukaryota       Rhizaria            Cercozoa Chlorarachniophyceae
    ## 3 sv21398 Eukaryota       Hacrobia          Haptophyta     Prymnesiophyceae
    ## 4 sv22340 Eukaryota      Alveolata         Apicomplexa     Gregarinomorphea
    ## 5 sv22777 Eukaryota Archaeplastida         Chlorophyta      Mamiellophyceae
    ## 6 sv22912 Eukaryota       Hacrobia Katablepharidophyta   Katablepharidaceae
    ## 7 sv23338 Eukaryota   Opisthokonta    Choanoflagellida    Choanoflagellatea
    ##               order                     family            genus species
    ## 1 Bacillariophyta_X Polar-centric-Mediophyceae      Chaetoceros    <NA>
    ## 2              <NA>                       <NA>             <NA>    <NA>
    ## 3      Prymnesiales       Chrysochromulinaceae Chrysochromulina    <NA>
    ## 4   Eugregarinorida               Porosporidae             <NA>    <NA>
    ## 5 Dolichomastigales         Dolichomastigaceae             <NA>    <NA>
    ## 6              <NA>                       <NA>             <NA>    <NA>
    ## 7        Craspedida                       <NA>             <NA>    <NA>

``` r
identical(eTax, lca.pr2.m)
```

    ## [1] TRUE

Checks out, so we're good! Alriiiight.
