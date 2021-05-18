taxmapper\_tutorial
================
D Catlett
5/15/2021

# Introduction

This vignette provides detailed examples to demonstrate the functionality of the *taxmapper* algorithm included with the ensembleTax package. For a more general demonstration of the ensembleTax package functionality/workflow, please go here: <https://github.com/dcat4/ensembleTax>

## The taxmapper algorithm

*taxmapper*'s purpose is to map a collection of taxonomic assignments onto a different taxonomic nomenclature (set of naming and ranking conventions).

It does this via rank-agnostic exact name matching. In other words, *taxmapper* doesn't care about the heirarchical structure of a taxonomic nomenclature, and assumes that a taxonomic name means the same thing regardless of which reference database that name is found in.

### Examples

To demonstrate the functionality of *taxmapper*, we'll create an artificial set of taxonomic assignments to mimic those you might obtain with, for example, dada2's *assignTaxonomy* function, as well as an artificial taxonomic nomenclature that mimic's those available in the ensembleTax R package.

So first, load the ensembleTax package, and create the artificial data sets:

``` r
library("ensembleTax")
packageVersion("ensembleTax")
```

    ## [1] '1.1.1'

``` r
# create a fake taxonomy table of ASVs and taxonomic assignments
fake.taxtab <- data.frame(ASV = c("CGTC", "AAAA"),
                          kingdom = c("Eukaryota", "Bacteria"), 
                          supergroup = c("Stramenopile", NA),
                          division = c("Ochrophyta", NA),
                          class = c("Diatomea", NA),
                          genus = c("Pseudo-nitzschia", NA))
# create a fake taxonomic nomenclature:
map2me <- data.frame(kingdom = c("Eukaryota"), 
                     largegroup = c("Stramenopile"),
                     division = c("Clade_X"),
                     class = c("Ochrophyta"),
                     order = c("Bacillariophyta"),
                     genus = c("Pseudonitzschia"))
# look at your artificial data:
fake.taxtab
```

    ##    ASV   kingdom   supergroup   division    class            genus
    ## 1 CGTC Eukaryota Stramenopile Ochrophyta Diatomea Pseudo-nitzschia
    ## 2 AAAA  Bacteria         <NA>       <NA>     <NA>             <NA>

``` r
map2me
```

    ##     kingdom   largegroup division      class           order           genus
    ## 1 Eukaryota Stramenopile  Clade_X Ochrophyta Bacillariophyta Pseudonitzschia

So we see we have a set of 2 ASVs with taxonomic assignments, and a taxonomic nomenclature that contains 1 taxonomic entry (we're trying to make a simple example here; most reference databases have thousands of taxonomic entries).

Now would be a good time to review the *taxmapper* documentation to get a sense of the different parameter spaces available. Here we'll try to demonstrate what these different parameters are doing.

#### Example 1: Strict exact name-matching and the "streamline" argument

To start, we'll run *taxmapper* with no exceptions, no format-ignoring, and no taxonomic synonyms, and we'll look at the different outputs you can expect based on the *streamline* argument:

``` r
mapped.tt.stmlin <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = NULL,
                       ignore.format = FALSE,
                       synonym.file = NULL,
                       streamline = TRUE)
mapped.tt.stmlin
```

    ##    ASV   kingdom   largegroup division      class order genus
    ## 1 AAAA      <NA>         <NA>     <NA>       <NA>  <NA>  <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta  <NA>  <NA>

``` r
mapped.tt.no.stmlin <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = NULL,
                       ignore.format = FALSE,
                       synonym.file = NULL,
                       streamline = FALSE)
mapped.tt.no.stmlin
```

    ## [[1]]
    ##     kingdom   supergroup   division    class            genus tax2map2_kingdom
    ## 1 Eukaryota Stramenopile Ochrophyta Diatomea Pseudo-nitzschia        Eukaryota
    ##   tax2map2_largegroup tax2map2_division tax2map2_class tax2map2_order
    ## 1        Stramenopile           Clade_X     Ochrophyta             NA
    ##   tax2map2_genus
    ## 1             NA
    ## 
    ## [[2]]
    ## [1] "Pseudo-nitzschia" "Diatomea"         "Bacteria"        
    ## 
    ## [[3]]
    ##    ASV   kingdom   largegroup division      class order genus
    ## 1 AAAA      <NA>         <NA>     <NA>       <NA>  <NA>  <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta  <NA>  <NA>

We see that when *streamline = TRUE* we return a dataframe with the input ASV's and their mapped taxonomic assignments. This is intended for users who want to automate their ensembleTax workflow and move on with further analyses right away.

If you want to take a look "under the hood", setting *streamline = FALSE*, returns a 3-element list. \[\[1\]\] shows the input taxonomic assignments aligned with their mapped values (a sort of mapping "rubric"). Because *Bacteria* was not found in tax2map2, it does not have a taxonomy to map onto and is not included in the "rubric". \[\[2\]\] shows the taxonomic names that could not be mapped. We see that these are the names that were not found (or did not have exact matches to any name) in *tax2map2* (or "map2me" in this example). Finally, \[\[3\]\] contains the mapped input taxonomy table, which is identical to what was returned when *streamline = TRUE*.

We see that the "CGTC" ASV was mapped to *Ochrophyta*, despite the use of different ranking conventions in the input taxonomy table and the taxonomic nomenclature we're mapping onto. This illustrates the "rank-agnostic" part of *taxmapper*. The "AAAA" ASV is entirely unassigned in the mapped output because our *tax2map2* didn't include *Bacteria*. If you'd like to retain a high-level taxonomic assignment like *Bacteria* in this example, you can address that with the *exceptions* argument.

#### Example 2: The "exceptions" argument

Here we'll specify that we want to keep *Bacteria* assignments even though they aren't included in *tax2map2*:

``` r
mapped.tt.exc <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = c("Bacteria"),
                       ignore.format = FALSE,
                       synonym.file = NULL,
                       streamline = TRUE)
mapped.tt.exc
```

    ##    ASV   kingdom   largegroup division      class order genus
    ## 1 AAAA  Bacteria         <NA>     <NA>       <NA>  <NA>  <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta  <NA>  <NA>

And we see that instead of a completely unassigned "AAAA" ASV as we had above, we've now retained the *Bacteria* assignment in the mapped output.

#### Example 3: Incorporating taxonomic synonyms

Folks who study phytoplankton might recognize that *Diatomea* in our input taxonomy table and *Bacillariophyta* in the nomenclature we're mapping onto are taxonomic synonyms (both refer to the same class of phytoplankton, diatoms). *taxmapper* can search for taxonomic synonyms as well.

If you'd like to use a custom compilation of taxonomic synonyms, please see this vignette: <https://github.com/dcat4/ensembleTax/blob/master/how_to_add_synonyms.md>.

ensembleTax includes a collection of pre-compiled eukaryotic taxonomic synonyms. Let's have a look at whether *Diatomea* and *Bacillariophyta* are included in this pre-compiled data set:

``` r
# load ensembleTax's pre-compiled synonyms:
syn.df <- ensembleTax::synonyms_20200816
# pull rows with Diatomea (there's only 1)
diatom.synonyms <- syn.df[which(syn.df == "Diatomea", arr.ind=TRUE)[,'row'],]
# look at it:
diatom.synonyms
```

    ##      Name_1          Name_2 Name_3 Name_4 Name_5 Name_6 Name_7      References
    ## 10 Diatomea Bacillariophyta   <NA>   <NA>   <NA>   <NA>   <NA> Adl et al. 2012
    ##    Notes.References X X.1
    ## 10

They are. You can follow a similar procedure to check for synonyms for your favorite taxonomic name, or enhance our pre-compiled synonym collection by saving the above syn.df dataframe to a csv and adding in your own collections of synonyms.

Moving on, if we tell *taxmapper* to consult the pre-compiled taxonomic synonyms included with the ensembleTax package, we should be able to get more refined mapped taxonomic assignments in this example. We'll do this here with the *synonym.file = "default"* argument:

``` r
mapped.tt.syn <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = c("Bacteria"),
                       ignore.format = FALSE,
                       synonym.file = "default",
                       streamline = TRUE)
mapped.tt.syn
```

    ##    ASV   kingdom   largegroup division      class           order genus
    ## 1 AAAA  Bacteria         <NA>     <NA>       <NA>            <NA>  <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta Bacillariophyta  <NA>

Taking a look at this output, we see that the "CGTC" ASV has now been mapped to *Bacillariophyta*, despite the fact that it is called *Diatomea* in the fake reference database we used to generate our fake taxonomic assignments. So our inclusion of taxonomic synonyms has reduced the information lost in taxonomy mapping.

We have just one more parameter to check out...

#### Example 4: The "ignore.format" argument

You might have noticed in the examples above that our input taxonomy table includes an ASV assigned as *Pseudo-nitzschia*, while the nomenclature we're mapping to includes the same taxonomic name with no hyphen in the middle. This is where the ignore.format argument can be helpful:

``` r
mapped.tt.igfo <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = c("Bacteria"),
                       ignore.format = TRUE,
                       synonym.file = NULL,
                       streamline = TRUE)
mapped.tt.igfo
```

    ##    ASV   kingdom   largegroup division      class           order
    ## 1 AAAA  Bacteria         <NA>     <NA>       <NA>            <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta Bacillariophyta
    ##             genus
    ## 1            <NA>
    ## 2 Pseudonitzschia

We see that setting *ignore.format = TRUE* has circumvented the formatting issue, and now we retain more information in our mapped annotations since we're able to map *Pseudo-nitzschia* onto *Pseudonitzschia*. Other special symbols handled with *ignore.format = TRUE* include " " (single space), "\_", "-", "\[", "\]". It also reduces case sensitivity (attempts to map all-lower- and all-upper- case variants of a taxonomic name).

If you read the *ignore.format* documentation carefully, you may notice there are other circumstances where the *ignore.format* option doesn't work as cleanly. Here we'll show an example to illustrate.

If the special characters *ignore.format* handles are found in *tax2map2* rather than *tt*, the mapping won't work. We'll make a second fake.taxtab and map2me with the *Pseudonitzschia* variants swapped to demonstrate:

``` r
fake.taxtab2 <- fake.taxtab
fake.taxtab2[fake.taxtab2 == "Pseudo-nitzschia"] <- "Pseudonitzschia"
map2me2 <- map2me
map2me2[map2me2 == "Pseudonitzschia"] <- "Pseudo-nitzschia"

fake.taxtab2
```

    ##    ASV   kingdom   supergroup   division    class           genus
    ## 1 CGTC Eukaryota Stramenopile Ochrophyta Diatomea Pseudonitzschia
    ## 2 AAAA  Bacteria         <NA>       <NA>     <NA>            <NA>

``` r
map2me2
```

    ##     kingdom   largegroup division      class           order            genus
    ## 1 Eukaryota Stramenopile  Clade_X Ochrophyta Bacillariophyta Pseudo-nitzschia

``` r
mapped.tt.igfo2 <- taxmapper(tt = fake.taxtab2,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me2,
                       exceptions = c("Bacteria"),
                       ignore.format = TRUE,
                       synonym.file = NULL,
                       streamline = TRUE)
mapped.tt.igfo2
```

    ##    ASV   kingdom   largegroup division      class order genus
    ## 1 AAAA  Bacteria         <NA>     <NA>       <NA>  <NA>  <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta  <NA>  <NA>

This example illustrates that formatting is only being ignored for the taxonomic names we're mapping, and NOT for the taxonomic nomenclature we're mapping onto. This is an important limitation to keep in mind. If you find this problematic, you may consider further customization of the *tax2map2* data. We are considering more detailed manipulations of the nomenclatures supported by ensembleTax to circumvent this issue but for now we supply these exactly as they are supplied by the creators of the reference databases.

Here's a few more examples to demonstrate some capabilities of

#### Example 5: ambiguous "placeholder" names

One last example we need to look at considers ambiguous taxonomic names that are sometimes included in reference databases.

Let's make a small adjustment to our fake.taxtab to see how these are handled by *taxmapper*. We'll add a "Clade\_X" supergroup annotation to our prokaryotic ASV.

``` r
# create a new fake taxonomy table of ASVs and taxonomic assignments
fake.taxtab <- data.frame(ASV = c("CGTC", "AAAA"),
                          kingdom = c("Eukaryota", "Bacteria"), 
                          supergroup = c("Stramenopile", "Clade_X"),
                          division = c("Ochrophyta", NA),
                          class = c("Diatomea", NA),
                          genus = c("Pseudo-nitzschia", NA))

# look at your artificial data again:
fake.taxtab
```

    ##    ASV   kingdom   supergroup   division    class            genus
    ## 1 CGTC Eukaryota Stramenopile Ochrophyta Diatomea Pseudo-nitzschia
    ## 2 AAAA  Bacteria      Clade_X       <NA>     <NA>             <NA>

``` r
map2me
```

    ##     kingdom   largegroup division      class           order           genus
    ## 1 Eukaryota Stramenopile  Clade_X Ochrophyta Bacillariophyta Pseudonitzschia

Re-inspecting map2me shows that "Clade\_X" is also the name of a clade of Eukaryotic Stramenopiles. Ruh roh. This might introduce errors in the mapped taxonomic assignments since Clade\_X is a name found in both a Bacterial and Stramenopile lineage.

Let's see what happens when we run taxmapper:

``` r
mapped.tt.ambigtest <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = NULL,
                       ignore.format = TRUE,
                       synonym.file = NULL,
                       streamline = TRUE)
mapped.tt.ambigtest
```

    ##    ASV   kingdom   largegroup division      class           order
    ## 1 AAAA      <NA>         <NA>     <NA>       <NA>            <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta Bacillariophyta
    ##             genus
    ## 1            <NA>
    ## 2 Pseudonitzschia

We see that despite the fact that there was an exact name match, *taxmapper* has avoided making an incorrect annotation in the mapped output. *taxmapper* does this by checking the names to be mapped for taxonomic names that BEGIN with certain words. Here's the complete list of what it checks for: "Clade", "CLADE", "clade", "Group", "GROUP", "group", "Class", "CLASS", "class", "Subgroup", "SubGroup", "SUBGROUP", "subgroup", "Subclade", "SubClade", "SUBCLADE", "subclade", "Subclass", "SubClass", "SUBCLASS", "subclass", "incertae sedis", "INCERTAE SEDIS", "Incertae sedis", "Incertae Sedis", "incertae-sedis", "INCERTAE-SEDIS", "Incertae-sedis", "Incertae-Sedis", "incertae\_sedis", "INCERTAE\_-SEDIS", "Incertae\_sedis", "Incertae\_Sedis", "incertaesedis", "INCERTAESEDIS", "Incertaesedis", "IncertaeSedis", "unclassified", "UNCLASSIFIED", "Unclassified", "Novel", "novel", "NOVEL", "sp", "sp.", "spp", "spp.", "lineage", "Lineage", "LINEAGE"

So, what does *taxmapper* do when it encounters an ambiguous name like "Clade\_X"? It doesn't just discard the name. Instead, it finds the lowest rank with a non-ambiguous taxonomic name (a name that doesn't begin with a word in the list above), and appends that non-ambiguous name to the ambiguous name, separated by a "-". In our example above, this means *taxmapper* was searching for "Bacteria-Clade\_X" rather than just "Clade\_X", removing the ambiguity in taxonomic identity.

Here we'll add an annotation to our *tax2map2* (the map2me variable defined above) and see that, in some cases, we can use *ignore.format* to map the ambiguous "Clade\_X" name assigned to our "AAAA" ASV:

``` r
# add an entry in our tax2map2 that matches (but not exactly) one of our ASVs:
map2me <- rbind(map2me,
                c("Bacteria", "Bacteria_Clade_X", rep(NA, times = ncol(map2me)-2)))
map2me
```

    ##     kingdom       largegroup division      class           order
    ## 1 Eukaryota     Stramenopile  Clade_X Ochrophyta Bacillariophyta
    ## 2  Bacteria Bacteria_Clade_X     <NA>       <NA>            <NA>
    ##             genus
    ## 1 Pseudonitzschia
    ## 2            <NA>

``` r
# map again with ignore.format = FALSE.. the Bacteria will only map to Bacteria
mapped.tt.ambigtest2 <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = NULL,
                       ignore.format = FALSE,
                       synonym.file = NULL,
                       streamline = TRUE)
# confirm:
mapped.tt.ambigtest2
```

    ##    ASV   kingdom   largegroup division      class order genus
    ## 1 AAAA  Bacteria         <NA>     <NA>       <NA>  <NA>  <NA>
    ## 2 CGTC Eukaryota Stramenopile  Clade_X Ochrophyta  <NA>  <NA>

``` r
# now set ignore.format = TRUE.. we'll map to Bacteria Clade X:
mapped.tt.ambigtest3 <- taxmapper(tt = fake.taxtab,
                       tt.ranks = colnames(fake.taxtab)[2:ncol(fake.taxtab)],
                       tax2map2 = map2me,
                       exceptions = NULL,
                       ignore.format = TRUE,
                       synonym.file = NULL,
                       streamline = TRUE)
# confirm:
mapped.tt.ambigtest3
```

    ##    ASV   kingdom       largegroup division      class           order
    ## 1 AAAA  Bacteria Bacteria_Clade_X     <NA>       <NA>            <NA>
    ## 2 CGTC Eukaryota     Stramenopile  Clade_X Ochrophyta Bacillariophyta
    ##             genus
    ## 1            <NA>
    ## 2 Pseudonitzschia

To clarify what's going on here one last time, when *taxmapper* encountered the "Clade\_X" assignment for the "AAAA" ASV, it appended the next-lowest non-ambiguous taxonomic assignment ("Bacteria") and searched for an exact match to this now-non-ambiguous name ("Bacteria-Clade\_X"). When *ignore.format = FALSE*, "Bacteria-Clade\_X" was not an exact match to "Bacteria\_Clade\_X" (the hyphen and underscore are different). But when *ignore.format = TRUE*, *taxmapper* searched for various formatting variants of "Bacteria-Clade\_X", one of which is "Bacteria\_Clade\_X". This results in an exact match in *tax2map2* and a more refined mapped taxonomic annotation for this ASV.

You might notice that if an ambiguous name like "Clade\_X" is found in *tax2map2*, we will NOT be able to map onto this taxonomic assignment under any circumstances with the current implementation of *taxmapper*. The strategy *taxmapper* uses here is based on inspection of the database nomenclatures included in ensembleTax and our desire to preserve the nomenclatures employed by different reference databases as closely as possible. Again, we are considering more detailed manipulations of the nomenclatures supported by ensembleTax to circumvent this issue but for now we supply these as they are supplied by the creators of the reference databases.

And that brings us to the end of this vignette. Please let us know about issues that come up on the esembleTax Github issues tracker.
