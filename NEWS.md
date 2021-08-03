NEWS
================
D Catlett
08/03/2021

# ensembleTax v1.2.1

1.  pr2 taxonomy updated to v4.14.0
2.  UNITE taxonomy now available in package data for use with taxmapper
3.  RDP taxonomy updated to v18

# ensembleTax v1.1.1

1.  ensembleTax algorithm name changed to assign.ensembleTax to avoid confusion. Behavior is identical but calling the ensembleTax now results in an error
2.  ignore.format argument of taxmapper has undergone big changes. See documentation and taxmapper vignette for details.
3.  Bug fix for handling ambiguous placeholder names (eg, "Clade\_X", "incertae sedis", etc) in taxmapper. See documentation for details
4.  new taxmapper & assign.ensembleTax vignettes to better illustrate the functionality and behavior of both, included both on Github and in the R package
5.  incorporated initial taxonomic assignment and comparisons of individual method predictions with ensembleTax in package overview vignettes available on Github and in R package

# ensembleTax v1.0.2

Documentation updated for 3rd CRAN submission

# ensembleTax v1.0.1

License issue fixed for CRAN submission

# ensembleTax v1.0.0

Minor documentation cleaning and code-commenting to prepare for CRAN submission

# ensembleTax v0.0.51

Minor changes only - documentation cleaning and updating.

# ensembleTax v0.0.5

## Major changes

1.  Added and documented RDP and GreenGenes databases
2.  Major re-write of ensembleTax to prevent 'frankenstein' assignments when count.na = TRUE.

## Minor changes

1.  Removed redundant operations in ensembleTax
2.  Updated ensembleTax doc
3.  taxmapper has updated doc, allows NULL synonym.file, works with RDP and GG
4.  Cleaned up vignette
5.  Updated 2df functions to work with RDP and GG

## Other

1.  Corrected License
2.  added NEWS.md

# ensembleTax v0.0.4 and earlier

Dylan was a bad person then and didn't document them well.
