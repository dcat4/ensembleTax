## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library("ensembleTax")
library("Biostrings")

data("idtax.pr2.sample")
data("idtax.silva.sample")
data("bayes.sample")
data("rubric.sample")

head(idtax.pr2.sample)
head(idtax.silva.sample)
head(bayes.sample)
head(rubric.sample)

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

head(idtax.pr2.pretty)
head(idtax.silva.pretty)
head(bayes.pr2.pretty)


## -----------------------------------------------------------------------------

idtax.silva.mapped2pr2 <- taxmapper(idtax.silva.pretty,
                      tt.ranks = colnames(idtax.silva.pretty)[3:ncol(idtax.silva.pretty)],
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
eTax1 <- ensembleTax(xx, 
                     tablenames = names(xx), 
                     ranknames = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                     tiebreakz = "none", 
                     count.na=TRUE, 
                     assign.threshold = 0, 
                     weights=rep(1,length(xx)))
head(eTax1)


## -----------------------------------------------------------------------------
eTax2 <- ensembleTax(xx, 
                     tablenames = names(xx), 
                     ranknames = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                     tiebreakz = "none", 
                     count.na=FALSE, 
                     assign.threshold = 0, 
                     weights=c(2,1,1))
head(eTax2)

## -----------------------------------------------------------------------------
eTax3 <- ensembleTax(xx, 
                     tablenames = names(xx), 
                     ranknames = c("kingdom", "supergroup", "division","class","order","family","genus","species"),
                     tiebreakz = c("bayes-pr2"), 
                     count.na=TRUE, 
                     assign.threshold = 0, 
                     weights=c(1,1,2))
head(eTax3)

