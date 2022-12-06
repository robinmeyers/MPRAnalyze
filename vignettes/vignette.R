## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load and examine data, include=TRUE--------------------------------------
library(MPRAnalyze)
data("ChrEpi")
summary(ce.colAnnot)
head(ce.colAnnot)

## ----init object, include=TRUE------------------------------------------------
obj <- MpraObject(dnaCounts = ce.dnaCounts, rnaCounts = ce.rnaCounts, 
                  dnaAnnot = ce.colAnnot, rnaAnnot = ce.colAnnot, 
                  controls = ce.control)

## ----library size estimation--------------------------------------------------
## If the library factors are different for the DNA and RNA data, separate 
## estimation of these factors is needed. We can also change the estimation 
## method (Upper quartile by default)
obj <- estimateDepthFactors(obj, lib.factor = c("batch", "condition"),
                            which.lib = "dna", 
                            depth.estimator = "uq")
obj <- estimateDepthFactors(obj, lib.factor = c("condition"),
                            which.lib = "rna", 
                            depth.estimator = "uq")

## In this case, the factors are the same - each combination of batch and 
## condition is a single library, and we'll use the default estimation
obj <- estimateDepthFactors(obj, lib.factor = c("batch", "condition"),
                            which.lib = "both")

