library(MetaboAnnotation)
library(msdata)
library(Spectra)
library(testthat)

fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
pest_ms2 <- filterMsLevel(Spectra(fl), 2L)
pest_ms2 <- pest_ms2[c(808, 809, 945:955)]
load(system.file("extdata", "minimb.RData", package = "MetaboAnnotation"))
