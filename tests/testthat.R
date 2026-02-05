library(MetaboAnnotation)
library(MsDataHub)
library(Spectra)
library(testthat)

fl <- MsDataHub::PestMix1_DDA.mzML()
pest_ms2 <- filterMsLevel(Spectra(fl), 2L)
pest_ms2 <- pest_ms2[c(808, 809, 945:955)]
load(system.file("extdata", "minimb.RData", package = "MetaboAnnotation"))

test_check("MetaboAnnotation")
