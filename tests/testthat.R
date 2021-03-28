library(msdata)
library(Spectra)
library(testthat)

fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
pest_ms2 <- filterMsLevel(Spectra(fl), 2L)

library(MsBackendMassbank)
library(RMariaDB)
con <- dbConnect(MariaDB(), host = "host.docker.internal",
                 dbname = "MassBank", user = "massbank",
                 password = "massbank")

be <- backendInitialize(MsBackendMassbankSql(), dbcon = con)
massbank <- filterPolarity(Spectra(be), polarity = 1)
massbank <- setBackend(massbank, MsBackendDataFrame())
