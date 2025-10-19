#!/usr/local/bin/Rscript
library(MetaCycle)
meta2d(infile = "metadata.csv", outdir = "metaout", filestyle = "csv", timepoints = "Line1", minper = 1080,
       maxper = 1800, cycMethod = c("ARS", "JTK", "LS"),
       analysisStrategy = "auto", outputFile = TRUE,
       outIntegration = "both", adjustPhase = "predictedPer",
       combinePvalue = "fisher", weightedPerPha = FALSE, ARSmle = "auto",
       ARSdefaultPer = 1440, outRawData = FALSE, releaseNote = TRUE,
       outSymbol = "", parallelize = FALSE, nCores = 1, inDF = NULL)