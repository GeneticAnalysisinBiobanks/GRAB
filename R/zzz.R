## This file declares global variables to avoid R CMD check warnings
## for variables used in non-standard evaluation (NSE) in data.table, dplyr, etc.
## See ?globalVariables for details. This is required for CRAN compliance.
utils::globalVariables(c(
  "ID", "AF_ref", "ALT", "ALT.x", "ALT.y", "Annos", "Cov", "E_pos1", "E_pos2", "G_Cov", "G_GxE_Cov1",
  "G_GxE_Cov2", "G_pos1", "G_pos2", "GxE_Cov", "GxE_pos1", "GxE_pos2", "ID.x", "ID.y", "ID1", "ID2",
  "IndicatorColumn", "IndicatorVec", "Info", "MAF", "Outlier", "RA", "REF", "REF.x", "REF.y", "Resid",
  "SampleIDColumn", "StatVec", "SurvTimeColumn", "V", "V1", "V2", "Value", "VarCorr", "add.edges",
  "adjVarSVec", "all_of", "altFreq", "as_tibble", "betaWeights", "delete.edges", "get.data.frame",
  "get.edgelist", "get.vertex.attribute", "graph_from_data_frame", "index",
  "make_graph", "missingRate", "mr0", "mr1", "mst", "mu.int", "mu0", "mu1", "n", "obj.mainRegion",
  "posRow", "pval1Vec", "r0", "summarize", "sym", "."
))