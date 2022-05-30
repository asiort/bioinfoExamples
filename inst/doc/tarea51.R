## -----------------------------------------------------------------------------
doAll = FALSE

## ---- eval=doAll--------------------------------------------------------------
#  library(pacman)
#  library(asierortega)
#  pacman::p_load(EnrichmentBrowser)
#  pacman::p_load(Biobase)
#  pacman::p_load(SummarizedExperiment)

## ---- eval = doAll------------------------------------------------------------
#  data(gse115484, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  gse115484 = asierortega::unifyeset(gse115484, ID = "ENTREZID")

## ---- eval = doAll------------------------------------------------------------
#  se115484 = makeSummarizedExperimentFromExpressionSet(gse115484)
#  se115484 = probe2gene(se115484)

## ---- eval = doAll------------------------------------------------------------
#  se115484$GROUP = se115484$Tissue
#  levels(se115484$GROUP)[levels(se115484$GROUP)=="Gluteus medius"] = 0
#  levels(se115484$GROUP)[levels(se115484$GROUP)=="Liver"] = 1

## ---- eval = doAll------------------------------------------------------------
#  se115484 = deAna(expr = se115484)
#  head(rowData(se115484))

## ---- eval = doAll------------------------------------------------------------
#  sscKEGGgsc= getGenesets("ssc", db = "kegg")
#  head(names(sscKEGGgsc))
#  save(sscKEGGgsc, file = "sscKEGGgsc.rda")

## ---- eval = doAll------------------------------------------------------------
#  data(sscKEGGgsc, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  se115484.oraKEGG= sbea(method = "ora", se = se115484, gs = sscKEGGgsc,
#                        perm = 0, alpha = 0.05)
#  gsRanking(se115484.oraKEGG)

## ---- eval = doAll------------------------------------------------------------
#  pvalues = gsRanking(se115484.oraKEGG)$PVAL
#  pBH = p.adjust(pvalues, method = "BH")
#  pBonferroni = p.adjust(pvalues, method = "bonferroni")

## ---- eval = doAll------------------------------------------------------------
#  length(pvalues)
#  alpha = 0.05
#  significativos_p = which(pvalues < alpha)
#  length(significativos_p)
#  
#  significativos_BH = which(pBH < alpha)
#  length(significativos_BH)
#  
#  significativos_Bonferroni =  which(pBonferroni < alpha)
#  length(significativos_Bonferroni)

## ---- eval = doAll------------------------------------------------------------
#  cabeceras = c(colnames(gsRanking(se115484.oraKEGG)), "p.adjust_BH", "p.adjust_Bonferroni")
#  kegg_url = asierortega::kegg2url(gsRanking(se115484.oraKEGG)[, "GENE.SET"])
#  
#  ## Construimos dataframe
#  df_kegg = data.frame(gsRanking(se115484.oraKEGG)[, "GENE.SET"], gsRanking(se115484.oraKEGG)[, "NR.GENES"], gsRanking(se115484.oraKEGG)[, "NR.SIG.GENES"],  gsRanking(se115484.oraKEGG)[, "PVAL"], pBH, pBonferroni)
#  
#  df_kegg_hip = data.frame(kegg_url, gsRanking(se115484.oraKEGG)[, "NR.GENES"], gsRanking(se115484.oraKEGG)[, "NR.SIG.GENES"],  gsRanking(se115484.oraKEGG)[, "PVAL"], pBH, pBonferroni)
#  
#  colnames(df_kegg) = cabeceras
#  head(df_kegg)

## ---- eval = doAll, warning = FALSE-------------------------------------------
#  df_kegg_prueba = overex2df(se115484, sscKEGGgsc, hiperlink = FALSE)
#  df_hip_kegg_prueba = overex2df(se115484, sscKEGGgsc, hiperlink = TRUE)

## ---- eval=doAll--------------------------------------------------------------
#  makeReport(df_kegg_hip, name = "Informe_GroupAnalysis_GSE115484", path = "../docs/articles")

## ---- eval = doAll------------------------------------------------------------
#  sscGOBPgsc = getGenesets(org = "ssc",db = "go", onto="BP")
#  save(sscGOBPgsc, file = "sscGOBPgsc.rda")

## ---- eval = doAll------------------------------------------------------------
#  data(sscGOBPgsc, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  se115484.GOBPora = sbea(method = "ora", se = se115484, gs = sscGOBPgsc,
#                         perm = 0, alpha = 0.05)
#  gsRanking(se115484.GOBPora)

## ---- eval = doAll------------------------------------------------------------
#  data(PRJNA548871, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  PRJNA548871 = cleanNA(PRJNA548871, Ncolumns = 3)

## ---- eval = doAll------------------------------------------------------------
#  PRJNA548871 = probe2gene(PRJNA548871, from = "ENSEMBL", to = "ENTREZID")

## ---- eval = doAll------------------------------------------------------------
#  PRJNA548871$GROUP = PRJNA548871$treatment
#  levels(PRJNA548871$GROUP )[levels(PRJNA548871$GROUP )=="nonOIR"] = 0
#  levels(PRJNA548871$GROUP )[levels(PRJNA548871$GROUP )=="OIR"] = 1

## ---- eval = doAll------------------------------------------------------------
#  PRJNA548871 = deAna(expr = PRJNA548871)
#  head(rowData(PRJNA548871))

## ---- eval = doAll------------------------------------------------------------
#  mmuKEGGgsc= getGenesets("mmu",db = "kegg")
#  head(names(mmuKEGGgsc))
#  save(mmuKEGGgsc, file = "mmuKEGGgsc.rda")

## ---- eval = doAll------------------------------------------------------------
#  data(mmuKEGGgsc, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  PRJNA548871.oraKEGG= sbea(method = "ora", se = PRJNA548871, gs = mmuKEGGgsc,
#                        perm = 0, alpha = 0.05)
#  gsRanking(PRJNA548871.oraKEGG)

## ---- eval = doAll------------------------------------------------------------
#  pvalues = gsRanking(PRJNA548871.oraKEGG)$PVAL
#  pBH = p.adjust(pvalues, method = "BH")
#  pBonferroni = p.adjust(pvalues, method = "bonferroni")

## ---- eval = doAll------------------------------------------------------------
#  cabeceras_rna = c(colnames(gsRanking(PRJNA548871.oraKEGG)), "p.adjust_BH", "p.adjust_Bonferroni")
#  kegg_url_rna = asierortega::kegg2url(gsRanking(PRJNA548871.oraKEGG)[, "GENE.SET"])
#  
#  ## Construimos dataframe
#  df_kegg_rna = data.frame(gsRanking(PRJNA548871.oraKEGG)[, "GENE.SET"],  gsRanking(PRJNA548871.oraKEGG)[, "NR.GENES"],  gsRanking(PRJNA548871.oraKEGG)[, "NR.SIG.GENES"],  gsRanking(PRJNA548871.oraKEGG)[, "PVAL"], pBH, pBonferroni)
#  
#  df_kegg_hip_rna = data.frame(kegg_url_rna, gsRanking(PRJNA548871.oraKEGG)[, "NR.GENES"], gsRanking(PRJNA548871.oraKEGG)[, "NR.SIG.GENES"],  gsRanking(PRJNA548871.oraKEGG)[, "PVAL"], pBH, pBonferroni)
#  
#  colnames(df_kegg_rna) = cabeceras_rna
#  colnames(df_kegg_hip_rna) = cabeceras_rna
#  head(df_kegg_rna)

## ---- eval = doAll, warning = FALSE-------------------------------------------
#  data(PRJNA548871, package = "asierortega")
#  PRJNA548871 = cleanNA(PRJNA548871, Ncolumns = 3)
#  PRJNA548871 = probe2gene(PRJNA548871, from = "ENSEMBL", to = "ENTREZID")
#  PRJNA548871$GROUP = PRJNA548871$treatment
#  levels(PRJNA548871$GROUP )[levels(PRJNA548871$GROUP )=="nonOIR"] = 0
#  levels(PRJNA548871$GROUP )[levels(PRJNA548871$GROUP )=="OIR"] = 1
#  
#  df_kegg_rna_prueba = overex2df(PRJNA548871, mmuKEGGgsc, hiperlink = FALSE)
#  df_hip_kegg_rna_prueba = overex2df(PRJNA548871, mmuKEGGgsc, hiperlink = TRUE)

## ---- eval=doAll--------------------------------------------------------------
#  makeReport(df_kegg_hip_rna, name = "Informe_GroupAnalysis_PRJNA548871", path = "../docs/articles")

