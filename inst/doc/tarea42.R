## -----------------------------------------------------------------------------
doAll = FALSE

## ---- eval = doAll------------------------------------------------------------
#  library(pacman)
#  library(asierortega)
#  pacman::p_load(SummarizedExperiment)
#  pacman::p_load(DESeq2)
#  pacman::p_load(ggplot2)
#  pacman::p_load(ReportingTools)
#  pacman::p_load(pheatmap)

## ---- eval = doAll------------------------------------------------------------
#  data(PRJNA548871, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  PRJNA548871 = cleanNA(PRJNA548871, Ncolumns = 3)

## ---- eval = doAll------------------------------------------------------------
#  dds =  DESeqDataSet(se = PRJNA548871, design = ~ treatment)
#  cantidad_genes_pre = dim(counts(dds))[1]
#  keep = rowSums(counts(dds)) >= 10
#  dds = dds[keep,]
#  PRJNA548871_keep = PRJNA548871[keep,]
#  cantidad_genes_post = dim(counts(dds))[1]
#  (cantidad_genes_pre - cantidad_genes_post)

## ---- eval = doAll------------------------------------------------------------
#  dds$treatment <- relevel(dds$treatment, ref  = "nonOIR")
#  levels(dds$treatment)

## ---- eval = doAll, warning = FALSE-------------------------------------------
#  dds = DESeq(dds)

## ---- eval = doAll, warning = FALSE-------------------------------------------
#  res = lfcShrink(dds, coef=2)
#  pvalues = res$pvalue
#  pBH = p.adjust(pvalues, method = "BH")
#  pBonferroni = p.adjust(pvalues, method = "bonferroni")

## ---- eval = doAll------------------------------------------------------------
#  length(pvalues)
#  alpha = 0.05
#  significativos = which(pvalues < alpha)
#  length(significativos)
#  
#  significativos_BH = which(pBH < alpha)
#  length(significativos_BH)
#  
#  significativos_Bonferroni =  which(pBonferroni < alpha)
#  length(significativos_Bonferroni)

## ---- eval = doAll------------------------------------------------------------
#  cabeceras = c(colnames(rowData(PRJNA548871_keep)), "p.value", "p.adjust_BH", "p.adjust_Bonferroni")
#  ensembl_url = asierortega::ensembl2url(rowData(PRJNA548871_keep)[, "ENSEMBL"])
#  entrzid_url = asierortega::entrezid2url(rowData(PRJNA548871_keep)[, "ENTREZID"])
#  
#  ## Construimos dataframe
#  df51 = data.frame(rowData(PRJNA548871_keep)[, "ENSEMBL"], rowData(PRJNA548871_keep)[, "GENENAME"], rowData(PRJNA548871_keep)[, "ENTREZID"], rowData(PRJNA548871_keep)[, "SYMBOL"], pvalues, pBH, pBonferroni)
#  df51_hip = data.frame(ensembl_url, rowData(PRJNA548871_keep)[, "GENENAME"], entrzid_url, rowData(PRJNA548871_keep)[, "SYMBOL"], pvalues, pBH, pBonferroni)
#  
#  colnames(df51) = cabeceras
#  colnames(df51_hip) = cabeceras
#  head(df51)

## ---- eval = doAll, warning = FALSE-------------------------------------------
#  df51_funcion = sumex2dfDESeq2(PRJNA548871, hiperlink = FALSE)
#  df51_hip_funcion = sumex2dfDESeq2(PRJNA548871, hiperlink = TRUE)

## ---- eval=doAll--------------------------------------------------------------
#  makeReport(df51_hip, name = "Informe_ED_DESeq2_PRJNA548871", path = "../docs/articles")

## ---- eval = doAll------------------------------------------------------------
#  df51_sorted = df51[order(df51$p.value, decreasing = FALSE),]
#  head(df51_sorted)

## ---- eval = doAll------------------------------------------------------------
#  vsd = vst(dds, blind = FALSE)
#  ## Seleccionemos los 304 genes con la mayor varianza entre las muestras
#  topVarGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), 304)
#  mat = assay(vsd)[topVarGenes,]
#  mat = mat - rowMeans(mat)
#  df_anot = as.data.frame(colData(vsd)[, c("treatment")])
#  pheatmap(mat, annotation_col = df_anot, show_rownames = FALSE)

