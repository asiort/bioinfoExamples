## -----------------------------------------------------------------------------
doAll = FALSE

## ---- eval = doAll------------------------------------------------------------
#  library(pacman)
#  library(asierortega)
#  pacman::p_load(SummarizedExperiment)
#  pacman::p_load(edgeR)
#  pacman::p_load(ggplot2)
#  pacman::p_load(ReportingTools)

## ---- eval = doAll------------------------------------------------------------
#  data(PRJNA548871, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  PRJNA548871 = PRJNA548871[which(apply(rowData(PRJNA548871),1,function(y) sum(is.na(y))!=3)),]
#  PRJNA548871_funcion = cleanNA(PRJNA548871, Ncolumns = 3)

## ---- eval = doAll------------------------------------------------------------
#  Ldge = DGEList(counts = assay(PRJNA548871), group = colData(PRJNA548871)[, "treatment"])
#  attributes(Ldge)
#  head(Ldge$counts)
#  class(Ldge$counts)
#  head(Ldge$samples)
#  class(Ldge$samples)

## ---- eval = doAll, warning = FALSE-------------------------------------------
#  ggplot(Ldge$samples, aes(x = lib.size)) + geom_dotplot()

## ---- eval = doAll------------------------------------------------------------
#  keep = rowSums(cpm(Ldge)>1) >= 2
#  Ldge = Ldge[keep, keep.lib.sizes = FALSE]
#  PRJNA548871_keep = PRJNA548871[keep,]

## ---- eval = doAll------------------------------------------------------------
#  ## Estimación dispersión común
#  Ldge.c = estimateCommonDisp(Ldge)
#  ## Estimación dispersión por gen
#  Ldge.t =  estimateTagwiseDisp(Ldge.c)
#  ## Aplicamos exactTest
#  et.c = exactTest(Ldge.c)
#  et.t = exactTest(Ldge.t)

## ---- eval = doAll------------------------------------------------------------
#  topTags(et.c)
#  topTags(et.t)

## ---- eval = doAll------------------------------------------------------------
#  pvalues_c = et.c$table$PValue
#  pvalues_t =et.t$table$PValue
#  pBH_c = p.adjust(pvalues_c, method = "BH")
#  pBH_t = p.adjust(pvalues_t, method = "BH")
#  pBonferroni_c = p.adjust(pvalues_c, method = "bonferroni")
#  pBonferroni_t = p.adjust(pvalues_t, method = "bonferroni")

## ---- eval = doAll------------------------------------------------------------
#  length(pvalues_t)
#  alpha = 0.05
#  significativos_p_t = which(pvalues_t < alpha)
#  length(significativos_p_t)
#  
#  significativos_BH_t = which(pBH_t < alpha)
#  length(significativos_BH_t)
#  
#  significativos_Bonferroni_t =  which(pBonferroni_t < alpha)
#  length(significativos_Bonferroni_t)

## ---- eval = doAll------------------------------------------------------------
#  cabeceras_c = c(colnames(rowData(PRJNA548871_keep)), "p.value_CommonDisp", "p.adjust_BH", "p.adjust_Bonferroni")
#  ensembl_url = asierortega::ensembl2url(rowData(PRJNA548871_keep)[, "ENSEMBL"])
#  entrzid_url = asierortega::entrezid2url(rowData(PRJNA548871_keep)[, "ENTREZID"])
#  ## Construimos dataframe
#  df51_c = data.frame(rowData(PRJNA548871_keep)[, "ENSEMBL"], rowData(PRJNA548871_keep)[, "GENENAME"], rowData(PRJNA548871_keep)[, "ENTREZID"], rowData(PRJNA548871_keep)[, "SYMBOL"], pvalues_c, pBH_c, pBonferroni_c)
#  df51_c_hip = data.frame(ensembl_url, rowData(PRJNA548871_keep)[, "GENENAME"], entrzid_url, rowData(PRJNA548871_keep)[, "SYMBOL"], pvalues_c, pBH_c, pBonferroni_c)
#  
#  colnames(df51_c) = cabeceras_c
#  colnames(df51_c_hip) = cabeceras_c
#  head(df51_c)

## ---- eval = doAll------------------------------------------------------------
#  ## Sin hipervinculo
#  cabeceras_t = c(colnames(rowData(PRJNA548871_keep)), "p.value_TagwiseDisp", "p.adjust_BH", "p.adjust_Bonferroni")
#  ensembl_url = asierortega::ensembl2url(rowData(PRJNA548871_keep)[, "ENSEMBL"])
#  entrzid_url = asierortega::entrezid2url(rowData(PRJNA548871_keep)[, "ENTREZID"])
#  ## Construimos data.frame
#  df51_t = data.frame(rowData(PRJNA548871_keep)[, "ENSEMBL"], rowData(PRJNA548871_keep)[, "GENENAME"], rowData(PRJNA548871_keep)[, "ENTREZID"], rowData(PRJNA548871_keep)[, "SYMBOL"], pvalues_t, pBH_t, pBonferroni_t)
#  df51_t_hip = data.frame(ensembl_url, rowData(PRJNA548871_keep)[, "GENENAME"], entrzid_url, rowData(PRJNA548871_keep)[, "SYMBOL"], pvalues_t, pBH_t, pBonferroni_t)
#  
#  colnames(df51_t) = cabeceras_t
#  colnames(df51_t_hip) = cabeceras_t
#  head(df51_t)

## ---- eval = doAll------------------------------------------------------------
#  df51_c_sorted = df51_c[order(df51_c$p.value_CommonDisp, decreasing = FALSE),]
#  df51_t_sorted = df51_t[order(df51_t$p.value_TagwiseDisp, decreasing = FALSE),]
#  
#  head(df51_c_sorted)
#  head(df51_t_sorted)

## ---- eval = doAll------------------------------------------------------------
#  df51_c_funcion = sumex2df(PRJNA548871_keep, columnName = "treatment", method = "c", hiperlink = FALSE)
#  df51_t_funcion = sumex2df(PRJNA548871_keep, columnName = "treatment", method = "t", hiperlink = FALSE)
#  df51_c_hip_funcion = sumex2df(PRJNA548871_keep, columnName = "treatment", method = "c", hiperlink = TRUE)
#  df51_t_hip_funcion = sumex2df(PRJNA548871_keep, columnName = "treatment", method = "c", hiperlink = TRUE)

## ---- eval=doAll--------------------------------------------------------------
#  nombre_c = "Informe_ED_CommonDisp_PRJNA548871"
#  nombre_t = "Informe_ED_TagwiseDisp_PRJNA548871"
#  
#  
#  report_c = HTMLReport(shortName = nombre_c, title = nombre_c, reportDirectory = "../docs/articles")
#  report_t = HTMLReport(shortName = nombre_t, title = nombre_t, reportDirectory = "../docs/articles")
#  
#  publish(df51_c_hip, report_c)
#  publish(df51_t_hip, report_t)
#  
#  finish(report_c)
#  finish(report_t)

## ---- eval = doAll------------------------------------------------------------
#  makeReport(df51_c_funcion, name = "Reportdeprueba", path = "../docs/articles")

## ---- eval = doAll------------------------------------------------------------
#  data(PRJNA548871, package = "asierortega")
#  PRJNA548871 = PRJNA548871[which(apply(rowData(PRJNA548871),1,function(y) sum(is.na(y))!=3)),]
#  
#  counts0 = assay(PRJNA548871)
#  counts0 = data.frame(counts0)
#  df0 = reshape2::melt(counts0)
#  ini_fdens = ggplot(df0, aes(x=value, colour = variable)) + stat_ecdf()
#  ini_fdist = ggplot(df0, aes(x = value,colour = variable))+geom_density()
#  
#  ## Normalizamos con el metodo TMM
#  Ldge_norm = calcNormFactors(Ldge, method = "TMM")
#  Ldeg_norm = estimateCommonDisp(Ldge_norm)
#  
#  counts1 = Ldeg_norm$pseudo.counts
#  counts1 = data.frame(counts1)
#  df0 = reshape2::melt(counts1)
#  pos_fdens = ggplot(df0, aes(x = value,colour = variable)) + stat_ecdf()
#  pos_fdist = ggplot(df0, aes(x = value,colour = variable)) + geom_density()

## ---- eval = doAll------------------------------------------------------------
#  summary(Ldge$samples[,"norm.factors"])
#  summary(Ldge_norm$samples[,"norm.factors"])

## ---- eval = doAll------------------------------------------------------------
#  ## funcion de densidad de los conteos inicial
#  ini_fdens
#  ## funcion de distribucion de los conteos inicial
#  ini_fdist
#  ## funcion de densidad de los conteos nromalizado con TMM
#  pos_fdens
#  ## funcion de distribucion de los conteos normalizado con TMM
#  pos_fdens

