## -----------------------------------------------------------------------------
doAll = FALSE

## ---- eval = doAll------------------------------------------------------------
#  library(pacman)
#  library(asierortega)
#  library(tami)
#  pacman::p_load(ggplot2)
#  pacman::p_load(Biobase)
#  pacman::p_load(SummarizedExperiment)

## ---- eval = doAll------------------------------------------------------------
#  data(gse115484, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  gse115484 = asierortega::unifyeset(gse115484, ID = "ENTREZID")

## ---- eval = doAll------------------------------------------------------------
#  data(sscKEGGgsc, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
#  gse115484_self_mean = GeneSetTest(x = gse115484, y = "Tissue",
#           test = rowtmod, association = "statistic", correction = "BH",
#           GeneNullDistr = "randomization", GeneSetNullDistr = "self-contained",
#           alternative = "two-sided",nmax = 100, id = "ENTREZID",
#           gsc = sscKEGGgsc, descriptive = mean, foutput = "gse115484_self_mean")

## ---- eval = doAll------------------------------------------------------------
#  gse115484_self_maxmean = GeneSetTest(x = gse115484, y = "Tissue",
#           test = rowtmod, association = "statistic", correction = "BH",
#           GeneNullDistr = "randomization", GeneSetNullDistr = "self-contained",
#           alternative = "two-sided",nmax = 100, id = "ENTREZID",
#           gsc = sscKEGGgsc, descriptive = maxmean, foutput = "gse115484_self_maxmean")

## ---- eval = doAll------------------------------------------------------------
#  gse115484_self_median = GeneSetTest(x = gse115484, y = "Tissue",
#           test = rowtmod, association = "statistic", correction = "BH",
#           GeneNullDistr = "randomization", GeneSetNullDistr = "self-contained",
#           alternative = "two-sided",nmax = 100, id = "ENTREZID",
#           gsc = sscKEGGgsc, descriptive = median, foutput = "gse115484_self_median")

## ---- eval = doAll------------------------------------------------------------
#  set.seed(123)
#  gse115484_comp_mean = GeneSetTest(x = gse115484, y = "Tissue",
#           test = rowtmod, association = "statistic", correction = "BH",
#           GeneNullDistr = "randomization", GeneSetNullDistr = "competitive",
#           alternative = "two-sided",nmax = 100, id = "ENTREZID",
#           gsc = sscKEGGgsc, descriptive = mean, foutput = "gse115484_comp_mean")
#  
#  set.seed(123)
#  gse115484_comp_maxmean = GeneSetTest(x = gse115484, y = "Tissue",
#           test = rowtmod, association = "statistic", correction = "BH",
#           GeneNullDistr = "randomization", GeneSetNullDistr = "competitive",
#           alternative = "two-sided",nmax = 100, id = "ENTREZID",
#           gsc = sscKEGGgsc, descriptive = maxmean, foutput = "gse115484_comp_maxmean")
#  
#  set.seed(123)
#  gse115484_comp_median = GeneSetTest(x = gse115484, y = "Tissue",
#           test = rowtmod, association = "statistic", correction = "BH",
#           GeneNullDistr = "randomization", GeneSetNullDistr = "competitive",
#           alternative = "two-sided",nmax = 100, id = "ENTREZID",
#           gsc = sscKEGGgsc, descriptive = median, foutput = "gse115484_comp_median")

## ---- eval = doAll------------------------------------------------------------
#  db_name = "KEGG"
#  gse115484_self_mean_df = tidy(gse115484_self_mean)
#  gse115484_comp_mean_df = tidy(gse115484_comp_mean)
#  gse115484_self_maxmean_df = tidy(gse115484_self_maxmean)
#  gse115484_comp_maxmean_df = tidy(gse115484_comp_maxmean)
#  gse115484_self_median_df = tidy(gse115484_self_median)
#  gse115484_comp_median_df = tidy(gse115484_comp_median)
#  
#  colnames(gse115484_self_mean_df)[1] = db_name
#  colnames(gse115484_comp_mean_df)[1] = db_name
#  colnames(gse115484_self_maxmean_df)[1] = db_name
#  colnames(gse115484_comp_maxmean_df)[1] = db_name
#  colnames(gse115484_self_median_df)[1] = db_name
#  colnames(gse115484_comp_median_df)[1] = db_name

## ---- eval = doAll------------------------------------------------------------
#  alpha = 0.05
#  length(gse115484_comp_maxmean_df$adjp)
#  significativos_p_self_mean = which(gse115484_self_mean_df$adjp < alpha)
#  length(significativos_p_self_mean)
#  
#  significativos_p_comp_mean = which(gse115484_comp_mean_df$adjp < alpha)
#  length(significativos_p_comp_mean)
#  
#  significativos_p_self_maxmean = which(gse115484_self_maxmean_df$adjp < alpha)
#  length(significativos_p_self_maxmean)
#  
#  significativos_p_comp_maxmean = which(gse115484_comp_maxmean_df$adjp < alpha)
#  length(significativos_p_comp_maxmean)
#  
#  significativos_p_self_median = which(gse115484_self_median_df$adjp < alpha)
#  length(significativos_p_self_median)
#  
#  significativos_p_comp_median = which(gse115484_comp_median_df$adjp < alpha)
#  length(significativos_p_comp_median)

