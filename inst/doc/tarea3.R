## -----------------------------------------------------------------------------
doAll = TRUE

## ---- eval = doAll------------------------------------------------------------
library(pacman)
pacman::p_load(Biobase)
pacman::p_load(genefilter)
pacman::p_load(limma)
library(asierortega)

## ---- eval = doAll------------------------------------------------------------
data(gse115484, package = "asierortega")

## ---- eval = doAll------------------------------------------------------------
gse115484 = gse115484[which(!is.na(fData(gse115484)[,"ENTREZID"])),]
seleccion = match(unique(fData(gse115484)[,"ENTREZID"]),fData(gse115484)[,"ENTREZID"])
gse115484 = gse115484[seleccion,]
## Utilizando la funcion
gse115484_funcion = unifyeset(gse115484, ID = "ENTREZID")

## ---- eval= doAll-------------------------------------------------------------
tissue = pData(gse115484)[,"Tissue"]
tt = genefilter::rowttests(gse115484, tissue)
head(tt)
p_values = tt$p.value

## ----eval = doAll-------------------------------------------------------------
pBH = p.adjust(p_values, method = "BH")
pBonferroni = p.adjust(p_values, method = "bonferroni")

## ---- eval = doAll------------------------------------------------------------
alpha = 0.05
significativos_p = which(p_values < alpha)
length(significativos_p)

significativos_BH = which(pBH < alpha)
length(significativos_BH)

significativos_Bonferroni =  which(pBonferroni < alpha)
length(significativos_Bonferroni)

## ---- eval=doAll--------------------------------------------------------------
design = model.matrix(~0 + tissue)
colnames(design) = c("Liver", "Gluteus_medius")
head(design)

## ---- eval = doAll------------------------------------------------------------
fit = lmFit(gse115484, design)
(contrast.matrix = makeContrasts(dif = Liver - Gluteus_medius, levels = design))

## ---- eval = doAll------------------------------------------------------------
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

## ---- eval = doAll------------------------------------------------------------
p_value_limma = fit2$p.value
pBH_limma = p.adjust(p_value_limma, method = "BH")
pBonferroni_limma = p.adjust(p_value_limma, method = "bonferroni")

## ---- eval = doAll------------------------------------------------------------
topTable(fit2, coef=1, adjust="BH")

## ---- eval = doAll------------------------------------------------------------
topTable(fit2, coef=1, adjust="bonferroni")

## ---- eval = doAll------------------------------------------------------------
cabeceras = c(colnames(fData(gse115484)), "p.value", "p.adjust_BH", "p.adjust_Bonferroni")
df1 = cbind(fData(gse115484), tt$p.value, pBH, pBonferroni)
colnames(df1) = cabeceras
head(df1)

## ---- eval=doAll--------------------------------------------------------------
df1_prueba = eset2df(gse115484, method = 1)
head(df1_prueba)
df1_prueba_limma = eset2df(gse115484, method = 2)
head(df1_prueba_limma)

## ---- eval = doAll------------------------------------------------------------
df3 = topTable(fit2, coef=1, adjust="BH", number = length(fit2$coefficients))
df3_sorted = df3[order(df3$adj.P.Val, decreasing = TRUE),]
(df3_sorted_filtrado = na.omit(df3_sorted))
lista_genes = df3_sorted_filtrado$SYMBOL[1:200]
write.table(df3_sorted_filtrado$SYMBOL[1:200], file="salida_genes.txt", row.names = FALSE, col.names = FALSE)

