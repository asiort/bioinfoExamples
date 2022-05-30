## -----------------------------------------------------------------------------
doAll = FALSE

## ---- eval=doAll--------------------------------------------------------------
#  library(pacman)
#  pacman::p_load(GEOquery)
#  pacman::p_load(oligo)
#  pacman::p_load(Biobase)
#  pacman::p_load(porcine.db)
#  pacman::p_load(AnnotationDbi)

## ---- eval=doAll--------------------------------------------------------------
#  options(timeout = 360)
#  getGEOSuppFiles("GSE115484")
#  setwd("GSE115484")
#  system("tar -xvf GSE115484_RAW.tar")
#  system("rm GSE115484_RAW.tar")

## ---- eval=doAll--------------------------------------------------------------
#  path = getwd()
#  setwd("GSE115484")
#  cel_files = list.files(getwd(), pattern = "CEL" ,full.names = TRUE)
#  rawData = read.celfiles(cel_files)
#  system("rm *CEL.gz")
#  system("rm -fr ../GSE115484")

## ---- eval=doAll--------------------------------------------------------------
#  exprsMraw = exprs(rawData)
#  dim(exprsMraw)
#  colnames(exprsMraw)
#  row.names(exprsMraw)

## ---- eval=doAll--------------------------------------------------------------
#  MAplot(rawData)

## ---- eval=doAll--------------------------------------------------------------
#  hist(rawData, ylab="Density", xlab = "log Intensity",
#       main = "Estimadores de la densidad de las intensidades previo al preprocesado")

## ---- eval=doAll--------------------------------------------------------------
#  boxplot(rawData, col="white",
#          main="Comparacion Diagramas Cajas previo al preprocesado",
#          ylab="Log2 Intensity", xlab="Samples", names = FALSE)

## ---- eval=doAll--------------------------------------------------------------
#  dataRMA = rma(rawData)
#  exprM = exprs(dataRMA)

## ---- eval=doAll--------------------------------------------------------------
#  MAplot(dataRMA)

## ---- eval=doAll--------------------------------------------------------------
#  hist(dataRMA, ylab="Density", xlab = "log Intensity",
#      main = "Estimadores de la densidad de las intensidades posterior al preprocesado")

## ---- eval=doAll--------------------------------------------------------------
#  boxplot(dataRMA, col="white",
#          main="Comparacion Diagramas Cajas posterior al preprocesado",
#          ylab="Log2 Intensity", xlab="Samples", names = FALSE)

## ---- eval=doAll--------------------------------------------------------------
#  infoData = new('MIAME',
#             name='Mármol-Sánchez E, Gonzalez-Prendes R',
#             lab='Animal Genetics',
#             contact ='Emilio Marmol <emilio.marmol@cragenomica.es>',
#             title = 'Expression data from Microarray experiment in muscle and liver tissues',
#             abstract = 'Summary: In the current work we have compared the transcriptomic and eQTL profiles of the porcine skeletal muscle and            liver by using a data set of 103 Duroc pigs genotyped with the Illumina SNP60 BeadChip and with available microarray measurements            of gene expression for both tissues.
#    	
#             Overall design: 	Gluteus medius and liver samples were collected from 103 Duroc pigs (Lipgen population) after slaughtering and              immediately frozen at -80ºC in liquid nitrogen. Total RNA was extracted from both GM and liver samples and mRNA expression                   profiles were characterized by hybridization to the GeneChip Porcine arrays',
#             url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115484')
#  
#  finput = paste(getwd(), "/../inst/extdata/sample_attributes3.csv", sep = "")
#  pd0 = data.frame(Tissue = factor(c(rep("Gluteus medius", 103), rep("Liver", 103))), Age = rep("190 days", 206))
#  rownames(pd0) = colnames(exprM)
#  metadatos = data.frame(labelDescription = c("Tissue", "Age"), row.names = colnames(pd0))
#  datosfenotipo = new("AnnotatedDataFrame", data=pd0, varMetadata=metadatos)
#  gse115484 = new("ExpressionSet", exprs = exprM, phenoData=datosfenotipo,
#              experimentData=infoData, annotation = "pd.porcine")

## ----eval=doAll---------------------------------------------------------------
#  columns(porcine.db)
#  identificadores = select(porcine.db, keys = featureNames(gse115484),
#                    column = c("ENTREZID","REFSEQ","SYMBOL"), keytype = "PROBEID")
#  identificadores2 = match(featureNames(gse115484), identificadores[,"PROBEID"])
#  fData(gse115484) = identificadores[identificadores2,]
#  fData(gse115484)
#  fData(gse115484)[!is.na(fData(gse115484)[,"ENTREZID"]),]
#  save(gse115484, file = "gse115484.rda")

## ---- eval=doAll--------------------------------------------------------------
#  dim(fData(gse115484))[1]
#  a = dim(fData(gse115484))
#  dim(fData(gse115484)[!is.na(fData(gse115484)[,"ENTREZID"]),])[1]
#  
#  gse115484_limpio = gse115484[which(!is.na(fData(gse115484)[,"ENTREZID"])),]
#  seleccion_unica = match(unique(fData(gse115484_limpio)[,"ENTREZID"]),fData(gse115484_limpio)[,"ENTREZID"])
#  gse115484_limpio = gse115484_limpio[seleccion_unica,]
#  
#  b = dim(fData(gse115484_limpio))
#  (b[1])/(a[1])

