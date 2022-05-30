## ----echo=FALSE---------------------------------------------------------------
doAll = FALSE

## ---- eval=doAll--------------------------------------------------------------
#  library(pacman)
#  pacman::p_load(Rsamtools)
#  pacman::p_load(GenomicFeatures)
#  pacman::p_load(GenomicAlignments)
#  pacman::p_load(AnnotationDbi)
#  pacman::p_load(org.Mm.eg.db)

## ---- eval=doAll--------------------------------------------------------------
#  setwd("../../bamfiles")
#  sampleTable = read.table("bamfiles.txt")
#  dirActual =  paste(getwd(), "/", sep="")
#  fls = paste(dirActual, sampleTable[,1], sep = "")
#  bamList = BamFileList(fls, index = character(), yieldSize = 100000, obeyQname = TRUE)
#  gtfFile = "Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
#  txdb = makeTxDbFromGFF(gtfFile, format = "gtf")
#  genes = exonsBy(txdb, by="gene")
#  PRJNA548871 = summarizeOverlaps(features = genes, read=bamList,
#      mode="Union",
#      singleEnd=TRUE,    ## Porque no son lecturas apareadas
#      ignore.strand=TRUE,
#      fragments=FALSE)

## ---- eval=doAll--------------------------------------------------------------
#  SampleName = paste("GSM38909", 42:57, sep = "")
#  Run = paste("SRR93023", 11:26, sep = "")
#  treatment = factor(c(rep("OIR", 8), rep("nonOIR", 8)))
#  colData(PRJNA548871) = DataFrame(SampleName,Run , treatment)
#  metadata(PRJNA548871)=list("Experimenter name"="Clemens Lange",
#                             "Organization"="Uniklinik Freiburg",
#                             "Title"="Temporo-spatial distribution and transcriptional profile of retinal microglia in the oxygen-induced retinopathy mouse model",
#                             "URL"="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132731",
#                             "Abstract"= " Myeloid cells such as resident retinal microglia or infiltrating blood-derived macrophages accumulate in areas of retinal ischemia and neovascularization (RNV) and modulate neovascular eye disease. Their temporo-spatial distribution and biological function in this process, however, remain unclarified. We determined the extent of microglia proliferation and macrophage infiltration in areas of ischemia and RNV using Cx3cr1CreERT2:Rosa26-tdTomato mice and assessed the transcriptional profile of microglia in the oxygen-induced retinopathy (OIR) mouse model. We show that microglia are the predominant myeloid cell population in areas of RNV. Thirty percent of retinal microglia were EdU-positive indicating considerable expansion of local microglia. RNA-Seq revealed an enrichment of processes related to cell division and chemokine receptor binding. We propose that activated retinal MG alter their transcriptional profile, exhibit proliferative ability and are by far the most frequent myeloid cell population in areas of RNV in the OIR model thus presenting a potential target for future therapeutic approaches.")

## ---- eval=doAll--------------------------------------------------------------
#  columns(org.Mm.eg.db)
#  a = AnnotationDbi::select(org.Mm.eg.db,keys=rownames(PRJNA548871),
#                            columns = c("GENENAME","ENTREZID","ENSEMBL", "SYMBOL"), keytype = "ENSEMBL")
#  b = match(rownames(PRJNA548871),a[,"ENSEMBL"])
#  rowData(PRJNA548871) = a[b,]
#  PRJNA548871 = PRJNA548871[which(!is.na(rowData(PRJNA548871)[,"ENSEMBL"])),]
#  seleccion = match(unique(rowData(PRJNA548871)[,"ENSEMBL"]), rowData(PRJNA548871)[,"ENSEMBL"])
#  PRJNA548871 = PRJNA548871[seleccion,]
#  PRJNA548871
#  assay(PRJNA548871)

## ---- eval=doAll--------------------------------------------------------------
#  nullsum = apply(assay(PRJNA548871), 1, sum) == 0
#  PRJNA548871 = PRJNA548871[!nullsum,]
#  PRJNA548871
#  assay(PRJNA548871)

## ---- eval=doAll--------------------------------------------------------------
#  save(PRJNA548871, file = "PRJNA548871.rda")

