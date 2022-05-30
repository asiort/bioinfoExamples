#' @title GSE115484
#' 
#' @description
#' Sus scrofa
#' 
#' Expression data from Microarray experiment in muscle and liver tissues
#' 
#' Expression profiling by array.
#' 
#'In the current work we have compared the transcriptomic and eQTL profiles of the porcine skeletal muscle and
#'liver by using a data set of 103 Duroc pigs genotyped with the Illumina SNP60 BeadChip and with available 
#'microarray measurements of gene expression for both tissues.' 
#' @source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115484}
#'  
#'  \url{https://pubmed.ncbi.nlm.nih.gov/31234802/}
#' @examples
#' data(gse115484, package="asierortega")
#' @docType data
#' @keywords datasets
#' @format SummarizedExperiment
#' 
#' @name GSE115484
NULL
#' @title PRJNA548871
#' 
#' @description
#' Mus musculus
#' 
#' SummarizedExperiment data from RNA-seq experiment in Mus musculus retinal microglia tissue.
#' 
#' Temporo-spatial distribution and transcriptional profile of retinal microglia in the oxygen-induced retinopathy mouse model.
#' 
#'RNA-Seq revealed an enrichment of processes related to cell division and chemokine receptor binding. We propose that activated 
#retinal MG alter their transcriptional profile, exhibit proliferative ability and are by far the most frequent myeloid cell population 
#in areas of RNV in the OIR model thus presenting a potential target for future therapeutic approaches. ' 
#' @source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132731}
#'  
#'  \url{https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA548871}
#' @examples
#' data(PRJNA548871, package="asierortega")
#' @docType data
#' @keywords datasets
#' @format SummarizedExperiment
#' 
#' @name PRJNA548871
NULL
#' @title GSE115484_DE_analysis
#' 
#' @description
#' Sus scrofa
#' 
#' Expression data from Microarray experiment in muscle and liver tissues.
#' 
#' Expression profiling by array.
#' 
#'Differential expression analysis with the data generated in GSE115484 using genefilter and limma methods, tarea3.Rmd. 
#' @source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115484}
#'  
#'  \url{https://pubmed.ncbi.nlm.nih.gov/31234802/}
#' @docType data
#' @keywords datasets
#' @format analysis
#' 
#' @name GSE115484_DE_analysis
NULL
#' @title PRJNA548871_DE_analysis_edgeR
#' 
#' @description
#' Mus musculus
#' 
#' Expression profiling by high throughput sequencing
#' 
#' SummarizedExperiment data from RNA-seq experiment in Mus musculus retinal microglia tissue.
#' 
#'Differential expression edgeR classic analysis with the data generated in PRJNA548871, tarea41.Rmd.' 
#' @source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132731}
#'  
#' @docType data
#' @keywords datasets
#' @format analysis
#' 
#' @name PRJNA548871_DE_analysis_edgeR
NULL
#' @title PRJNA548871_DE_analysis_DESeq2
#' 
#' @description
#' Mus musculus
#' 
#' Expression profiling by high throughput sequencing
#' 
#' SummarizedExperiment data from RNA-seq experiment in Mus musculus retinal microglia tissue.
#' 
#'Differential expression DESeq2 analysis with the data generated in PRJNA548871, tarea42.Rmd.' 
#' @source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132731}
#'  
#' @docType data
#' @keywords datasets
#' @format analysis
#' 
#' @name PRJNA548871_DE_analysis_DESeq2
NULL
#' @title GeneSet_Analysis_GSE115484
#' 
#' @description
#' Sus scrofa
#' 
#' Expression data from Microarray experiment in muscle and liver tissues.
#' 
#' Expression profiling by array.
#' 
#'Gene set analysis with the data generated in GSE115484 using EnrichmentBrowser package (tarea51.Rmd) and GeneSetTest function from tami (tarea52.Rmd). 
#' @source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115484}
#'  
#'  \url{https://pubmed.ncbi.nlm.nih.gov/31234802/}
#' @docType data
#' @keywords datasets
#' @format analysis
#' 
#' @name GeneSet_Analysis_GSE115484
NULL
#' @title GeneSet_Analysis_PRJNA548871
#' 
#' @description
#' Mus musculus
#' 
#' Expression profiling by high throughput sequencing
#' 
#' SummarizedExperiment data from RNA-seq experiment in Mus musculus retinal microglia tissue.
#' 
#'Gene set analysis with the data generated in GSE115484 using EnrichmentBrowser package (tarea51.Rmd).
#' @source
#'  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132731}
#'  
#' @docType data
#' @keywords datasets
#' @format analysis
#' 
#' @name GeneSet_Analysis_PRJNA548871
#' NULL