#' URL's from entrez identifiers
#' @description  
#' It makes the url's from entrez identifiers
#' @param id ENTREZ identifiers
#' @return 
#' The corresponding URL's in the database 
#' @export
#' @family URL generation
entrezid2url = function(id)
  ifelse(id == "NA",NA,
         paste("<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
               id,"'>",id,"</a>",sep=""))

#' URL's from Ensembl identifiers
#' @description  
#' It makes the url's from Ensembl identifiers
#' @param id Ensembl identifiers
#' @return 
#' The corresponding URL's in the database 
#' @export
#' @family URL generation
ensembl2url = function(id)
  ifelse(id == "NA",NA,
         paste("<a href='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
               id,"'>",id,"</a>",sep=""))

#' Obtain a data.frame from expression set with p.values using rowttests or limma.
#' @description  
#' Generates a dataframe with the fData information of the expression set
#' and adds the p-value and adjusted p-value information using rowttests of
#' genefilter or limma library.
#' @name eset2df
#' @return 
#' data.frame 
#' @usage 
#' eset2df(eset, method = 1)
#' @param eset Is the ExpressionSet object
#' @param method 1 To determinate the p-values with genefilter library, and 2 to determinate the p-values with limma.
#' @export
#' @family data.frame generation
eset2df = function(eset, method=1){
           if(method==1){
           library(genefilter)
           tt = genefilter::rowttests(eset, pData(eset)[,1])
           pBH = p.adjust(tt$p.value, method = "BH")
           pBonferroni = p.adjust(tt$p.value, method = "bonferroni")
           cabeceras = c(colnames(fData(eset)), "p.value", "p.adjust_BH", "p.adjust_Bonferroni")
           df1 = cbind(fData(eset), tt$p.value, pBH, pBonferroni)
           colnames(df1) = cabeceras
           return(df1)}
          else{
            library(limma)
            library(stringr)
            pData(eset)$Tissue = str_replace(pData(eset)$Tissue, " ", "_") ## En caso de que una variable tenga espacios
            design = model.matrix(~0 + pData(eset)[,1])
            colnames(design) = c(as.character(unique(pData(eset)[,1])[1]), as.character(unique(pData(eset)[,1])[2]))  
            fit = lmFit(eset, design)
            contrast.matrix = makeContrasts(dif = Liver - Gluteus_medius, levels = design)
            fit2 = contrasts.fit(fit, contrast.matrix)
            fit2 = eBayes(fit2)
            p_value_limma = fit2$p.value
            p_value_limma_df=data.frame(p_value_limma)
            pBH_limma = p.adjust(p_value_limma, method = "BH")
            pBonferroni_limma = p.adjust(p_value_limma, method = "bonferroni")
            cabeceras = c(colnames(fData(eset)), "p.value", "p.adjust_BH", "p.adjust_Bonferroni")
            df1_limma = cbind(fData(eset), p_value_limma_df$dif, pBH_limma, pBonferroni_limma)
            return(df1_limma)}
}

#' Remove the rows with NA of the SummarizeExperiment.
#' @description  
#' Removes rows that have as many NAs in the identifiers as the specified value by argument.
#' The number of identifiers of the SummarizedExperiment should be introduced as an argument,
#' in this way we can remove the rows that only have NA values in the different identifiers. 
#' @name cleanNA
#' @return 
#' SummarizedExperiment 
#' @usage 
#' cleanNA(sumex, Ncolumns = 3)
#' @param sumex is the SummarizeExperiment
#' @param Ncolumns The number of identifiers of the SummarizedExperiment.
#' @export
#' @family NA clean
cleanNA = function(sumex, Ncolumns=3){
        library(SummarizedExperiment)
        sumex = sumex[which(apply(rowData(sumex),1,function(y) sum(is.na(y))!=Ncolumns)),]
        return(sumex)
}

#' Obtain a data.frame from summarized experiment with p.values using edgeR.
#' @description  
#' Generates a dataframe with the rowData information of the SummarizeExperiment
#' and adds the p-value and adjusted p-value information using edgeR.
#' The column names of the rowData in the Summarized Experiment must contain the information
#' at least of ENSEMBL, GENENAME, SYMBOL and ENTREZID. 
#' @name sumex2df
#' @return 
#' data.frame 
#' @usage 
#' sumex2df(sumex, columnName="treatment", method="c", hiperlink=FALSE)
#' @param sumex Is the SummarizeExperiment object
#' @param columnName The name of the column with the information of the groups (treatment and control).
#' Default value will be "treatment".
#' @param method is the function that will use to obtain the p-values. 
#' "c" to apply the function of Common Dispersion.
#' "t" to apply the function of Tagwise Dispersion. 
#' As default value will apply the Common Dispersion function.
#' @param hiperlink TRUE will include the corresponding hyperlinks to the ENSEMBL and ENTREZID identifiers. 
#' FALSE to just include the identifiers without the corresponding hyperlinks. 
#' @export
#' @family dataframe generation
sumex2df = function(sumex, columnName="treatment", method="c", hiperlink=FALSE){
    library(SummarizedExperiment)
    library(edgeR)
    library(asierortega)
    Ldge = DGEList(counts = assay(sumex), group = colData(sumex)[, columnName])
    keep = rowSums(cpm(Ldge)>1) >= 2
    Ldge = Ldge[keep, keep.lib.sizes = FALSE]
    sumex = sumex[keep,]
    Ldge.c = estimateCommonDisp(Ldge) 
    if(method=="c"){
      et = exactTest(Ldge.c)
      cabeceras = c(colnames(rowData(sumex)), "p.value_CommonDisp", "p.adjust_BH", "p.adjust_Bonferroni")
    }
    else{
      Ldge.t =  estimateTagwiseDisp(Ldge.c)
      et = exactTest(Ldge.t)
      cabeceras = c(colnames(rowData(sumex)), "p.value_TagwiseDisp", "p.adjust_BH", "p.adjust_Bonferroni")
    }
      pvalues = et$table$PValue
      pBH = p.adjust(pvalues, method = "BH")
      pBonferroni = p.adjust(pvalues, method = "bonferroni")
    if(hiperlink==TRUE){
      ensembl = asierortega::ensembl2url(rowData(sumex)[, "ENSEMBL"])
      entrezid = asierortega::entrezid2url(rowData(sumex)[, "ENTREZID"])
    }
    else{
      ensembl = rowData(sumex)[, "ENSEMBL"]
      entrezid = rowData(sumex)[, "ENTREZID"]
    }
      df = data.frame(ensembl, rowData(sumex)[, "GENENAME"], entrezid, rowData(sumex)[, "SYMBOL"], pvalues, pBH, pBonferroni)
      colnames(df) = cabeceras
      return(df)
  }

#' Generates an HTML report of a given data.frame.
#' @description  
#'Generates an HTML report of a given given data.frame in a exact path given as argument. 
#' @name makeReport
#' @return 
#' Generates an HTML report.
#' @usage 
#' makeInform(df, name, path=".")
#' @param df Is the dataframe
#' @param name Name of the inform/file
#' @param path Is the path where you want to generate the HTML inform. As default value will be current working directory.
#' @export
#' @family report
makeReport = function(df, name, path="."){
  library(ReportingTools)
  report = HTMLReport(shortName = name, title = name, reportDirectory = path)
  publish(df, report)
  finish(report)
}

#' Obtain a data.frame from summarized experiment with p.values using DESeq2.
#' @description  
#' Generates a dataframe with the rowData information of the SummarizeExperiment
#' and adds the p-value and adjusted p-value information using DESeq2 differential
#' expression analysis. This function is specially built for the following dataset: PRJNA548871.
#' The column names of the rowData in the Summarized Experiment must contain the information
#' at least of ENSEMBL, GENENAME, SYMBOL and ENTREZID. 
#' @name sumex2dfDESeq2
#' @return 
#' data.frame 
#' @usage 
#' sumex2df(sumex, hiperlink=FALSE)
#' @param sumex Is the SummarizeExperiment object
#' @param hiperlink TRUE will include the corresponding hyperlinks to the ENSEMBL and ENTREZID identifiers. 
#' FALSE to just include the identifiers without the corresponding hyperlinks. 
#' @export
#' @family dataframe generation
sumex2dfDESeq2 = function(sumex, hiperlink=FALSE){
    library(SummarizedExperiment)
    library(DESeq2)
    library(asierortega)
    dds =  DESeqDataSet(se = PRJNA548871, design = ~ treatment)
    keep = rowSums(counts(dds)) >= 10
    dds = dds[keep,]
    sumex = sumex[keep,]
    dds$treatment <- relevel(dds$treatment, ref  = "nonOIR")
    dds = DESeq(dds)
    res = lfcShrink(dds, coef=2)
    pvalues = res$pvalue
    pBH = p.adjust(pvalues, method = "BH")
    pBonferroni = p.adjust(pvalues, method = "bonferroni")
    cabeceras = c(colnames(rowData(sumex)), "p.value", "p.adjust_BH", "p.adjust_Bonferroni")
    if(hiperlink==TRUE){
      ensembl = asierortega::ensembl2url(rowData(sumex)[, "ENSEMBL"])
      entrezid = asierortega::entrezid2url(rowData(sumex)[, "ENTREZID"])
    }
    else{
      ensembl = rowData(sumex)[, "ENSEMBL"]
      entrezid = rowData(sumex)[, "ENTREZID"]
    }
    df = data.frame(ensembl, rowData(sumex)[, "GENENAME"], entrezid, rowData(sumex)[, "SYMBOL"], pvalues, pBH, pBonferroni)
    colnames(df) = cabeceras
    return(df)
}

#' Unify the probes thar correspond to the same gene.
#' @description  
#' Unify the probes thar correspond to the same gene.
#' It unify all the probes to correspond to the same gene identifier that is 
#' introduced by argument, as default value will be ENTREZID. 
#' @name unifyeset
#' @return 
#' ExpressionSet
#' @usage 
#' unifyeset(eset, ID = ENTREZID)
#' @param eset is the ExpressionSet
#' @param ID Is the identifier that the probe have been annotated.
#' @export
#' @family NA clean
unifyeset = function(eset, ID="ENTREZID"){
  eset = eset[which(!is.na(fData(eset)[,ID])),]
  seleccion = match(unique(fData(eset)[,ID]),fData(eset)[,ID])
  eset = eset[seleccion,]
  return(eset)
}

#' URL's from KEGG pathway identifiers
#' @description  
#' It makes the url's from KEGG pathway identifiers
#' @param id KEGG pathway identifiers
#' @return 
#' The corresponding URL's in the database 
#' @export
#' @family URL generation
kegg2url = function(id)
  ifelse(id == "NA",NA,
         paste("<a href='https://www.genome.jp/kegg-bin/show_pathway?",
               id,"'>",id,"</a>",sep=""))

#' Obtain a data.frame from summarized experiment and group of genes using Fisher's test.
#' @description  
#' Generates a dataframe with the gsRanking information of the SummarizeExperiment and the group of genes. 
#' The SummarizedExperiment must have in the colData a column that is called 'GROUP' and must have '0' and '1' values.
#' Must have ENTRZID identifiers.
#' @name overex2df
#' @return 
#' data.frame 
#' @usage 
#' overex2df(sumex, hiperlink=FALSE)
#' @param sumex Is the SummarizeExperiment object.
#' @param grpgenes Group of genes generated  getGenesets function.
#' @param hiperlink TRUE will include the corresponding hyperlinks to KEGG. 
#' @export
#' @family dataframe generation
overex2df = function(sumex, grpgenes, hiperlink=FALSE){
  library(EnrichmentBrowser)
  sumex = deAna(expr = sumex) 
  sumex.oraKEGG= sbea(method = "ora", se = sumex, gs = grpgenes,
                         perm = 0, alpha = 0.05)
  cabeceras = c(colnames(gsRanking(sumex.oraKEGG)), "p.adjust_BH", "p.adjust_Bonferroni")
  if(hiperlink==TRUE){
    kegg_url = asierortega::kegg2url(gsRanking(sumex.oraKEGG)[, "GENE.SET"])
  }
  else{
    kegg_url = gsRanking(sumex.oraKEGG)[, "GENE.SET"]
  }
  df = data.frame(kegg_url, gsRanking(sumex.oraKEGG)[, "NR.GENES"],  
                 gsRanking(sumex.oraKEGG)[, "NR.SIG.GENES"],  
                 gsRanking(sumex.oraKEGG)[, "PVAL"], pBH, pBonferroni)
  colnames(df) = cabeceras
  return(df)
}