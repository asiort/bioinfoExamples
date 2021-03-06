% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{sumex2dfDESeq2}
\alias{sumex2dfDESeq2}
\title{Obtain a data.frame from summarized experiment with p.values using DESeq2.}
\usage{
sumex2df(sumex, hiperlink=FALSE)
}
\arguments{
\item{sumex}{Is the SummarizeExperiment object}

\item{hiperlink}{TRUE will include the corresponding hyperlinks to the ENSEMBL and ENTREZID identifiers. 
FALSE to just include the identifiers without the corresponding hyperlinks.}
}
\value{
data.frame
}
\description{
Generates a dataframe with the rowData information of the SummarizeExperiment
and adds the p-value and adjusted p-value information using DESeq2 differential
expression analysis. This function is specially built for the following dataset: PRJNA548871.
The column names of the rowData in the Summarized Experiment must contain the information
at least of ENSEMBL, GENENAME, SYMBOL and ENTREZID.
}
\seealso{
Other dataframe generation: 
\code{\link{overex2df}()},
\code{\link{sumex2df}()}
}
\concept{dataframe generation}
