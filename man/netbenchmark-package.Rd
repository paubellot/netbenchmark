\name{netbenchmark-package}
\alias{netbenchmark-package}
\docType{package}
\title{Benchmarking of several inference networks methods}
\description{
  For a given list of network inference algorithms, \code{netbenchmark} performs a benchmark 
  between them. It makes use of five different big gene datasources, it relies on a random subsampling of 
  each one of the datasource and noise addition in order to generate the datasets. 
  This package is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.
}
\author{
Pau Bellot, Catharina Olsen, Patrick E Meyer, with contributions from Alexandre Irrthum

Maintainer: Pau Bellot <pau.bellot@upc.edu>
}

\references{
Pau Bellot, Catharina Olsen, Philippe Salembier, Albert Oliveras-Verges, 
and Patrick E Meyer. NetBenchmark: A Bioconductor Package for Reproducible
Benchmarks of Gene Regulatory Network Inference. Submitted, 2015.
}

\keyword{ package }
 
\examples{
    \dontrun{
        AUPR20.list<-netbenchmark(datasources.names=c("syntren300",
            "rogers1000"),datasets.num=7)
        AUPR20.300exp.list<-netbenchmark(methods=c("aracne.wrap","mrnet.wrap",
            "GeneNet.wrap"),datasources.names=c("syntren300","rogers1000"),
        experiments=300,global.noise=10,noiseType="lognormal")
        AUPR20.n30.list<-netbenchmark(methods=c("all.fast","mrnet.wrap",
            "Genie3.wrap"),eval="AUROC",local.noise=30)
    }
}
