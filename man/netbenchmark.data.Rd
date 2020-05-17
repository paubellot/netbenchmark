\name{netbenchmark.data}
\alias{netbenchmark.data}
\title{Benchmarking of several network inference algorithms for your own data}
\usage{
netbenchmark.data(methods = "all.fast", data = NULL, true.net = NULL,
  eval = "AUPR", no.topedges = 20, sym = TRUE, plot = FALSE,
  verbose = TRUE)
}
\arguments{
\item{methods}{A vector of characters containing the names of network i
inference algorithms wrappers to be compared (default: "all.fast").}

\item{data}{data.frame containing the data. Each row should contain a
microarray experiment and each column a gene (default: NULL).}

\item{true.net}{matrix containg underlying network in the form of adjacency
matrix (default: NULL).}

\item{eval}{The name of the evaluation metric among the following ones:
"no.truepos", "AUROC" or "AUPR" (default : "AUPR").}

\item{no.topedges}{Float specifying the percentage number of links to be
considered in the evaluation (default: 20)}

\item{sym}{Logical specifying if the evaluation is symmetric (default: TRUE)
- see \code{\link{evaluate}}}

\item{plot}{(default: FALSE)}

\item{verbose}{Logical specifying if the code should provide a log about
what the function is doing (default: TRUE).}
}
\description{
Benchmarking of several network inference algorithms for your own data
}
\value{ 
    \code{netbenchmark.data} returns a list with three elements. 
    \enumerate{
        \item A data.frame which is the result table of the selected measure.
        \item A data.frame which contains the CPU Time Used (in seconds) by 
        the algorithm to infer the network.
        \item A list containing the mean precision recall curves of the 
        different algorithms for each datasource.
    }
    Each of these data.frame will have the same number of columns as methods 
    provided by the user and an additional one for a random method.
}
\details{
    The argument \code{methods} accepts "all.fast" and "all" 
    (case insensitive) as a parameters:
    \itemize{
        \item "all.fast" performs network inference with "aracne", "c3net", 
        "clr", "GeneNet", "mutual ranking", "mrnet" (and registered
        methods with \code{\link{RegisterWrapper}}.)
        \item "all" performs network inference with "aracne", "c3net", "clr", 
        "GeneNet", "Genie3", "mutual ranking", "mrnet", "mrnetb" (and 
        registered methods with \code{\link{RegisterWrapper}}.)
    } 
    All the measures only evaluates the first \code{no.topedges} \% of the 
    possible links inferred by each algorithm at each dataset.
}
\author{
    Pau Bellot, Catharina Olsen and Patrick E Meyer 
    Maintainer: Pau Bellot <pau.bellot@upc.edu>
}
\seealso{ 
    \code{\link{netbenchmark}}, \code{\link{evaluate}},
    \code{\link{comp.metr}}
}    
\examples{
    Data <- grndata::getData(datasource.name="toy")
    top20.aupr <- netbenchmark.data(methods="all",data = Data[[1]],
                               true.net= Data[[2]])
}
