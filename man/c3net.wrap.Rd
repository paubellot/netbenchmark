\name{c3net.wrap}
\alias{c3net.wrap}
\title{c3net wrapper function}
\description{
    Default wrapper function for the C3net network 
    inference algorithm
}
\usage{
    c3net.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer 
    the network. Columns contain variables and rows contain samples.}
}
\details{
    The Conservative Causal Core NETwork (C3NET) consists of two main steps. 
    The first step is  the elimination of non-significant edges,
    and the second step selects for each gene the edge among the
    remaining ones with maximum mutual information value. C3NET does not aim 
    at inferring the entire network underlying gene regulation but mainly 
    tries to recover the core structure.
}
\value{
    \code{c3net.wrap} returns a matrix which is the weighted adjacency 
    matrix of the network inferred by c3net algorithm. 
    The Mutual Information threshold is set to 0 
    - see \code{\link[c3net]{c3net}}.
}
\references{
    Altay, G\"{o}kmen, and Frank Emmert-Streib. "Inferring the conservative 
    causal core of gene regulatory networks." 
    BMC Systems Biology 4.1 (2010): 132.
}
\seealso{
    \code{\link{netbenchmark}},  
    \code{\link{evaluate}}, \code{\link[c3net]{c3net}}
}
\examples{
    # Data
    data <- runif(100)
    dim(data) <- c(10,10)
    # Inference
    net <- c3net.wrap(as.data.frame(data))
}
\keyword{ misc }
