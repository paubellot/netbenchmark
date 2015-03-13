\name{rate}
\alias{rate}
\alias{netbenchmark_rate}
\title{Inference Validation}
\description{
    \code{rate} compares the infered network to the true underlying network 
    for all the sorted predictions provided and appends the resulting 
    confusion matrices to the returned object.
}
\usage{
    rate(PredEdgeList, GSEdgeList, ngenes, sym)
}
\arguments{
    \item{PredEdgeList}{The inferred network in the form of a EdgeList.}
    \item{GSEdgeList}{The true underlying in the form of a EdgeList.}
    \item{ngenes}{Integrer denoting the number of total genes 
    in the network.}
    \item{sym}{Logical specifying if the evaluation is symmetric 
    (default: TRUE) - see \code{\link{evaluate}}.}
}
\value{
    A matrix of numerics with the contingency table for each link in 
    \code{PredEdgeList}.
}
\author{
    Pau Bellot
}
\seealso{ 
    \code{\link{netbenchmark}}, \code{\link{evaluate}}, 
    \code{\link{comp.metr}}
} 
\examples{
    # Data
    net <- matrix(0,10,10)
    net[sample(1:100,20)] <- 1
    # Simulated Inference
    inf <- net+matrix(rnorm(100,sd=0.5),10,10)
    table <- evaluate(inf,net)
}
\keyword{ misc }

