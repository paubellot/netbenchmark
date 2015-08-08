\name{GeneNet.wrap}
\alias{GeneNet.wrap}
\title{GeneNet wrapper function}
\description{
    Default wrapper function for the GeneNet network inference algorithm
}
\usage{
    GeneNet.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network. Columns contain variables and rows contain samples.}
}
\details{
    GeneNEt uses an heuristic for learning statistically a causal network. 
    It relies on a conversion of a network inferred through correlation 
    into a partial correlation graph. Then, a partial ordering of the nodes 
    is assigned by means of a multiple testing of the log-ratio of 
    standardized partial variances. This allows identifying a directed 
    acyclic causal network as a sub-graph of the partial correlation network.
}
\value{
    \code{GeneNet.wrap} The function returns a matrix which is the weighted 
    adjacency matrix of the network inferred by GeneNet algorithm. 
    The shrinkage method used to estimate the partial correlation matrix is 
    "static". - see \code{\link[GeneNet]{ggm.estimate.pcor}}. 
    The probability threshold is set to 0.8. 
    - see \code{\link[GeneNet]{ggm.estimate.pcor}}. 
}
\references{
    Opgen-Rhein, Rainer, and Korbinian Strimmer. "Inferring gene dependency 
    networks from genomic longitudinal data: a functional data approach." 
    RevStat 4.1 (2006): 53-65.

    Opgen-Rhein, Rainer, and Korbinian Strimmer. "Using regularized dynamic 
    correlation to infer gene dependency networks from time-series 
    microarray data." Proceedings of the 4th International Workshop on 
    Computational Systems Biology (WCSB 2006), Tampere. Vol. 4. 2006.

    Sch\"{a}fer, Juliane, and Korbinian Strimmer. "A shrinkage approach 
    to large-scale covariance matrix estimation and implications for 
    functional genomics." Statistical applications in genetics and 
    molecular biology 4.1 (2005): 32.
}

\seealso{ 
    \code{\link{netbenchmark}}, \code{\link{evaluate}}, 
    \code{\link[GeneNet]{GeneNet-package}}
}
\examples{
    # Data
    data <- grndata::getData(datasource.name = "toy",FALSE)
    # Inference
    net <- GeneNet.wrap(data)
}
\keyword{ misc }
