\name{mrnetb.wrap}
\alias{mrnetb.wrap}
\title{mrnetb wrapper function}
\description{
    Default wrapper function for the MRNETB network inference algorithm
}
\usage{
    mrnetb.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network. Columns contain variables and rows contain samples.}
}
\details{
    \code{\link[minet]{mrnetb}} takes the mutual information matrix as input 
    in order to infer the network using the maximum relevance/minimum 
    redundancy criterion combined with a backward elimination and a 
    sequential replacement - see references. 
    This method is a variant of mrnet.
}
\value{
    \code{mrnetb.wrap} returns a matrix which is the weighted adjacency 
    matrix of the network inferred by mrnetb algorithm.
    The wrapper uses the "spearman" correlation 
    (can be used with continuous data) 
    to estimate the entropy - see \code{\link[minet]{build.mim}}.
}

\references{
    Patrick E. Meyer, Daniel Marbach, Sushmita Roy and Manolis Kellis. 
    Information-Theoretic Inference of Gene Networks Using 
    Backward Elimination. The 2010 International Conference on 
    Bioinformatics and Computational Biology. 

    Patrick E. Meyer, Kevin Kontos, Frederic Lafitte and Gianluca Bontempi. 
    Information-theoretic inference of large transcriptional regulatory
    networks. EURASIP Journal on Bioinformatics and Systems Biology, 2007.
}
\seealso{ 
    \code{\link{netbenchmark}},
    \code{\link{evaluate}}, \code{\link[minet]{mrnet}}
}
\examples{
    # Data
    data <- grndata::getData(datasource.name = "toy",FALSE)
    # Inference
    net <- mrnetb.wrap(data)
}
\keyword{ misc }
