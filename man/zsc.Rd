\name{zsc}
\alias{zsc}
\alias{netbenchmark_zsc}
\title{
    Z-score c++ function
}
\description{
    Z-score c++ function
}
\usage{
    zsc(x)
}
\arguments{
    \item{x}{Numeric matrix with the microarray dataset to infer the network. 
    Columns contain variables and rows contain samples.}
}
\value{
    A matrix of numerics with the inferred adjacency matrix.
}
\references{
    Prill, Robert J., et al. "Towards a rigorous assessment of systems 
    biology models: the DREAM3 challenges." PloS one 5.2 (2010): e9202.
}
\seealso{
    \code{\link{netbenchmark}}
}
\examples{
    # Data
    data <- runif(100)
    dim(data) <- c(10,10)
    # Inference
    net <- zsc(data)
}
\keyword{ misc }
