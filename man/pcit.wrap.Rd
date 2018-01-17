\name{pcit.wrap}
\alias{pcit.wrap}
\title{pcit wrapper function}
\description{
    Default wrapper function for the pcit network inference algorithm
}
\usage{
    pcit.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network. Columns contain variables and rows contain samples.}
}
\details{
    The Partial Correlation coefficient with Information
    Theory (PCIT) algorithm, combines the concept of
    partial correlation coefficient with information theory to identify 
    significant gene-to-gene associations. 

    For every trio of genes in \eqn{X_i}, \eqn{X_j} and \eqn{X_l}, the three 
    first-order partial correlation 
    coefficients are computed. These coefficients
    indicate the strength of the linear relationship between \eqn{X_i} and
    \eqn{X_j} that is uncorrelated with \eqn{X_l}, being therefore a measure
    of conditional independence. Then, the average ratio of partial to direct 
    correlation is computed in order to obtain the tolerance level to be used
    as the local threshold for eliminating non-significant associations.
}
\value{
    \code{pcit.wrap} returns a matrix which is the weighted adjacency matrix 
    of the network inferred by pcit algorithm.
}
\references{
    Reverter, Antonio, and Eva KF Chan. "Combining partial 
    correlation and an information theory approach to the reversed engineering 
    of gene co-expression networks." Bioinformatics 24.21 (2008): 2491-2497.
}
\seealso{
    \code{\link{netbenchmark}}, \code{\link{evaluate}}, 
    \code{\link[PCIT]{pcit}}
}
\examples{
    # Data
    data <- grndata::getData(datasource.name = "toy",FALSE)
    # Inference
    net <- pcit.wrap(data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
