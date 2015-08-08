\name{zscore.wrap}
\alias{zscore.wrap}
\title{
    Zscore wrapper function
}
\description{
    Z-score wrapper function.
}
\usage{
    zscore.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network. Columns contain variables and rows contain samples.}
}
\details{
    Zscore is a method that assumes interventional data, more concretely 
    knockout experiments that leads to a change in other genes. The 
    assumption  is that the knocked-out gene \eqn{i} in the experiment 
    \eqn{k} affects  more strongly to the genes that it regulates than  
    the others, the effect of the gene \eqn{i} over the gene \eqn{j} 
    is captured with the Zscore \eqn{z_{ij}}:
    \deqn{ 
        z_{ij}=|\frac{x_{jk}-\mu_j}{\sigma_j}|
    }
    \eqn{mu_j} and \eqn{\sigma_j} are respectively the mean and standard 
    deviation of the empirical distribution of the gene \eqn{j}.
}
\value{
    \code{zscore.wrap} returns a matrix which is the weighted adjacency 
    matrix of the network inferred by Zscore algorithm.
}
\references{
    Prill, Robert J., et al. "Towards a rigorous assessment of systems 
    biology models: the DREAM3 challenges." PloS one 5.2 (2010): e9202.
}
\examples{
    # Data
    data <- grndata::getData(datasource.name = "toy",FALSE)
    # Inference
    net <- zscore.wrap(data)
}
\keyword{ misc }
