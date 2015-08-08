\name{clr.wrap}
\alias{clr.wrap}
\title{CLR wrapper function}
\description{
    Default wrapper function for the CLR network inference algorithm
}
\usage{
    clr.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer 
    the network. Columns contain variables and rows contain samples.}
}
\details{
    The Context Likelihood or Relatedness network (CLR) method derives 
    a score that is associated to the empirical distribution of the 
    mutual information  values, in practice the score between gene 
    \eqn{X_i}{Xi} and gene \eqn{X_j}{Xj} is defined as follows 
    \eqn{z_{ij}=\sqrt{z^2_i+z^2_j}}{z_{ij}=\sqrt{z^2_i+z^2_j}}, where: 
    \deqn{
        z_i=max \left( 0,\frac{I(X_i;X_j)-\mu_i}{\sigma_i} \right)
    }{z_i=max( 0,(I(Xi;Xj)-\mu_i)/\sigma_i )}
    
    \eqn{\mu_i}{\mu_i} and \eqn{\sigma_i}{\sigma_i} are 
    respectively the mean and standard deviation of the empirical 
    distribution of the mutual information between both genes.
}
\value{
    \code{clr.wrap} returns a matrix which is the weighted adjacency 
    matrix of the network inferred by CLR algorithm.
    The wrapper uses the "spearman" correlation 
    (can be used with continuous data) to estimate the 
    entropy - see \code{\link[minet]{build.mim}}. 
}
\references{
    Faith, Jeremiah J., et al. "Large-scale mapping and validation of 
    Escherichia coli transcriptional regulation from a compendium of 
    expression profiles." PLoS biology 5.1 (2007): e8.
}  

\seealso{ \code{\link{netbenchmark}},  
\code{\link{evaluate}}, \code{\link[minet]{clr}}}
\examples{
    # Data
    data <- grndata::getData(datasource.name = "toy",FALSE)
    # Inference
    net <- clr.wrap(data)
}
\keyword{ misc }
