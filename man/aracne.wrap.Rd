\name{aracne.wrap}
\alias{aracne.wrap}
\title{aracne wrapper function}
\description{
    Default wrapper function for the aracne network inference algorithm
}
\usage{
    aracne.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network.Columns contain variables and rows contain samples.}
}
\details{
    The motivation of the Algorithm for the Reconstruction of 
    Accurate Cellular NEtworks (ARACNE) is that many similar measures 
    between variables may be the result of indirect effects. In order 
    to delete the indirect effect the algorithm relies on the 
    ``Data Processing Inequality'', this process removes the weakest link 
    in every triplet of connected variables.
}
\value{
    \code{aracne.wrapper} returns a matrix which is the weighted adjacency 
    matrix of the network inferred by aracne algorithm.
    The wrapper uses the "spearman" correlation 
    (can be used with continuous data) to estimate the 
    entropy - see \code{\link[minet]{build.mim}}
}
\references{
    Margolin, Adam A., et al. "ARACNE: an algorithm for the reconstruction 
    of gene regulatory networks in a mammalian cellular context." 
    BMC Bioinformatics 7.Suppl 1 (2006): S7.

    Meyer, Patrick E., Frederic Lafitte, and Gianluca Bontempi. 
    "minet: AR/Bioconductor package for inferring large transcriptional 
    networks using mutual information." BMC Bioinformatics 9.1 (2008): 461.
}  
\seealso{
    \code{\link{netbenchmark}}, \code{\link{evaluate}},
    \code{\link[minet]{aracne}}
}
\examples{
    # Data
    data <- runif(100)
    dim(data) <- c(10,10)
    # Inference
    net <- aracne.wrap(as.data.frame(data))
}
\keyword{ misc }
