\name{Genie3.wrap}
\alias{Genie3.wrap}
\title{Genie3 wrapper function}
\description{Default wrapper function for the Genie3 network 
inference algorithm}
\usage{
    Genie3.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network. Columns contain variables and rows contain samples.}
}
\details{      
    GEne Network Inference with Ensemble of trees (Genie3) algorithm 
    uses the Random Forests feature selection technique to solve 
    a regression problem for each of the genes in the network. 
    In each of the regression problems, the expression pattern of the target
    gene should be predicted from the expression patterns of all 
    transcription factors. 
    The importance of each transcription factor in the prediction of the 
    target gene is taken as an indication of an apparent regulatory link. 
    Then these candidate regulatory links are aggregated over all genes to
    generate a ranking for the whole network.
}
\value{
    \code{Genie3.wrap} returns a matrix which is the weighted adjacency 
    matrix of the network inferred by Genie3 algorithm. 500 trees are used 
    in ensemble for each target gene. 
}
\references{
    Irrthum, Alexandre, Louis Wehenkel, and Pierre Geurts. "Inferring 
    regulatory networks from expression data using tree-based methods." 
    PloS one 5.9 (2010): e12776. 
    
    Breiman, Leo. "Random forests." Machine learning 45.1 (2001): 5-32.
}

\seealso{ \code{\link{netbenchmark}}, \code{\link{evaluate}}}
\examples{
    # Data
    data <- runif(100)
    dim(data) <- c(10,10)
    # Inference
    net <- Genie3.wrap(as.data.frame(data))
}
\keyword{ misc }
