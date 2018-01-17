\name{mrnet.wrap}
\alias{mrnet.wrap}
\title{mrnet wrapper function}
\description{
    Default function for the MRNET network inference algorithm
}
\usage{
    mrnet.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network. Columns contain variables and rows contain samples.}
}
\details{
    The MRNET approach consists in repeating a MRMR feature selection 
    procedure for each variable of the dataset.
    The MRMR method starts by selecting the variable \eqn{X_i}{Xi} having 
    the highest mutual information with the target \eqn{Y}{Y}. 
    In the following steps, given a set \eqn{\mathcal{S}}{S} of selected 
    variables, the criterion updates \eqn{\mathcal{S}}{S} by choosing the 
    variable \eqn{X_k}{Xk} that maximizes 
    \eqn{ I(X_k;Y) - \frac{1}{|\mathcal{S}|}\sum_{X_i \in \mathcal{S}} 
    I(X_k;X_i)}{ I(Xk;Y) - mean(I(Xk;Xi)), Xi in S.}\cr
    The weight of each pair \eqn{X_i,X_j}{Xi,Xj} will be the maximum score 
    between the one computed when \eqn{X_i}{Xi} is the target and the one 
    computed when \eqn{X_j}{Xj} is the target.
}
\value{
    \code{mrnet.wrap} returns a matrix which is the weighted adjacency  
    matrix of the network inferred by MRNET algorithm. 
    The wrapper uses the "spearman" correlation (can be used with continuous 
    data) to estimate the entropy - see \code{\link[minet]{build.mim}}.
}
\references{
    Patrick E. Meyer, Kevin Kontos, Frederic Lafitte and Gianluca Bontempi. 
    Information-theoretic inference of large transcriptional regulatory 
    networks. EURASIP Journal on Bioinformatics and Systems Biology, 2007.

    Patrick E. Meyer, Frederic Lafitte and Gianluca Bontempi. 
    minet: A R/Bioconductor Package for Inferring Large Transcriptional 
    Networks Using Mutual Information. BMC Bioinformatics, Vol 9, 2008.

    H. Peng, F.long and C.Ding. 
    Feature selection based on mutual information: Criteria of 
    max-dependency, max relevance and min redundancy.
    IEEE transaction on Pattern Analysis and Machine Intelligence, 2005.
}
\seealso{ 
    \code{\link{netbenchmark}},  
    \code{\link{evaluate}}, \code{\link[minet]{mrnet}}
}
\examples{
     # Data
    data <- grndata::getData(datasource.name = "toy",FALSE)
    # Inference
    net <- mrnet.wrap(data)
}
\keyword{ misc }
