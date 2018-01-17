\name{evaluate}
\alias{evaluate}
\title{Inference Evaluation}
\description{
    \code{evaluate} compares the inferred network to the true 
    underlying network for several threshold values and appends the 
    resulting confusion matrices to the returned object.     
}
\usage{
    evaluate(inf.net,true.net,sym=TRUE,extend=0)
}
\arguments{
    \item{inf.net}{An adjacency matrix representing the inferred 
    network.}
    \item{true.net}{An adjacency matrix representing the true 
    underlying network.}
    \item{sym}{Logical, make a symmetric evaluation 
    (default = TRUE).}
    \item{extend}{Integer, specifying the desired number of links to 
    extend in the network (default=0)}
}
\value{
    \code{evaluate} returns a matrix with four columns representing 
    \code{TP,FP,TN,FN}.
    These values are computed for each of the predicted links that 
    should be sorted. Thus, each row of the returned object contains 
    the confusion matrix as a function of the cutoff in the edgelist.      
}
\details{
    The first edgelist network \code{inet} is compared to the true 
    underlying network, \code{tnet}, in order to compute the 
    metrics of the performance.
    If extend is specified, \code{extend} links that network 
    \code{inet} has set to 0 are added to the inferred network 
    randomly at the end of the edgelist.
}
\seealso{ \code{\link{netbenchmark}} }
\examples{
    # Inference
    inf.net <- cor(syntren300.data)
    #Evaluate
    table <- evaluate(inf.net, syntren300.net)
    table.nosym <- evaluate(inf.net, syntren300.net,sym=FALSE)
}
\keyword{misc}
