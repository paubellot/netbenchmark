\name{noise.bench}
\alias{noise.bench}
\title{
    Noise sensitivity test
}
\description{
    For a given vector of character of the names of wrapper functions that 
    compute a network inference methods, \code{noise.bench} performs a noise 
    sensitivity test.
    It makes use of different big gene datasets adding Gaussian noise with 
    different intensity to evaluate the performance of the methods. 
}
\usage{
    noise.bench(methods = "all.fast", datasources.names = "all",
        eval = "AUPR", no.topedges = 20,
        datasets.num = 3, local.noise = seq(0, 100, len = 3),
        global.noise = 0, noiseType = "normal", sym = TRUE,
        seed = NULL, verbose = TRUE)
}
\arguments{
    \item{methods}{A vector of characters containing the names of network 
    inference algorithms wrappers to be compared (default: "all.fast").}
    \item{datasources.names}{A vector of characters containing the names of 
    network datasets to be included in the benchmark (default: "all").}
    \item{eval}{The name of the evaluation metric among the following ones: 
    "no.truepos", "AUROC" or "AUPR" (default : "AUPR") 
    - see \code{\link{evaluate}}.}
    \item{datasets.num}{Number of repetitions in the noise evaluation, for 
    each method and each dataset and each noise intensity (default: 5).}
    \item{no.topedges}{Float specifying the percentage  number of links to be 
    considered in the evaluation (default: 20).}
    \item{local.noise}{Vector specifying the desired percentage of local 
    noise to be added at each of the subsampled datasets 
    (default: seq(0, 100, len = 3)).}
    \item{global.noise}{Vector specifying the desired percentage of global 
    noise to be added at each of the subsampled datasets (default: 0).}
    \item{noiseType}{Character specifying the type of the noise to be added:
    "normal" (default: "normal").}  
    \item{sym}{Logical specifying if the evaluation is symmetric
    (default: TRUE) - see \code{\link{evaluate}}.}
    \item{seed}{A single value, interpreted as an integer to specify seeds,
    useful for creating simulations that can be reproduced 
    (default: \code{NULL}) - see \code{\link[base]{set.seed}}.}
    \item{verbose}{Logical specifying if the code should provide a log about
    what the function is doing (default: TRUE).}
}
\details{
    The argument \code{methods} accepts "all.fast" and "all" 
    (case insensitive) as a parameters:
    \itemize{
        \item "all.fast" performs network inference with "aracne", "c3net", 
        "clr", "GeneNet", "mutual ranking", "mrnetb", "pcit"
        \item "all" performs network inference with "aracne", "c3net", "clr", 
        "GeneNet", "Genie3", "mutual ranking", "mrnet", "mrnetb", "pcit" 
    } 
    It evaluates the first \code{no.topedges} \% of the possible links 
    inferred by each algorithm at each dataset.
}
\value{
    \code{noise.bench} returns a list with three elements: 
    \enumerate{
        \item A data.frame which is the result table containing the number of 
        true positives as an evaluation measure. 
        It evaluates each algorithm specified at \code{methods} at each one 
        of the specified \code{datasources.names} with the \code{local.noise} 
        and  \code{global.noise} specified. For each combination the 
        algorithms are evaluated \code{datasets.num} times and their results 
        are averaged.
        \item A data.frame which is the corresponding pvalue table of the 
        corresponding statistical test for each one of the \code{datasets.num} 
        between the best algorithm and the others.
        \item The seed of the random number generators that allows the 
        replication of the results.
    }
}
\author{
    Pau Bellot and Patrick Meyer
}
\seealso{ \code{\link{netbenchmark}}, \code{\link{experiments.bench}}}
\examples{
    results <- noise.bench(datasources.names="toy",
    datasets.num=2,methods="all.fast")
}
\keyword{ misc }
