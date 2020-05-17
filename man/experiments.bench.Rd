\name{experiments.bench}
\alias{experiments.bench}
\title{
    Noise sensitivity benchmark
}
\description{
    For a given vector of character of the names of wrapper functions that 
    compute a network inference methods, \code{experiments.bench} performs a 
    number of experiments sensitivity test. 
    It makes use of five different big gene datasets subsampling them to 
    generate different datasets.num of the network with different number 
    of experiments.
}
\usage{
    experiments.bench(methods = "all.fast", datasources.names = "all",
        experiments = c(20, 50, 150), eval = "AUPR",
        no.topedges = 20, datasets.num = 3, local.noise = 20,
        global.noise = 0, noiseType = "normal", sym = TRUE,
        seed = NULL, verbose= TRUE)
}
\arguments{
    \item{methods}{A vector of characters containing the names of network 
    inference algorithms wrappers to be compared (default: "all.fast").}
    \item{datasources.names}{A vector of characters containing the names of 
    network datasets to be included in the benchmark (default: "all").}
    \item{experiments}{A vector to set the number of experiments to test 
    the methods (default=c(20,50,150)).}
    \item{eval}{The name of the evaluation metric among the following ones: 
    "no.truepos", "AUROC" or "AUPR" (default : "AUPR") 
    - see \code{\link{evaluate}}.}
    \item{no.topedges}{Float specifying the percentage  number of links to be 
    considered in the evaluation (default: 20).}
    \item{datasets.num}{Number of repetitions in the noise evaluation, for 
    each method and each dataset and each noise intensity (default: 3).}
    \item{local.noise}{Integer specifying the desired percentage of local 
    noise to be added at each of the subsampled datasets (default: 20) 
    - see \code{\link{datasource.subsample}}.}
    \item{global.noise}{Integer specifying the desired percentage of global 
    noise to be added at each of the subsampled datasets (default: 20) 
    - see \code{\link{datasource.subsample}}.}
    \item{noiseType}{Character specifying the type of the noise to be added: 
    "normal" or "lognormal" (default: "normal") 
    - see \code{\link{datasource.subsample}}.}
    \item{sym}{Logical specifying if the evaluation is symmetric 
    (default: TRUE) - see \code{\link{evaluate}}.}
    \item{seed}{A single value, interpreted as an integer to specify seeds, 
    useful for creating simulations that can be reproduced 
    (default: \code{NULL}) - see \code{\link[base]{set.seed}}.}
    \item{verbose}{Logical specifying if the code should provide a log 
    about what the function is doing (default: TRUE).}
}
\details{
    The argument \code{methods} accepts "all.fast" and "all" 
    (case insensitive) as a parameters:
    \itemize{
        \item "all.fast" performs network inference with "aracne", "c3net", 
        "clr", "GeneNet", "mutual ranking", "mrnet"
        \item "all" performs network inference with "aracne", "c3net", 
        "clr", "GeneNet", "Genie3", "mutual ranking", "mrnet", "mrnetb"    
        } 
    It evaluates the first \code{no.topedges} \% of the possible links 
    inferred by each algorithm at each dataset.

    Two different types of noises are added independently:
    \itemize{
        \item "Local": the standard deviation of the noise is different 
        for each variable.  
        \code{local.noise} specifies the percentage for each 
        variable ( \eqn{ \pm 20 \% }).
        \item "Global": the standard deviation of the noise is the same 
        for the whole dataset. 
        \code{global.noise} specifies the percentage of the mean standard 
        deviation of all the variables ( \eqn{ \pm 20 \% }).
    }
    The distribution of noise is set with \code{noiseType}, it is possible 
    to choose between "normal" (\code{\link[stats]{rnorm}}) and "lognormal" 
    (\code{\link[stats]{rlnorm}}). The argument \code{noiseType} can be a 
    single character, this specifies the same distribution for both "Local" 
    and "Global" noise, it also can be a vector of characters with two 
    elements, the former specifies the distribution of "Local" noise and 
    the later the distribution of "Global" noise.
}
\value{
    \code{experiments.bench} returns a list with three elements:
    \enumerate{
        \item A data.frame which is the result table containing the number of 
        true positives as an evaluation measure. 
        It evaluates each algorithm specified at \code{methods} at each one 
        of the specified \code{datasources.names} with different noise 
        intensities.
        \item A data.frame which is the corresponding pvalue table of the 
        corresponding statistical test 
        for each one of the \code{datasets.num} between the best algorithm 
        and the others. 
        \item The seed of the random number generators that allows the 
        replication of the results.
    }
}
\author{
    Pau Bellot and Patrick Meyer
}
\seealso{ 
    \code{\link{netbenchmark}}, \code{\link{noise.bench}}
}
\examples{
    results <- experiments.bench(datasources.names="toy",
    datasets.num=2,methods="all.fast",experiments=c(20,40))
}
\keyword{ misc }