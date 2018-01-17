\name{datasource.subsample}
\alias{datasource.subsample}
\title{Subsampling datasource procedure}
\description{
    \code{datasource.subsample} picks randomly the specified amount of 
    samples from the original datasource and also adds noise to the 
    subsampled dataset if it is specified.
}
\usage{
    datasource.subsample(datasource,experiments=NA,datasets.num=5,
        local.noise=20,global.noise=0,noiseType="normal",
        samplevar=TRUE, seed = NULL)
}
\arguments{
    \item{datasource}{data.frame where columns contain variables and 
    rows contain experiments.} 
    \item{experiments}{Integer specifying the number of experiments 
    that for performing the subsampling of datasources (default: NA).}
    \item{datasets.num}{Integer specifying the number of datasets to 
    be generated for each of the selected original datasources 
    (default: 5).}
    \item{local.noise}{Integer specifying the desired percentage of 
    local noise to be added at each of the subsampled datasets 
    (default: 20).}
    \item{global.noise}{Integer specifying the desired percentage of 
    global noise to be added at each of the subsampled datasets 
    (default: 0).}
    \item{noiseType}{Character specifying the type of the noise to be 
    added: "normal" (default: "normal").}
    \item{samplevar}{Logical specifying if the datasets should have 
    variability in the number of experiments between them 
    (default: TRUE).}
    \item{seed}{A single value, interpreted as an integer to specify seeds, 
    useful for creating simulations that can be reproduced 
    (default: \code{NULL}) - see \code{\link[base]{set.seed}}.}	
}
\value{
    \code{datasource.subsample} returns a list with 
    \code{datasets.num} 
    elements, each one of objects contains a data.frame of the 
    subsampled dataset with the amount of Gaussian noise 
    specified that would contain the same number of variables. 
}
\details{
    If the argument \code{experiments} is NA, the value 
    \code{experiments} will be calculated automatically in order to 
    have \code{datasets.num} smaller datasets that does not have the 
    same experiment twice inside each dataset.
    Each of the subsampled datasets \code{experiments} would have a 
    number of experiments around \code{experiments} \eqn{ \pm 20 \% } 
    that would be chosen randomly among the original the original 
    number of experiments without replacement. 

    If the argument \code{experiments} is a number, the number of 
    \code{datasets.num} is calculated automatically.
    If the number of specified \code{experiments} is greater or equal 
    than the original number of experiments, then only a replicate 
    will be generated and the subsampled dataset would have the same 
    dimensions as the original one but the experiments will be 
    unsorted randomly.
    
    
    Two different types of noises could be added, that are specified 
    with the argument \code{noiseType}:
    \itemize{
        \item "local": the variance of the noise is different for each 
        variable and it is the percentage
        specified of the variance of each variable 
        ( \eqn{ \pm 20 \% }).
        \item "Globlal": the variance of the noise is the same for the 
        whole datasource, it is the percentage specified of the mean 
        variance of all the variables ( \eqn{ \pm 20 \% }).
    }
}
\seealso{ \code{\link{netbenchmark}} }

\examples{
    # Subsample
    data.list.1 <- datasource.subsample(syntren300.data)
    data.list.2 <- datasource.subsample(syntren300.data,
        local.noise=10)
    # Inference
    inf.net.1 <- cor(data.list.1[[1]])
    inf.net.2 <- cor(data.list.2[[4]])
}
\keyword{misc}

