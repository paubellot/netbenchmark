\name{RegistrerWrapper}
\alias{RegistrerWrapper}

\title{
    Wrapper registration routine
}
\description{
    This function allows the registration of a wrapper function to the
    all.fast or all methods of \code{\link{netbenchmark}}.
    After registring it \code{wrapper.name} function will belong 
    to all.fast or all methods during the R session. 
}
\usage{
    RegistrerWrapper(wrapper.name=NULL,all.fast=TRUE)
}
\arguments{
  \item{wrapper.name}{The character (vector) of wrapper names (default: NULL).}
  \item{all.fast}{Logical indicating if the \code{wrapper.name} should be 
  added to all.fast or all methods (default: TRUE).}
}

\value{
    Displays a message if the registration could be performed or not.
}

\author{
    Pau Bellot, Catharina Olsen and Patrick E Meyer 
    Maintainer: Pau Bellot <pau.bellot@upc.edu>
}

\seealso{
    \code{\link{netbenchmark}}
}

\examples{
    # Define a wrapper function
    Spearmancor <- function(data){
        cor(data,method="spearman")
    }
    # Register it to all.fast methods
    RegistrerWrapper("Spearmancor")
    # Register it to all methods
    RegistrerWrapper("Spearmancor", all.fast=FALSE)
}
\keyword{misc}
