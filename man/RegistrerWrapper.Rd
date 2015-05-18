\name{RegisterWrapper}
\alias{RegisterWrapper}
\alias{UnregisterWrapper}
\title{
    Wrapper (un)registration routine
}
\description{
    These function allows the registration and unregistration of a wrapper     
    function to the all.fast or all methods of \code{\link{netbenchmark}}.
    After registring it \code{wrapper.name} function will belong 
    to all.fast or all methods during the R session. 
    Unregistring the \code{wrapper.name} function will remove it from
    all.fast or all methods during the R session.
}
\usage{
    RegisterWrapper(wrapper.name=NULL,all.fast=TRUE)
    UnregisterWrapper(wrapper.name=NULL,all.fast=TRUE)
}
\arguments{
  \item{wrapper.name}{The character (vector) of wrapper names 
  (default: NULL).}
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
    \dontrun{
        # Register it to all.fast methods
        RegisterWrapper("Spearmancor")
        # Register it to all methods
        RegisterWrapper("Spearmancor", all.fast=FALSE)
        # Unregister it from all.fast methods
        UnregisterWrapper("Spearmancor")
        # Unregister it from all methods
        UnregisterWrapper("Spearmancor", all.fast=FALSE)
    }
    
}
\keyword{misc}
