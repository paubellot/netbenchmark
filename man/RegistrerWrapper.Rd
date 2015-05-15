\name{RegistrerWrapper}
\alias{RegistrerWrapper}
\alias{UnregistrerWrapper}
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
    RegistrerWrapper(wrapper.name=NULL,all.fast=TRUE)
    UnregistrerWrapper(wrapper.name=NULL,all.fast=TRUE)
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
        RegistrerWrapper("Spearmancor")
        # Register it to all methods
        RegistrerWrapper("Spearmancor", all.fast=FALSE)
        # Unregister it from all.fast methods
        UnregistrerWrapper("Spearmancor")
        # Unregister it from all methods
        UnregistrerWrapper("Spearmancor", all.fast=FALSE)
    }
    
}
\keyword{misc}
