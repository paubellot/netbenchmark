\name{mutrank.wrap}
\alias{mutrank.wrap}
\title{
    Mutual Rank wrapper function
}
\description{
    A wrapper function for mutual rank.
}
\usage{
    mutrank.wrap(data)
}
\arguments{
    \item{data}{Numeric matrix with the microarray dataset to infer the 
    network. Columns contain variables and rows contain samples.}
}
\value{
    \code{mutrank.wrap} returns a matrix which is the weighted adjacency 
    matrix of the network inferred by Mutual Rank algorithm.
}
\references{
    Obayashi, Takeshi, and Kengo Kinoshita. "Rank of correlation coefficient 
    as a comparable measure for biological significance of gene 
    coexpression." DNA research 16.5 (2009): 249-260.
}
\examples{
     # Data
    data <- grndata::getData(datasource.name = "toy",FALSE)
    # Inference
    net <- mutrank.wrap(data)
}
\keyword{ misc }
