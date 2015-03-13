\name{comp.metr}
\alias{comp.metr}
\alias{pr.compute}
\alias{roc.compute}
\alias{fscore}
\alias{aupr}
\alias{auroc}
\alias{pr.plot}
\alias{mean.pr}
\alias{roc.plot}
\title{Compute metrics}
\description{
    A group of functions to plot precision-recall and ROC curves and 
    to compute f-scores from the matrix returned by the 
    \code{\link{evaluate}} function.
}
\usage{
    fscore(table, beta=1)
    auroc(table,k=-1)
    aupr(table,k=-1)
    pr.plot(table,device=-1,...)
    roc.plot(table,device=-1,...)
}
\arguments{
    \item{table}{This is the matrix returned by the \code{evaluate} 
    function where columns contain the confusion matrix 
    \code{TP,FP,TN,FN} values. - see \code{\link{evaluate}}.}
    \item{beta}{Numeric used as the weight of the recall in the 
    f-score formula - see details. The default value of this argument 
    is -1, meaning precision as important as recall.}
    \item{k}{Numeric used as the index to compute the area under the 
    curve until that point- see details. The default value of this 
    argument is -1, meaning that the whole area under the curve is 
    computed}
    \item{device}{The device to be used. This parameter allows the 
    user to plot precision-recall and receiver operating 
    characteristic curves for various inference algorithms on the 
    same plotting window - see examples.}
    \item{...}{ Arguments passed to \code{plot}.}
}
\value{
    The function \code{roc.plot} (\code{pr.plot}) plots the ROC-curve 
    (PR-curve) and returns the device associated with the plotting 
    window.
    
    The function \code{auroc} (\code{aupr}) computes the area under 
    the ROC-curve (PR-curve) using the trapezoidal approximation 
    until point k. 
    
    The function \code{fscore} returns fscores according to the 
    confusion matrices contained in the 'table' argument - see details.
}

\details{
    A confusion matrix contains FP,TP,FN,FP values.
    \itemize{
        \item "true positive rate" 
        \eqn{tpr = \frac{TP}{TN+TP}}{tpr = TP/(TN+TP)}
        \item "false positive rate" 
        \eqn{fpr = \frac{FP}{FN+FP}}{fpr = FP/(FN+FP)}
        \item "precision" \eqn{p = \frac{TP}{FP+TP}}{p = TP/(FP+TP)}
        \item "recall" \eqn{r = \frac{TP}{TP+FN}}{r = TP/(TP+FN)}
        \item "f-beta-score" 
        \eqn{F_\beta = (1+\beta) \frac{p r} {r + \beta p}}
        {Fbeta = (1+beta) * p*r/(r + beta*p)}
    }
}
\seealso{\code{\link{evaluate}}, \code{\link{plot}}}
\examples{
    # Inference
    Net <- cor(syntren300.data)
    # Validation
    tbl <-  evaluate(Net,syntren300.net)
    # Plot PR-Curves
    max(fscore(tbl))
    dev <- pr.plot(tbl, col="green", type="l")
    aupr(tbl)
    idx <- which.max(fscore(tbl))
}
\keyword{misc}