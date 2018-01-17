Genie3.wrap <- function(data){
    net <- .get.weight.matrix(t(data),trace=FALSE,nb.trees=500)
    return(net);
}

# R implementation of GEne Network Inference with Ensemble of Trees (GENIE3)

# Modification history:
#   19 July 2010:   - first version
#   13 August 2010: - fixed formatting of link list table
#       - fixed warnings when K="all"
#       - added choice of importance measure
#       - added possibility to fix random number generator seed
#   5 March 2014:   - Integration with netbenchmark package
#   Author:
#       Alexandre Irrthum
#       alexandre.irrthum@ulg.ac.be


# get.weight.matrix: compute weighted adjacency matrix of inferred network

# Parameters (required):
#   -- expr.matrix: gene expression matrix as returned by read.expr.matrix

# Parameters (optional):
#   -- K: choice of number of input genes randomly selected as candidates  
#       at each node 
#     Must be one of
#       - "sqrt" for square root of total number of input genes (default)
#       - "all" for total number of input genes (minus 1)
#       -  an integer 
#   -- nb.trees: number of trees in ensemble for each target gene 
#       (default 1000)
#     -- input.idx: subset of genes used as input genes (default all genes)
#   Must be either
#       - a vector of indices, e.g. c(1,5,6,7), or
#       - a vector of gene names, e.g. c("at_12377", "at_10912")
#   -- importance.measure: Type of variable importance measure
#   Must be one of
#       - "IncNodePurity" for importance measure based on decrease of 
#           residual sum of squares (default)
#       - "%IncMSE" for importance measure obtained by permutation of 
#           OOB data
#   -- seed: random number generator seed for replication of analyses
#       (default NULL means the seed is not reset)
#   -- trace: index of currently computed gene is reported (default TRUE)
#   -- All additional parameters are passed to the randomForest function
#       (see randomForest manual for more info)

# Returns:
# weighted adjacency matrix of inferred network.
# element w_ij (row i, column j) gives the importance of the link
# from regulatory gene i to target gene i

.get.weight.matrix <- function(expr.matrix, K="sqrt", nb.trees=1000, 
    input.idx=NULL, importance.measure="IncNodePurity", seed=NULL, 
    trace=TRUE, ...)
{
    # set random number generator seed if seed is given
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # to be nice, report when parameter importance.measure is not 
    # correctly spelled
    if (importance.measure != "IncNodePurity" && 
        importance.measure != "%IncMSE") {
            stop("Parameter importance.measure must be 
            \"IncNodePurity\" or \"%IncMSE\"")
    }
    # transpose expression matrix to (samples x genes)
    expr.matrix <- t(expr.matrix)
    # normalize expression matrix
    expr.matrix <- apply(expr.matrix, 2, 
        function(x) { (x - mean(x)) / sd(x) } )
    # setup weight matrix
    num.samples <- dim(expr.matrix)[1]
    num.genes <- dim(expr.matrix)[2]
    gene.names <- colnames(expr.matrix)
    weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
    rownames(weight.matrix) <- gene.names
    colnames(weight.matrix) <- gene.names
    # get number of input genes, names of input genes
    if (is.null(input.idx)) {
        num.input.genes <- num.genes
        input.gene.names <- gene.names
    } else {
        num.input.genes <- length(input.idx)
        # input gene indices given as integers
        if (is.numeric(input.idx)) {
            input.gene.names <- gene.names[input.idx]
            # input gene indices given as names
        } else {
            input.gene.names <- input.idx
            # for security, abort if some input gene name is not in gene names
            missing.gene.names <- setdiff(input.gene.names, gene.names)
            if (length(missing.gene.names) != 0) {
                for (missing.gene.name in missing.gene.names) {
                    message(paste("Gene ", missing.gene.name,
                        " was not in the expression matrix\n", sep=""))
                }
                stop("Aborting computation")
            }
        }
    }
    # set mtry
    if (class(K) == "numeric") {
        mtry <- K
    } else if (K == "sqrt") {
        mtry <- round(sqrt(num.input.genes))
    } else if (K == "all") {
        mtry <- num.input.genes-1
    } else {
        stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
    }
    if (trace) {
        message(paste("Starting RF computations with ", nb.trees,
            " trees/target gene,\nand ", mtry,
            " candidate input genes/tree node\n",
            sep=""))
        flush.console()
    }
    # compute importances for every target gene
    for (target.gene.idx in seq(from=1, to=num.genes)) {
        if (trace) {
            message(paste("Computing gene ", target.gene.idx, "/", num.genes,
                "\n", sep=""))
            flush.console()
        }
        target.gene.name <- gene.names[target.gene.idx]
        # remove target gene from input genes
        these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
        x <- expr.matrix[,these.input.gene.names]
        y <- expr.matrix[,target.gene.name]
        rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, 
            importance=TRUE, ...)
        im <- importance(rf)[,importance.measure]
        im.names <- names(im)
        weight.matrix[im.names, target.gene.name] <- im
    }
    return(weight.matrix / num.samples)
}

# utility function to convert linear index to (row,col) index for matrix
.lin.to.square <-  function(i, nrow) {
    col <- ((i - 1) %/% nrow) + 1
    row <- ((i - 1) %% nrow) + 1
    return(c(row, col))
}