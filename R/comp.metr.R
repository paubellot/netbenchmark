#This file is an extension of a file that belongs to
#minet: Mutual Information NETworks, <http://minet.meyerp.com>
#a package that implements various algorithms for inferring mutual 
#information networks from data.

pr.plot <- function( table, device=-1, ... )
{
    table <- as.data.frame(table)
    names(table) <- sapply(names(table),tolower)
    pr <- minet::pr(table)
    if(device==-1) {
        dev.new()
        plot(pr$r,pr$p, xlab="recall",
            ylab="precision",
            main="PR-Curve",
            xlim=0:1,ylim=0:1,...)
    }else{
        dev.set(device)
        points(pr$r,pr$p, xlab="recall",
            ylab="precision", 
            main="PR-Curve",
            xlim=0:1,ylim=0:1,... )
    }
    dev.cur()
}

roc.plot <- function( table, device=-1, ... )
{
    table <- as.data.frame(table)
    names(table) <- sapply(names(table),tolower)
    roc <- minet::rates(table)
    if(device==-1) {
        dev.new()
        plot( roc$fpr,roc$tpr,
            xlab="FP rate", 
            ylab="TP rate",
            main="ROC-Curve",
            xlim=0:1,ylim=0:1,...)
    }else{
    dev.set(device)
    points( roc$fpr,roc$tpr,
            xlab="FP rate",
            ylab="TP rate",
            main="ROC-Curve",
            xlim=0:1,ylim=0:1,... )
    }
    lines( 0:1, 0:1, col="black" )
    dev.cur()
}

fscore <- function(table,beta=1)
{
    table <- as.data.frame(table)
    names(table) <- sapply(names(table),tolower)
    res <- minet::fscores(table)
    return(res)
}

auroc <- function(table,k=-1)
{ 
    table <- as.data.frame(table)
    names(table) <- sapply(names(table),tolower)
    roc <- minet::rates(table)
    if(k!=-1){
        if(k>length(roc$fpr)){
            warning("k is greater than length")
            k <- length(roc$fpr)
        }  
    }else{
        k <- length(roc$fpr)
    }
    return(pracma::trapz(roc$fpr[1:k],roc$tpr[1:k]))
}

aupr <- function(table,k=-1)
{
    table <- as.data.frame(table)
    names(table) <- sapply(names(table),tolower)
    pr <- minet::pr(table)
    if(k!=-1){
        if(k>length(pr$r)){
            warning("k is greater than length")
            k <- length(pr$r)
        }  
    }else{
        k <- length(pr$r)
    }
    return(pracma::trapz(pr$r[1:k],pr$p[1:k]))
}

.results.plot <- function(table){
    dev.new();
    par(mfrow=c(length(table),1)) 
    for(n in seq_along(table)){
        boxplot(table[[n]])
        title(main=paste("Datasource: ",as.character(n),sep=""))
    }
}

.get.pr.plot <- function(m.pr,dataset.name="",...)
{
    dev.new() 
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    nmeths <- dim(m.pr$pre)[2]
    col <- rainbow(nmeths)
    for(i in seq_len(nmeths)){
        if(i==1){
            plot(m.pr$rec[,i],m.pr$pre[,i], xlab="recall",
                ylab="precision", 
                main=paste(dataset.name,"PR-Curve"),
                xlim=0:1,ylim=0:1,col=col[i],...)
        }
        points(m.pr$rec[,i],m.pr$pre[,i],col=col[i],...)
    }
    legend("topright",inset=c(-0.3,0),colnames(m.pr$pre),col=col,
        lty=rep(1,nmeths),lwd=rep(2.5,nmeths))
}

.get.pr <- function(tp.mat,tp)
{
    s <- dim(tp.mat)
    pre <- tp.mat/matrix(rep(1:s[1],s[2]),s[1])
    rec <- tp.mat/tp
    list("pre"=pre,"rec"=rec)
}

.pgfplots.export <- function(m.pr,dataset.name="",dir,points=2000)
{
    s <- dim(m.pr$pre)
    id <- round(seq(1,s[1],length.out=points))
    names <- colnames(m.pr$pre)
    names <- sapply(names, file_path_sans_ext)
    pwd <- getwd() 
    for(i in seq_len(s[2])){
        aux <- matrix(0,points,2)
        colnames(aux) <- c("pre","rec")
        aux[,2] <- m.pr$rec[id,i]
        aux[,1] <- m.pr$pre[id,i]
        setwd(dir)
        write.table(aux,file=paste(dataset.name,"-",names[i],".txt",sep=""),
            row.names = FALSE,quote = FALSE)
    }
    setwd(pwd)
}