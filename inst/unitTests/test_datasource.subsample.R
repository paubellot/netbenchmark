test_datasource.subsample <- function()
{
    l<-datasource.subsample(datasource=syntren300.data,
        local.noise=0,experiments = 40,datasets.num = 5)
    d<-sapply(l,dim)
    checkEquals(d[1,1],length(unique(rownames(l[[1]]))))
    checkEquals(d[1,2],length(unique(rownames(l[[2]]))))
    checkEquals(d[1,3],length(unique(rownames(l[[3]]))))
    checkEquals(d[1,4],length(unique(rownames(l[[4]]))))
    checkEquals(d[1,5],length(unique(rownames(l[[5]]))))
}