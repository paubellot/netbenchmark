RegistrerWrapper<-function(wrapper.name=NULL,all.fast=TRUE){
    tst<-paste("exists('",wrapper.name,"')",sep="")
    if(eval(parse(text=tst))){
        if(all.fast){
            aux<-c(Fast,wrapper.name)
            assign("Fast",aux, envir=.GlobalEnv)
            message(paste(wrapper.name, "added to all.fast methods"))
            
        }else{
            aux<-c(All,wrapper.name)
            assign("All",aux, envir=.GlobalEnv)
            message(paste(wrapper.name, "added to all methods"))
        }
    }else{
        stop(paste(wrapper.name,"does not exist!"))
    }
}
