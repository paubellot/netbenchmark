#include <Rcpp.h>
#include <math.h>       /* sqrt */
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix zsc(NumericMatrix x) 
{
    int nrow = x.nrow(), ncol = x.ncol(), pmin, i, j;
    double min, tmp;
    NumericVector mean(ncol),var(ncol),t(ncol);
    NumericMatrix z(nrow,ncol),Z(ncol,ncol);
    for(i = 0; i < nrow; i++){
        for(j = 0; j < ncol; j++){
            if(i<ncol){
                Z(i,j)=0;
            } 
            z(i,j)=0;
        }
    }
    for(j = 0; j < ncol; j++){
        t(j)=0;
        mean(j)=0;
        for(i = 0; i < nrow; i++){
            mean[j]+=x(i,j);
        } 
        mean[j]/=nrow;
    }
    for(j = 0; j < ncol; j++){
        var[j]=0;
        for(i = 0; i < nrow; i++){
            tmp = (x(i,j)-mean[j]); 
            var[j] += tmp*tmp;
        }
        var[j]/=(nrow-1);
        var[j]=sqrt(var[j]);
    }
    for(j = 0; j < ncol; j++){
        for(i = 0; i < nrow; i++){
            z(i,j)=(x(i,j) - mean[j])/var[j];
        }
    }
    for (i = 0; i < nrow; i++) {
        min = x(i,0);
        pmin = 0;
        for (j = 1; j < ncol; j++) {
            if(x(i, j)<min){
                min=x(i, j);
                pmin=j;
            }
        }
        t(pmin)++;
        for (j = 0; j < ncol; j++) {
            Z(pmin,j)=Z(pmin,j)+z(i,j);
        }
    }
    for(i = 0; i< ncol; i++){
        if(t(i)!=0){
            for (j = 0; j < ncol; j++) {
                Z(i,j)= Z(i,j)/t(i);
            }
        }else{
            for (j = 0; j < ncol; j++) {
                Z(i,j)= 0;
            }
        }
    }
    return Z;
}
