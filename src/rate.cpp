#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <set>          // std::set
#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>
#include <sstream>      // std::istringstream
#include <stdio.h>
#include <string.h> 
#include <cstdlib>
#include <limits>
#include <vector>       // std::vector

using namespace Rcpp;

/* 
* Evaluation
* ------------------------------------------------------------------------ 
* Copyright (C)
* Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
*
* Pau Bellot Pujalte <pau.bellot@upc.edu>
* December 2013
* ------------------------------------------------------------------------
* This function evaluates a prediction of a network as a list of confidence 
* scores between genes and the Golden Standard:
*   Parameters: 
*       EdgeList Prediction File [ TF(string) TG(string) confidence(float)]
*       GS: Golden Standard [ TF(string) TG(string) Weight(float)]
* ------------------------------------------------------------------------
*/

void printvec(std::vector<double> vec){
    for(int i=0;i<vec.size();i++){
        std::cout<<i<<" "<< vec[i]<<std::endl;
    }
}

// [[Rcpp::export]]
NumericMatrix rate(CharacterMatrix PredEdgeList,CharacterMatrix GSEdgeList,
    int ngenes,int sym)
{
    int  npos, nl=0;
    double  TP=0,FP=0,FN,TN;
    npos = GSEdgeList.nrow();
    int N=((ngenes*ngenes)-ngenes)/sym-npos;
    std::pair <std::string,std::string> auxPair;
    std::set <std::pair <std::string,std::string> > GSset;
    std::set <std::pair <std::string,std::string> >::iterator it; 
    std::string s1, s2;
    for(int i=0;i<npos;i++){
        s1=GSEdgeList(i,0);
        s2=GSEdgeList(i,1);
        auxPair=std::make_pair(s1,s2);
        GSset.insert(auxPair); 
    }
    nl=PredEdgeList.nrow();
    std::vector<double> P(nl,0);
    int j,k;
    bool aux;
    for(int i=0;i<nl;i++){
        std::cout<<"eval "<<i<<std::endl;
        for(j=1;j+i<nl;j++){
            if(PredEdgeList(i+j,2)!=PredEdgeList(i,2)){
                break;
            }
        }
        if(j!=1){
            std::cout<<"= val from "<<i<<" to "<< i+(j-1)<<std::endl;
            std::cout<<"# = el "<< j<<std::endl;
            double mP=0;
            for(k=0;k<j;k++){
                s1=PredEdgeList(i+k,0);
                s2=PredEdgeList(i+k,1);
                auxPair=std::make_pair(s1,s2);
                it=GSset.find(auxPair);
                if(it!=GSset.end()){
                   mP=mP+1;
                }
            }
            std::cout<<"mP "<< mP<<std::endl;
            mP=mP/j;
            for(k=0;k<j;k++){
                P[i+k]=mP;
            }
            i=i+j-1;
            std::cout<<"seting i to "<< i<<std::endl;
        }else{
            s1=PredEdgeList(i,0);
            s2=PredEdgeList(i,1);
            auxPair=std::make_pair(s1,s2);
            it=GSset.find(auxPair);
            if(it!=GSset.end()){
                P[i]=1;
            }else{
                P[i]=0;
            }
        }
    }
    printvec(P); 
    NumericMatrix results(nl,4); // matrix with TP,FP,TN,FN
    for(int i=0;i<nl;i++){
        TP=TP+P[i];
        FP=FP+(1-P[i]);
        FN=npos-TP;
        TN=N-FP;
        results(i,0)=TP;
        results(i,1)=FP;
        results(i,2)=TN;
        results(i,3)=FN;
    }
    return results;
}
