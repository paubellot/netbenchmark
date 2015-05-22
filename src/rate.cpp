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
    for(int i=0;i<nl;i++){
        for(j=1;j+i<nl;j++){
            if(PredEdgeList(i+j,2)!=PredEdgeList(i,2)){
                break;
            }
        }
        if(j!=1){
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
            mP=mP/j;
            for(k=0;k<j;k++){
                P[i+k]=mP;
            }
            i=i+j-1;
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
