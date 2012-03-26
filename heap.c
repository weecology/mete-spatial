/*
Purpose: To compute the HEAP model, Eq. 4.15 in Harte's book
Arguments:
n -- the number of individuals
A -- the area
no -- the total number of individuals in the pool
Ao -- the total area of the pool
prob -- the probability of observing n individuals (what we want back)

*/

#include <R.h>
#include <math.h> 

void HEAP(int *n, double *A, int *no, double *Ao, double *prob){
    double rHEAP(int n, double A, int no, double Ao) ;
    *prob = rHEAP(*n,*A,*no,*Ao) ;
} 

double rHEAP(int n, double A, int no, double Ao){
    int i, q ;  
    double total = 0;
    i = log2(Ao / A) ; 
    if(i == 1){
        return 1.0 / (no + 1) ; 
    }
    else{
        A = A * 2 ;
        for (q = n ; q < (no + 1) ; q++){
            total += rHEAP(q,A,no,Ao) / (q + 1) ;
        }
        return total ; 
    }
}