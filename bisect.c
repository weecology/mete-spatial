/*
Purpose: To compute probabilities related to Conlisk et al.'s 
Bisection model

Arguments:
n -- the number of individuals
A -- the area
no -- the total number of individuals in the pool
Ao -- the total area of the pool
psi -- aggregation parameter {0,1}
prob -- a probability value
i -- the number of bisections = log2(Ao/A)

References:
Conlisk et al. 2007. A new class of models of spatial distribution.
  Ecological Monographs 77:269–284. 

*/

#include <R.h>
#include <math.h> 

void piBisect(int *n, double *A, int *no, double *Ao, double *psi, double *prob){
    /*
    Computes the probability of observing n individuals in a cell 
    of area A
    */
    double piBisectRecur(int n, double A, int no, double Ao, double psi) ;
    *prob = piBisectRecur(*n,*A,*no,*Ao,*psi) ;
} 

double getF(double a, int n){
     /*
     Eq. 7 in Conlisk et al. (2007)
     */
     double out = 1 ; 
     if(n != 0){
         int i ; 
         for(i = 1 ; i < (n + 1) ; i ++){
             out *= (a + i - 1) / i ; 
         }
     }
     return out ; 
}

double piBisectRecur(int n, double A, int no, double Ao, double psi){
    /*
    Theorem 2.3 in Conlisk et al. (2007)
    */
    int i, q ;  
    double total = 0 ;
    double a ; 
    i = log2(Ao / A) ; 
    a = (1 - psi) / psi ; 
    if(i == 1){
        return (getF(a,n) * getF(a,no-n)) / getF(2*a,no) ; 
    }
    else{
        A = A * 2 ;
        for(q = n ; q < (no + 1) ; q++){
            total += piBisectRecur(q,A,no,Ao,psi) * ((getF(a,n)*getF(a,q-n)) / getF(2*a,q)) ;
        }
        return total ; 
    }
}
