/*
Purpose: To compute probabilities related to Harte's HEAP model.

Arguments:
n -- the number of individuals
A -- the area
no -- the total number of individuals in the pool
Ao -- the total area of the pool
prob -- a probability value
i -- the number of bisections = log2(Ao/A)

References:
Harte, J. 2007. Towards a mechanistic basis of unified theory of 
spatial structure in ecological communities at multiple scales. in 
Scaling Biodiversity. Eds. Storch et al. 

Harte, J. 2011. Maximum Entropy and Ecology.

*/

#include <math.h> 

void piHEAP(int *n, double *A, int *no, double *Ao, double *prob){
    /*
    Computes the probability of observing n individuals in a cell 
    of area A
    */
    double piRecur(int n, double A, int no, double Ao) ;
    *prob = piRecur(*n,*A,*no,*Ao) ;
} 

double piRecur(int n, double A, int no, double Ao){
    /*
    Eq. 4.15 in Harte (2011)
    */
    int i, q ;  
    double total = 0;
    i = log2(Ao / A) ; 
    if(i == 1){
        return 1.0 / (no + 1) ; 
    }
    else{
        A = A * 2 ;
        for(q = n ; q < (no + 1) ; q++){
            total += piRecur(q,A,no,Ao) / (q + 1) ;
        }
        return total ; 
    }
}

double lambda(int i, int no){ 
    /* Eq. 6.4 in Harte (2007) */
    double piRecur(int n, double A, int no, double Ao) ; 
    double A, lamb ;
    if(no == 0){
        lamb = 1 ;
    }
    else if(i == 0) {
        lamb = 1 ; 
    }
    else{
        A = 1/pow(2,i) ;
        lamb = 1 - piRecur(0,A,no,1) ;
    }   
    return(lamb) ; 
}

void chiHEAP(int *i, int *j, int *no, double *prob){
    /*
    Computes the probability of observing a species in two quadrats of
    bisected i times and of seperation order j.
    */
    double chiRecur(int i, int j, int no) ; 
    *prob = chiRecur(*i,*j,*no) ; 
}

double chiRecur(int i, int j, int no){
    /* Eq. 6.9 and 6.10 in Harte (2007) */
    double total = 0 ; 
    int m ;
    if(no == 1){
        total += 0 ; 
    }
    else{
        if(j == 1){
            for(m = 1 ; m < no ; m++){
                total += lambda(i-1,m) * lambda(i-1,no-m) / (no + 1) ;
            }
         }  
         else{
             i-- ;
             j-- ; 
             for(m = 2 ; m < (no + 1); m++){
                 total += chiRecur(i,j,m) / (no + 1) ; 
             }  
         }
    }  
    return(total) ;
}

