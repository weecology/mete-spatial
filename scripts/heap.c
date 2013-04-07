/*
Purpose: To compute probabilities related to Harte's HEAP model and Conlisk et 
al.'s Bisection model

Arguments:
n -- the number of individuals
A -- the area of the quadrat
n0 -- the total number of individuals in the pool
A0 -- the total area of all the quadrats under consideration 
prob -- a probability value
area_ratio -- spatial scale = A0/A 

References:

Conlisk et al. 2007. A new class of models of spatial distribution.
  Ecological Monographs 77:269?284. 

Harte, J. 2007. Towards a mechanistic basis of unified theory of 
  spatial structure in ecological communities at multiple scales. in 
  Scaling Biodiversity. Eds. Storch et al. 

Harte, J. 2011. Maximum Entropy and Ecology.
*/

#include <math.h> 

/* functions for Harte's HEAP model begin here */

void heap_prob(int *n, int *A, int *n0, int *A0, double *prob){
    /*
    Computes the probability of observing n individuals in a cell 
    of area A
    */
    double heap_recur(int n, int A, int n0, int A0) ;
    *prob = heap_recur(*n, *A, *n0, *A0) ;
} 

double heap_recur(int n, int A, int n0, int A0){
    /*
    Eq. 4.15 in Harte (2011)
    */
    int area_ratio, q ;  
    double total = 0;
    area_ratio = A0 / A ; 
    if(area_ratio == 2){
        return 1.0 / (n0 + 1) ; 
    }
    else{
        A = A * 2 ;
        for(q = n ; q < (n0 + 1) ; q++){
            total += heap_recur(q, A, n0, A0) / (q + 1) ;
        }
        return total ; 
    }
}

double lambda(int i, int n0){ 
    /* Eq. 6.4 in Harte (2007) 
    i -- the number of bisections
    n0 -- the total number of individuals
    */
    double heap_recur(int n, int A, int n0, int A0) ; 
    double lamb ;
    int A0 ; 
    if(n0 == 0){
        lamb = 1 ;
    }
    else if(i == 0) {
        lamb = 1 ; 
    }
    else{
        A0 = pow(2, i) ;
        lamb = 1 - heap_recur(0, 1, n0, A0) ;
    }   
    return(lamb) ; 
}

void chi_heap(int *i, int *j, int *n0, double *prob){
    /*
    Computes the probability of observing a species in two quadrats of
    bisected i times and of seperation order j.
    */
    double chi_recur(int i, int j, int n0) ; 
    *prob = chi_recur(*i, *j, *n0) ; 
}

double chi_recur(int i, int j, int n0){
    /* Eq. 6.9 and 6.10 in Harte (2007) */
    double total = 0 ; 
    int m ;
    if(n0 == 1){
        total += 0 ; 
    }
    else{
        if(j == 1){
            for(m = 1 ; m < n0 ; m++){
                total += lambda(i - 1, m) * lambda(i - 1, n0 - m) / (n0 + 1) ;
            }
         }  
         else{
             i-- ;
             j-- ; 
             for(m = 2 ; m < (n0 + 1); m++){
                 total += chi_recur(i, j, m) / (n0 + 1) ; 
             }  
         }
    }  
    return(total) ;
}


/* functions for Conlisk et al. (2007)'s bisection model begin here */

double calc_F(double a, int n){
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

double single_prob(int n, int A, int n0, int A0, double psi){
    /* Single division model of Conlisk et al. (2007) Theorem 1.3 */
    double a ; 
    a = (1 - psi) / psi ; 
    return (calc_F(a, n) * calc_F(a, n0 - n)) / calc_F(2 * a, n0) ;
}

void single_cdf(int *A, int *n0, int *A0, double *psi, double *cdf){
    /* cumulative density function for the single single model of 
    Conlisk et al. (2007) */
    int n ; 
    for(n = 0 ; n < (*n0 + 1) ; n++){
        if(n == 0){
            cdf[n] = single_prob(n, *A, *n0, *A0, *psi) ;
        }
        else{
            cdf[n] = cdf[n-1] + single_prob(n, *A, *n0, *A0, *psi) ;
        }
    }
}

void bisect_prob(int *n, int *A, int *n0, int *A0, double *psi, double *prob){
    /*
    Computes the probability of observing n individuals in a cell 
    of area A
    */
    double bisect_recur(int n, int A, int n0, int A0, double psi) ;
    *prob = bisect_recur(*n, *A, *n0, *A0, *psi) ;
} 

double bisect_recur(int n, int A, int n0, int A0, double psi){
    /*
    Theorem 2.3 in Conlisk et al. (2007)
    */
    int area_ratio, q ;  
    double total = 0 ;
    double a ; 
    area_ratio = A0 / A ; 
    a = (1 - psi) / psi ; 
    if(area_ratio == 2){
        return single_prob(n, A, n0, A0, psi) ; 
    }
    else{
        A = A * 2 ;
        for(q = n ; q < (n0 + 1) ; q++){
            total += bisect_recur(q, A, n0, A0, psi) * single_prob(n, A, q, A0, psi) ;
        }
        return total ; 
    }
}