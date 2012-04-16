## Author: Dan Mcglinn
## Description: This script contains functions for the computing probabilities 
## related to the following models: HEAP, bisection model, and METE. These
## probabilities primarily are related to the spatial abundance distribution, but
## not exclusively.  
## $Id$ 

piBin = function(n,A,no,Ao){
  p = A/Ao
  (factorial(no)/(factorial(n)*factorial(no-n))) * p^n * (1 - p)^(no - n)
}

g = function(no,Mo){
  ## Harte book Eq. 4.10 pg. 93
  ## returns the number of ways no indistinguishable individuals can be 
  ## arranged into Mo cells
  ## The algebraic form of the function is: 
  ## factorial(Mo + no - 1) / (factorial(no) * factorial(Mo - 1))
  ## But the following function makes numerical shortcuts so that fewer
  ## terms must be multiplied together
  num = 1
  den = 1
  termDiff = no - Mo
  nDrop = ifelse(termDiff >= 0,termDiff + 1,0)
  if(nDrop > 0){
    for(i in 1:(no - nDrop)){
      num = num * (no + Mo - i)
    }  
    den = factorial(no - nDrop)
  }
  else{
    i = 0 
    while(i < no){
      num = num * (no + Mo - 1 - i)
      i = i + 1
    }
    den = factorial(no)
  }
  return(num/den)
}

piLap = function(n,A,no,Ao){
  ## Generalized Laplace
  ## Harte book Eq. 4.11 pg. 93
  ## returns the probability that n individuals are located in a randomly
  ## chosen cell that is Mo times smaller than the area it is embedded within.
  ## no is the total number of individuals in the larger area.
 Mo = Ao / A
 sapply(n,function(x)g(no - x,Mo - 1) / g(no,Mo))
}

loadBisect = function(...){
  if(!is.loaded("piBisect")){
    OS = Sys.info()['sysname']
    if(OS == 'Linux')
      dyn.load('bisect.so')
    else
      dyn.load('bisect.dll')
  } 
}

getF = function(a,n){
  ## Conlisk et al. (2007)
  ## Eq. 7
  ## this is a computationaly efficient way to compute 
  ## gamma(a+n) / (gamma(a)*gamma(n+1))
  if(n == 0)
    out = 1
  else
    out = prod(sapply(1:n,function(i) (a+i-1)/i ))
  return(out)
}  

piSingle = function(n,A,no,Ao,psi,c=2){
  ## Single division model
  ## Conclisk et al. (2007)
  ## Theorem 1.3
  if(psi <= 0 | psi >= 1){
    out = 0
  }
  else{
    a = (1 - psi) / psi
    out = (getF(a,n) * getF((c-1)*a,no-n)) / getF(c*a,no) 
  }  
  return(out)
}

cdfSingle = function(A,no,Ao,psi,fast=TRUE){
  if(psi <= 0 | psi >= 1){  
    out = rep(0,no+1)
  }
  else{
    if(fast){
      loadBisect()
      out = .C("cdfSingle",A=as.double(A),no=as.integer(no),Ao=as.double(Ao),
               psi=as.double(psi),cdf=as.double(rep(0,no+1)))$cdf
    }  
    else{
      out = rep(0,no+1)
      for(n in 0:no){
       if(n == 0)
           out[n+1] = piSingle(n,A,no,Ao,psi)
        else
          out[n+1] = out[n] + piSingle(n,A,no,Ao,psi)
      }
    }  
  }  
  return(out)             
}

randSingle = function(no,psi,size=1){
  rands = runif(size)
  cdf = cdfSingle(1,no,2,psi)
  xvals = sapply(rands, function(u) which(order(c(cdf,u)) == (no + 2)) - 1)
  return(xvals)
}

piBisect = function(n,A,no,Ao,psi,fast=TRUE){
  ## Bisection model
  ## Conlisk et al. (2007)
  ## Theorem 2.3
  ## psi is an aggregation parameter {0,1}
  ## Note that when psi = 0.5 that the Bisection Model = HEAP Model
  if(psi <= 0 | psi >= 1){  
    out = 0
  }
  else{
    if(fast){
      loadBisect()
      out = sapply(n,function(x){
           .C("piBisect",n=as.integer(x),A=as.double(A),no=as.integer(no),
              Ao=as.double(Ao),psi=as.double(psi),prob=as.double(0))$prob})
    }
    else{
      i = log2(Ao/A)
      if(i == 1)
        out = piSingle(n,A,no,Ao,psi)
      else{
        A = A*2 
        out = sum(sapply(n:no,function(q)
                  piBisect(q,A,no,Ao,psi,fast) * piSingle(n,A,q,Ao,psi) ))
      }
    }
  }  
  return(out)
}

loadHEAP = function(...){
  if(!is.loaded("piHEAP")){
    OS = Sys.info()['sysname']
    if(OS == 'Linux')
      dyn.load('heap.so')
    else
      dyn.load('heap.dll')
  }
}

piHEAP = function(n,A,no,Ao,fast=TRUE){
  ##HEAP model 
  ##Harte book Eq. 4.16 pg. 93
  ##The source code is in the file heap.c
  if(fast){
    loadHEAP()
    out = sapply(n,function(x){
          .C("piHEAP",n=as.integer(x),A=as.double(A),no=as.integer(no),
             Ao=as.double(Ao),prob=as.double(0))$prob})
  }
  else{
    i = log2(Ao/A)
    if(i == 1)
      out = 1/(no+1)
    else{
      A = A*2 
      out = sum(sapply(n:no,function(q) piHEAP(q,A,no,Ao,fast) / (q + 1)))
    }
  }
  return(out)
}

piHEAP2 = function(n,A,no,Ao){
  ## Eq.30 in Harte et al. (2005)
  ## this equation works only when the term 'no' is relatively small,
  ## when 'no' is a large number then rounding problems with computing the 
  ## number of combinations break the equation
  i = log2(Ao/A)
  out = sum(sapply(n:no, function(q)((-1)^(n+q))*((q+1)^-i)*
                                    exp(lchoose(no,q)+lchoose(q,n))))
  return(out)
}

piHEAP3 = function(n,A,no,Ao){
  ## Eq.30 in Harte et al. (2005)
  ## this equation works only when the term 'no' is relatively small,
  ## when 'no' is a large number then rounding problems with computing the 
  ## number of combinations break the equation
  i = log2(Ao/A)
  out = sapply(n:no, function(q) lchoose(no,q) + lchoose(q,n) - i*log((q+1)) ) 
  return(out)
}

piNegBi = function(n,A,no,Ao,k=1){
  ## Harte book Eq. 4.16 pg. 95
  nbar = no*A/Ao
  (factorial(n + k - 1) / (factorial(n)*factorial(k-1))) * (nbar / (k + nbar))^n * (k / (k+nbar))^k
}


D = function(j,L=1){
  ## Distance calculation for a golden rectangle, L x L(2^.5)
  ## From Ostling et al. (2004) pg. 130
  ## j: seperation order
  ## L: width of rectangle of area Ao  
  d = L/2^((j:1)/2)
  return(d)
}

lambda = function(i,no){ 
  ## Scaling Biodiveristy Chp. Eq. 6.4, pg.106 
  ## i: number of bisections
  if (i == 0)
    lambda = 1
  if (i != 0){
    A = 1/2^i
    lambda = 1 - piHEAP(0,A,no,1)
  }
  return(lambda)
}

chiHEAP = function(i,j,no){
  ## calculates the commonality function for a given degree of bisection (i) at 
  ## orders of seperation (j)
  ## Scaling Biodiveristy Chp. Eq. 6.10, pg.113  
  ## i: number of bisections
  ## j: order of seperation
  if(no == 1){
    out = 0
  }
  else{
    if(j == 1){
      out = (no + 1)^-1 *
            sum(sapply(1:(no-1),function(m) lambda(i-1,m) * lambda(i-1,no-m)))
    }  
    else{
      i = i-1
      j = j-1
      out = (no + 1)^-1 * sum(sapply(2:no,function(m) chiHEAP(i,j,m)))
    }
  }  
  return(out)
}

sorHEAP = function(A,no,Ao,fast=TRUE){
  ## Computes sorensen's index for a given spatial grain (A) at 
  ## all possible seperation distances 
  ## Scaling Biodiveristy Chp. Eq. 6.10, pg.113  
  ## source code in the file heap.c
  i = log2(Ao/A)
  d = D(log2(Ao/A))
  chi = matrix(NA,nrow=length(no),ncol=length(d))
  lambda = chi
  for(s in seq_along(no)){
    if(fast){
      if(!is.loaded("chiHEAP")){
        OS = Sys.info()['sysname']
        if(OS == 'Linux')
          dyn.load('heap.so')
        else
          dyn.load('heap.dll') 
      }
      chi[s,] = sapply(1:i, function(j){
                .C("chiHEAP",i=as.integer(i),j=as.integer(j),
                   no=as.integer(no[s]),prob=as.double(0))$prob })
    }                   
    else{
      chi[s,] = sapply(1:i, function(j) chiHEAP(i,j,no[s]))
    }  
    lambda[s,] = sapply(1:i, function(j) lambda(i,no[s]))
  }
  sor = apply(chi,2,sum)/apply(lambda,2,sum)
  out = data.frame(Dist = d, Sor = sor)
  return(out)
}

negllBisectVector = function(psi){
  ## single species neg log likelihood function
  ## Global variable:
  ##   dat: vector of abundance for each quadrat
  ## Local variable:
  ##   psi: aggregation parameter to be estimated
  tab = table(dat)
  n = as.numeric(names(tab))
  no = sum(dat)
  if(no == 1){
    warning('It is not appropriate to attempt to estimate psi when no = 1,
    because all psi values are equally likely')
  }  
  A = 1
  Ao = length(dat)
  freq = as.numeric(tab)
  negll = -sum(sapply(1:length(n),function(k) freq[k] * 
                      log(piBisect(n[k],A,no,Ao,psi))))
  return(negll)
}

negllBisectMatrix = function(psi){
  ## Multiple species neg log likelihood function
  ## Global variable:
  ##   dat: site x species matrix of abundance
  ## Local variable:
  ##   psi: aggregation parameter to be estimated
  noAll = colSums(dat)
  noTab = table(noAll)
  noUni = as.numeric(names(noTab))
  noFreq = as.numeric(noTab)
  negll = 0
  A = 1
  Ao = nrow(dat)
  for(i in seq_along(noUni)){
    no = noUni[i]
    sp = which(noAll == no)
    x = as.vector(dat[,sp])
    tab = table(x)
    n = as.numeric(names(tab))
    freq = as.numeric(tab)
    negll = negll - (noFreq[i] * sum(sapply(1:length(n),function(k) freq[k] *
                                            log(piBisect(n[k],A,no,Ao,psi)))))
  }  
  return(negll)
}



