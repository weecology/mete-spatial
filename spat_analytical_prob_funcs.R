
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

piHEAP = function(n,A,no,Ao){
  i = log2(Ao/A)
  if(i == 1)
    out = 1/(no+1)
  else{
    A = A*2 
    out = sum(sapply(n:no,function(q) piHEAP(q,A,no,Ao) / (q + 1)))
  }
  return(out)
}

piNegBi = function(n,A,no,Ao,k=1){
  ## Harte book Eq. 4.16 pg. 95
  nbar = no*A/Ao
  (factorial(n + k - 1) / (factorial(n)*factorial(k-1))) * (nbar / (k + nbar))^n * (k / (k+nbar))^k
}

piMETE = function(n,A,no,Ao){
  ## special case when A = Ao / 2
  ## Harge Book Eq. 7.51 pg. 159
  ## approximation of pi when A << Ao
  ## Harte book Eq. 7.53 pg. 160
  if(Ao / A == 2){
    1 / (1 + no) 
  }
  else{
    nbar = no*A/Ao
    (nbar/(1+nbar))^n / (1+nbar)
  }
}

