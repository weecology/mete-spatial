
setwd('/home/danmcglinn/maxent/spat')
source('spat_analytical_prob_funcs.R')

ll.Bisec = function(n,A,no,Ao,psi){
  sum(log(sapply(n,function(x) piBisec(x,A,no,Ao,psi))))
}

neg.ll.Bisec = function(psi){
  -sum(log(sapply(dat,function(x) piBisec(x,A,no,Ao,psi))))
}

dat = rpois(2^3,5)
no = sum(dat)
Ao = length(dat)
ll.Bisec(dat,1,no,Ao,.5)
ll.Bisec(dat,1,no,Ao,.001)

psi = seq(.001,.99,length.out=10)
plot(psi,sapply(psi,function(j) ll.Bisec(dat,1,no,Ao,j)),type='o')

?nlm
A = 1
nlm(neg.ll.Bisec,.25)

