
piBin = function(n,A,no,Ao){
  p = A/Ao
  (factorial(no)/(factorial(n)*factorial(no-n))) * p^n * (1 - p)^(no - n)
}

g = function(no,Mo){
  ## Harte book Eq. 4.10 pg. 93
  ## returns the number of ways no indistinguishable individuals can be 
  ## arranged into Mo cells
  factorial(Mo + no - 1) / (factorial(no) * factorial(Mo - 1))
}

piLap = function(n,A,no,Ao){
  ## Harte book Eq. 4.11 pg. 93
  ## returns the probability that n individuals are located in a randomly
  ## chosen cell that is Mo times smaller than the area it is embedded within.
  ## no is the total number of individuals in the larger area.
  Mo = Ao / A
  g(no - n, Mo - 1) / g(no, Mo)
}

piHEAP = function(n,A,no,Ao){
  ## Harte book Eq. 4.15 pg. 94
  i = log2(Ao/A)
  if(length(n) > 1){
    sapply(1:length(n),function(j){
      sum(sapply(n[j]:no, function(q) piLap(q,Ao/2^(i-1),no,Ao) / (q + 1)))
    })
  }  
  else{
    sum(sapply(n:no,function(q) piLap(q,Ao/2^(i-1),no,Ao) / (q + 1)))
  }
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

piMETEiter = function(n,A,no,Ao){
  ## iterative METE approach, downscaling only
  ## special case when A = Ao / 2 applies
  ## Harge Book Eq. 7.51 pg. 159
  i = log2(Ao/A)
  if(length(n) > 1){
    sapply(1:length(n),function(j){
      sum(sapply(n[j]:no, function(q) piMETE(q,Ao/2^(i-1),no,Ao) / (q + 1)))
    })
  }  
  else{
    sum(sapply(n:no,function(q) piMETE(q,Ao/2^(i-1),no,Ao) / (q + 1)))
  }
}

## generate some comparitive frequency distributions
## example from pg 96 fig 4.1
setwd('c:/users/white lab/documents/lab data/maxent/spat')
pdf('analytical_prob_funcs.pdf',width=14,height=7)
par(mfrow=c(1,2))
no = 100
n = 0:5
Ao = 64
A = 1
out = matrix(NA,nrow=6,ncol=length(n))
out[1,] = piBin(n,A,no,Ao)
out[2,] = piLap(n,A,no,Ao)
out[3,] = piHEAP(n,A,no,Ao)
out[4,]= piMETE(n,A,no,Ao)
out[5,]= piMETEiter(n,A,no,Ao)
out[6,] = piNegBi(n,A,no,Ao)


plot(n,out[1,],ylim=range(out,na.rm=TRUE),type='n',ylab='Probabiliy',
     main='No = 100, A = Ao/64')
for(i in 1:4)
  lines(n,out[i,],col=i,type='l',lwd=2)
for(i in 5:6)
 points(n,out[i,],col=i,pch=19,cex=1.5)
legend('topright',c('bin','lap','heap','mete','meteiter','negbin(k=1)'),col=1:6,
       lwd=c(rep(3,4),NA,NA),lty=c(rep(1,4),NA,NA),pch=c(rep(NA,4),19,19),cex=2,bty='n')
## lap, mete and negbi are equivalent
## heap and meteiter are equivalent

no = 5
n = 0:5
Ao = 4
A = 1

out = matrix(NA,nrow=6,ncol=length(n))
out[1,] = piBin(n,A,no,Ao)
out[2,] = piLap(n,A,no,Ao)
out[3,] = piHEAP(n,A,no,Ao)
out[4,]= piMETE(n,A,no,Ao)
out[5,]= piMETEiter(n,A,no,Ao)
out[6,] = piNegBi(n,A,no,Ao)

plot(n,out[1,],ylim=range(out,na.rm=TRUE),type='n',ylab='Probability',
     main='No = 5, A = Ao/4')
for(i in 1:4)
  lines(n,out[i,],col=i,type='l',lwd=2)
for(i in 5:6)
 points(n,out[i,],col=i,pch=19,cex=1.5)
legend('topright',c('bin','lap','heap','mete','meteiter','negbin(k=1)'),col=1:6,
       lwd=c(rep(3,4),NA,NA),lty=c(rep(1,4),NA,NA),pch=c(rep(NA,4),19,19),cex=2,bty='n')
## mete and negbi are equivalent
## heap and meteiter are equivalent
dev.off()
