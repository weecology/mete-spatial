
setwd('/home/danmcglinn/maxent/spat')
source('spat_analytical_prob_funcs.R')

ll.Bisect = function(n,A,no,Ao,psi){
  sum(log(sapply(n,function(x) piBisect(x,A,no,Ao,psi))))
}

ll.Bisect2 = function(n,A,no,Ao,psi){
  tab = table(n)
  n = as.numeric(names(tab))
  freq = as.numeric(tab)
  sum(sapply(1:length(n),function(i) freq[i] * log(piBisect(n[i],A,no,Ao,psi))))
}

neg.ll.Bisect = function(psi){
  -sum(log(sapply(dat,function(x) piBisect(x,A,no,Ao,psi))))
}

neg.ll.Bisect2 = function(psi){
  tab = table(dat)
  n = as.numeric(names(tab))
  freq = as.numeric(tab)
  -sum(sapply(1:length(n),function(i) freq[i] * log(piBisect(n[i],A,no,Ao,psi))))
}

bci = read.csv('./data/bci_comms.csv')
unique(bci$grain)
bci = bci[bci$grain==3906,]
bci = bci[,-(1:3)]
bci = bci[,order(apply(bci,2,sum),decreasing=TRUE)]
hist(bci[,40])
dat = bci[,40]
#dat = sample(bci[,10],2^5)
dat = dat[sample(length(dat),8)]
hist(dat)
A  = 1
Ao = length(dat)
no = sum(dat)
ll.Bisect2(dat,1,no,Ao,.5)
ll.Bisect2(dat,1,no,Ao,.001)

##nlm
start = proc.time()[3]
est = nlm(neg.ll.Bisect2,.5)
end = proc.time()[3]
nlmtime = end - start
##
fp = function(x) {print(x) ; neg.ll.Bisect2(x)}
start = proc.time()[3]
est2 = optimize(fp,c(0,1),tol=1e-3)
end = proc.time()[3]
optizetime = end - start
##
system.time(optimize(neg.ll.Bisect,c(0,1),tol=1e-3))

plot(table(dat),main=paste('psi est',est$est,sep=' '))
lines(0:6,sapply(0:6,function(n) piBisect(n,1,no,Ao,est$est))*no,type='o',col='red')
lines(0:6,sapply(0:6,function(n) piBisect(n,1,no,Ao,.5))*no,type='o',col='blue')
lines(0:6,sapply(0:6,function(n) piBisect(n,1,no,Ao,.75))*no,type='o',col='green')


dat = rpois(2^3,5)
no = sum(dat)
Ao = length(dat)
ll.Bisect(dat,1,no,Ao,.5)
ll.Bisect(dat,1,no,Ao,.001)

psi = seq(.001,.99,length.out=10)
plot(psi,sapply(psi,function(j) ll.Bisect(dat,1,no,Ao,j)),type='o')

?nlm
A = 1
nlm(neg.ll.Bisect,.25)

dat = randSingle(10,.75,2^4)
no = sum(dat)
A = 1
Ao = length(dat)
ll.Bisect(dat,1,no,Ao,.5)
ll.Bisect(dat,1,no,Ao,.001)
ll.Bisect(dat,1,no,Ao,.75)

nlm(neg.ll.Bisect,.25)



