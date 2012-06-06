
setwd('/home/danmcglinn/maxent/spat')
source('spat_analytical_prob_funcs.R')


piHEAP(5,1,5,4)
(1/(5+1))^2

piHEAP(4,1,5,4)
(1/(5+1))*(1/(5+1))+(1/(5+1)*(1/(4+1)))

piHEAP(3,1,5,4)
(1/(5+1))*(1/(5+1))+(1/(5+1)*(1/(4+1))) + (1/(5+1)*(1/(3+1)))




logFactorialTab = read.delim('logfactorial_lookup_table.txt',header=F)
factHash = hash(keys = logFactorialTab[,1], values = logFactorialTab[,2])




serp = read.csv('./data/serp_comms.csv')
dat = as.matrix(serp[serp$grain==16,-(1:3)])
colSums(dat)

no = 2^(1:9) 
notimes = NA
for(i in seq_along(no)){
  notimes[i] = system.time(piBisect2(0,1,no[i],256,.25))[3]
  print(no[i])
}

no = no[1:8]
plot(log2(no),log2(notimes))
mod = lm(log2(notimes) ~ log2(no))
abline(mod)
round(2^(predict(mod,newdata = data.frame(no = colSums(dat))))/(60),3)
round(2^(predict(mod,newdata = data.frame(no = colSums(dat))))/(60^2),3)
round(2^(predict(mod,newdata = data.frame(no = colSums(dat))))/(60^2 * 24),3)



dat = dat[,10:11]

fp = function(x) {print(x) ; negllBisectMatrix(x)}

start = proc.time()[3]
est = optimize(fp,c(0,1),tol=1e-3)
end = proc.time()[3]

write.csv(est,file='./mle/serp_comm_est.csv')
paste('Optimization took',round(end - start),'seconds')

