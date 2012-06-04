## Author: Dan McGlinn
## Purpose: to create site x species matrices for the serpentine dataset
## which was surveyed in 1998 on a 64 m^2 plot.
## Metadata: http://socrates.berkeley.edu/~hartelab/MaxEnt.html

setwd('/home/danmcglinn/maxent/spat')

source('spat_sim_vario_func.R')

dat = read.csv('./data/serpentine_data.csv',header=T)


i_bisections = c(7:1)
n_quadrats = 2^i_bisections
domain = c(1, 17, 1, 17) # defined in number of quadrats here

# prepare data for the make_comm_matrix function
abu = as.vector(as.matrix(dat))
S = ncol(dat)
spnum = rep(1:S, each=256)
x = rep(rep(1:16, each=16), S)
y = rep(rep(1:16, times=16), S)

input_dat = data.frame(spnum, x, y, abu)[abu > 0,]

comms1 = make_comm_matrix(input_dat$spnum, S, input_dat[,2:3], n_quadrats, domain,
                         input_dat$abu)

abu = as.vector(as.matrix(dat))
S = ncol(dat)
spnum = rep(1:S, each=256)
y = rep(rep(1:16, each=16), S)
x = rep(rep(1:16, times=16), S)

input_dat = data.frame(spnum, x, y, abu)[abu > 0,]

comms2 = make_comm_matrix(input_dat$spnum, S, input_dat[,2:3], n_quadrats, domain,
                         input_dat$abu)

write.csv(comms, file='./data/serp_comms.csv',row.names=FALSE)


      
  