## $Id$
## Conlisk et al. (2007) fig3
## Question: does the multivariate probability of a particular landscape
## depend upon if one first bisects vertically or horizontally, this code
## demonstrates that this is the case at least according to conlisk's example
setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

dat = c(0,1,0,0,3,0,11,0,2,4,6,3,5,10,4,0)

dat = round(runif(64) * 10)

i_bisections = 0:6
n_quadrats = 2^i_bisections
xdim = ydim = sqrt(length(dat))

domain = c(1, xdim + 1, 1, ydim + 1) # defined in number of quadrats here

# prepare data for the make_comm_matrix function
abu = as.vector(as.matrix(dat))
S = 1
spnum = rep(1:S, each=length(dat))
x = rep(rep(1:xdim, each=xdim), S)
y = rep(rep(1:ydim, times=ydim), S)

input_dat = data.frame(spnum, x, y, abu)[abu > 0, ]

comms1 = make_comm_matrix(input_dat$spnum, S, input_dat[ ,2:3], n_quadrats, domain,
                          input_dat$abu)

abu = as.vector(as.matrix(dat))
S = 1
spnum = rep(1:S, each=length(dat))
y = rep(rep(1:xdim, each=xdim), S)
x = rep(rep(1:ydim, times=ydim), S)

input_dat = data.frame(spnum, x, y, abu)[abu > 0, ]
comms2 = make_comm_matrix(input_dat$spnum, S, input_dat[ ,2:3], n_quadrats,
                          domain, input_dat$abu)

abu1 = as.numeric(comms1[ , 4])
abu2 = as.numeric(comms2[ , 4])
prod(c(1/50, 1/16, 1/35, 1/5, 1/12, 1/22, 1/14, 1/2, 1/4, 1, 1/12, 1/7, 1/16, 1/10, 1/5))
prod(ifelse(abu1>0, 1/(abu1 + 1), 1))
prod(ifelse(abu2>0, 1/(abu2 + 1), 1))

## the mulivariate probability appears to depend on if one bisects vertically or
## horizontally first... 
## possibily the solution should be to average thse two probs
0.5 * (prod(ifelse(abu1>0, 1/(abu1 + 1), 1)) + prod(ifelse(abu2>0, 1/(abu2 + 1), 1)))


