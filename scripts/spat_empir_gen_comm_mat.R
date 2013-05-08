## Purpose: to create site x species matrices for census 7
## which was surveyed in 2010 for each spatial scale of interest
## Metadata: bci50ha.doc 

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

## read in data from the 2010 census (i.e. census 7)

dat = read.csv('./data/filtered_data/bci_census7_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$Latin)))
dat$spnum = match(dat$Latin,uniSpeciesNames)
S = max(dat$spnum) 

i_bisections = c(13, 11, 9, 7, 5, 3)
n_quadrats = 2^i_bisections
domain = c(0,1000,0,500) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$gx, dat$gy), n_quadrats, domain)

write.csv(comms, file='./data/bci_comms.csv', row.names=FALSE)

##------------------------------------------------------------------------------

## Purpose: to create site x species matrices for census 3 of cocoli
## which was surveyed in 1998 for each spatial scale of interest
## Metadata: cocoli_files.doc 
## There is a visual diagram of the Sherman plot in the following publication,
## which is fairly close to the Cocoli plot except for two 40 m strips 
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

## read in data from the 1998 census (i.e. census 3)

dat = read.csv('./data/filtered_data/cocoli_census3_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$spcode)))
dat$spnum = match(dat$spcode,uniSpeciesNames)
S = max(dat$spnum) 

plot(dat$x,dat$y,pch='.')

## split plot into two 200 x 100 m rectangles
dat1 = dat[dat$y>=100,]
dat2 = dat[dat$y<100,]

plot(dat1$x,dat1$y,col='red',pch='.',ylim=range(dat$y),xlim=range(dat$x))
points(dat2$x,dat2$y,col='blue',pch='.')

i_bisections = c(13, 11, 9, 7, 5)
n_quadrats = 2^i_bisections
domain1 = c(0,100,100,300) # spatial domain in meters defined here
domain2 = c(0,200,0,100) 

## generate a site x species matrix for each spatial scale

comms1 = make_comm_matrix(dat$spnum, S, cbind(dat$x, dat$y), n_quadrats,
                          domain1)
comms2 = make_comm_matrix(dat$spnum, S, cbind(dat$x, dat$y), n_quadrats,
                          domain2)

## output results
write.csv(comms1,file='./data/cocoli1_comms.csv',row.names=F)
write.csv(comms2,file='./data/cocoli2_comms.csv',row.names=F)

##------------------------------------------------------------------------------

## Purpose: to create site x species matrices for census 3 of sherman
## which was surveyed in 1999 for each spatial scale of interest
## Metadata: sherman_files.doc 
## There is a visual diagram of the Sherman plot in the following publication,
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

## read in data from the 1999 census (i.e. census 3)

dat = read.csv('./data/filtered_data/sherman_census3_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$spcode)))
dat$spnum = match(dat$spcode,uniSpeciesNames)
S = max(dat$spnum) 

plot(dat$x,dat$y,pch='.')
abline(v=140,col='red')
abline(v=240,col='red')
abline(h=40,col='red')
abline(h=140,col='red')
abline(h=440,col='red')

## split plot into three quadrats (two 200 x 100 m rectangles and one 140 x 140m)
dat1 = dat[dat$x >= 140 & dat$y >= 240, ]  ## one of the 200 x 100 m rectangles
dat2 = dat[dat$x >= 140 & dat$y < 240, ] ## the other 200 x 100 m rectangle
dat3 = dat[dat$x < 140, ]  ## the 140 x 140 square

plot(dat1$x,dat1$y,col='red',pch='.',ylim=range(dat$y),xlim=range(dat$x))
points(dat2$x,dat2$y,col='blue',pch='.')
points(dat3$x,dat3$y,col='green3',pch='.')

i_bisections = c(13, 11, 9, 7, 5)
n_quadrats = 2^i_bisections
domain1 = c(140,240,240,440) # spatial domain in meters defined here
domain2 = c(140,240,40,240)

## generate a site x species matrix for each spatial scale
comms1 = make_comm_matrix(dat1$spnum, S, cbind(dat1$x, dat1$y), n_quadrats,
                          domain1)
comms2 = make_comm_matrix(dat2$spnum, S, cbind(dat2$x, dat2$y), n_quadrats,
                          domain2)

## now for dat3 the square plot
i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain3 = c(0,140,0,140)

comms3 = make_comm_matrix(dat3$spnum, S, cbind(dat3$x, dat3$y), n_quadrats,
                          domain3)

## output results
write.csv(comms1,file='./data/sherman1_comms.csv',row.names=F)
write.csv(comms2,file='./data/sherman2_comms.csv',row.names=F)
write.csv(comms3,file='./data/sherman3_comms.csv',row.names=F)

##------------------------------------------------------------------------------

## Purpose: to create site x species matrices for the crosstimbers dataset
## which was surveyed in 1998 on a 64 m^2 plot.
## Metadata: CrosstimbersMasterDataSheet1-20-12.xlsx

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/cross1998_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$sp)))
dat$spnum = match(dat$sp, uniSpeciesNames)
S = max(dat$spnum) 

## for square quadrats the lengths are in meters defined here
i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0, 200, 0, 200) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$x, dat$y), n_quadrats, domain)

write.csv(comms, file='./data/cross_comms.csv', row.names=FALSE)

##------------------------------------------------------------------------------

## Purpose: to create site x species matrices for the ferp woody plants dataset
## Metadata: http://ferp.ucsc.edu/ferp/about/using/

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/ferp_2007_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$Latin)))
dat$spnum = match(dat$Latin, uniSpeciesNames)
S = max(dat$spnum) 

range(dat$gx) ## max 200
range(dat$gy) ## max 300
## change into a single 150 x 300 quadrat
trim = (200 - 150) / 2

i_bisections = c(13, 11, 9, 7, 5)
n_quadrats = 2^i_bisections
domain = c(trim, 200 - trim, 0, 300) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$gx, dat$gy), n_quadrats, domain)

write.csv(comms, file='./data/ferp_comms.csv', row.names=FALSE)

##------------------------------------------------------------------------------

## Purpose: to create site x species matrices for the oosting tree dataset
## Metadata: http://esapubs.org/archive/ecol/E088/162/metadata.htm

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/oosting_trees_1990_filtered.csv')

S = length(unique(dat$CODE))
dat$spnum = as.integer(dat$CODE)

i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0, 256.00001, 0, 256.00001)

comms = make_comm_matrix(dat$spnum, S, cbind(dat$TX, dat$TY), n_quadrats, domain)

## output results
write.csv(comms,file='./data/oosting_comms.csv',row.names=F)

##------------------------------------------------------------------------------
## Purpose: to create a site by species matrix for the the Luquillo forest
## dataset for analysis of spatial biodiversity patterns, this survey is
## from 2006/2007
## Metadata: http://luq.lternet.edu/data/luqmetadata119
## see the following page for the status codes defined:
## http://luq.lternet.edu/data/variable/codes

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/luquillo_census4_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$SPCODE)))
dat$spnum = match(dat$SPCODE, uniSpeciesNames)
S = max(dat$spnum) 

range(dat$GX) ## max 320
range(dat$GY) ## max 500
## change into a single 250 x 500 quadrat
trim = (320 - 250) / 2

i_bisections = c(13, 11, 9, 7, 5, 3) 
n_quadrats = 2^i_bisections
domain = c(trim, 320 - trim, 0, 500) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$GX, dat$GY), n_quadrats, domain)

write.csv(comms, file='./data/luquillo_comms.csv', row.names=FALSE)

##------------------------------------------------------------------------------
## Purpose: to create site x species matrices for the serpentine dataset
## which was surveyed in 1998 on a 64 m^2 plot.
## Metadata: http://socrates.berkeley.edu/~hartelab/MaxEnt.html

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/raw_data/serpentine_data.csv')

i_bisections = c(8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(1, 17, 1, 17) # defined in number of quadrats here

# prepare data for the make_comm_matrix function
abu = as.vector(as.matrix(dat))
S = ncol(dat)
spnum = rep(1:S, each=256)
x = rep(rep(1:16, each=16), S)
y = rep(rep(1:16, times=16), S)

input_dat = data.frame(spnum, x, y, abu)[abu > 0,]

comms = make_comm_matrix(input_dat$spnum, S, input_dat[,2:3], n_quadrats, domain,
                         input_dat$abu)

## fix grain names, so that they are in units of m2
comms[ , 1] = as.numeric(comms[ , 1]) * 0.25

write.csv(comms, file='./data/serp_comms.csv',row.names=FALSE)

##------------------------------------------------------------------------------
## Purpose: to create site x species matrics for the NC plots f
## Metadata: is in the raw data files and in ~/datasets/NC_LTREB/Maps/Readme_maps.txt 

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat_names = paste('m', c('04','07','12','13',91:94), sep='')
plc_names = c('graveyard','landsend','rocky','bormann',
              'woodbridge','baldmnt','bryan','bigoak')
yrs = c(92, 93, 90, 93, 91, 91, 91, 93)

file_names = paste(dat_names, '_', plc_names, '_19', yrs, '_filtered.csv', sep='') 

dat = vector('list', length(file_names))
for (i in seq_along(file_names))
  dat[[i]] = read.csv(file.path('./data/filtered_data', file_names[i]))
names(dat) = plc_names

comms = vector('list', length(file_names))
names(comms) = plc_names

######### graveyard
i = 1 
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 100
range(tmp$Y) ## max 100

i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0, 100, 0, 100) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

######### landsend
i = 2 ## 
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 130
range(tmp$Y) ## max 65

i_bisections = c(13, 11, 9, 7, 5)
n_quadrats = 2^i_bisections
domain = c(0, 130, 0, 65) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

######### rocky
i = 3
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 120
range(tmp$Y) ## max 120

i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0, 120, 0, 120) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

######### bormann
i = 4
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 140
range(tmp$Y) ## max 140

i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0, 140, 0, 140) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

######### wood bridge
i = 5
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 71
range(tmp$Y) ## max 71

i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0, 71, 0, 71) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

######### baldmnt
i = 6
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 50
range(tmp$Y) ## max 100

i_bisections = c(11, 9, 7, 5)
n_quadrats = 2^i_bisections
domain = c(0, 50, 0, 100) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

######### bryan
i = 7
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 185
range(tmp$Y) ## max 92.5

i_bisections = c(13, 11, 9, 7, 5)
n_quadrats = 2^i_bisections
domain = c(0, 185, 0, 92.5) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

######### bigoak
i = 8
tmp = dat[[i]]
plot(tmp$X, tmp$Y) 

uniSpeciesNames = as.character(sort(unique(tmp$SPEC)))
tmp$spnum = match(tmp$SPEC, uniSpeciesNames)
S = max(tmp$spnum) 

range(tmp$X) ## max 200
range(tmp$Y) ## max 100

i_bisections = c(13, 11, 9, 7, 5)
n_quadrats = 2^i_bisections
domain = c(0, 200, 0, 100) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale
comms[[i]] = make_comm_matrix(tmp$spnum, S, cbind(tmp$X, tmp$Y), n_quadrats, domain)

#########
## output the comms object to seperate files.
for (i in seq_along(comms)) {
  write.csv(comms[[i]], file=paste('./data/', plc_names[i], '_comms.csv', sep=''),
            row.names=F)
}


