## $Id: $
###################################################
##Author: Dan McGlinn
##Date: 10/25/10
##Purpose: to provide necessary support functions
##to simulate a community influenced by neutral
##and niche processes and then examine the spatial 
##pattern of associations in the resulting community
###################################################
##nuetral.vp functions
census<-function (sim, snap = length(sim$snaps), type = "census") 
{
    pop <- sim$snaps[[snap]]
    dim(pop) <- c(sim$p$S, sim$p$M, sim$p$M)
    if (type == "census") {
        output <- t(pop[, , 1])
        for (i in 2:sim$p$M) output <- rbind(output, t(pop[, 
            , i]))
    }
    else if (type == "richness") 
        output <- apply(pop, MARGIN = 2:3, FUN = function(x) sum(as.logical(x)))
    else if (type == "abundance") 
        output <- apply(pop, MARGIN = 2:3, FUN = sum)
    else stop("Invalid type")
    return(output)
}

spat.dep<-function(x,g,sp.dep=FALSE){
 ##create spatially autocorrelated comm but still species independent
 ##use Palmer & van der Maarel method for transects (may be applied to quads although strangely)
 ##use transect cords 'tcords'
 #x<-matrix(runif(length(x),min=-1,max=1),ncol=ncol(x),nrow=nrow(x))
 y<-x ##x is a site x sp matrix
 p<-apply(y,2,sum)/dim(y)[1]
 for(k in 1:g){ ##
  for(j in 1:nrow(y)){ ##sites
   if(j == 1){
    y[j,]<-y[nrow(y),]+y[1,]+y[2,]
   }
   else{
    if(j == nrow(y)){
     y[j,]<-y[j-1,]+y[j,]+y[1,]
     }
    else{
     y[j,]<-y[j-1,]+y[j,]+y[j+1,]
 }}}}
 ##convert to 0's and 1's according to prob occupancy  
 y<-sapply(1:ncol(y),function(x)ifelse(y[,x]>=quantile(y[,x],1-p[x]),1,0))
 if(sp.dep){
  ##relection algorithm
  for(i in 1:ncol(y)){
   if(runif(1)>=.5){
    y[,i]<-y[N:1,i] ##reflected
  }} 
 }
 else{
  ##eliminate inter-specific autocorr
  for(i in 1:ncol(y)){
   ind<-sample(nrow(y),1):nrow(y)
   if(length(ind)<nrow(y))
    ind<-c(ind,1:(nrow(y)-length(ind)))
   y[,i]<-y[ind,i]
 }}
 y
}


comp.rand<-function(mat){
 ##completely random null model
 rmat<-matrix(mat[sample(prod(dim(mat)))],ncol=ncol(mat),nrow=nrow(mat))
 rexpv<-sum(apply(rmat,2,var))
 rexpv
}

site.rand<-function(mat){
 ##site totals are randomized, species totals are fixed
 N<-nrow(mat)
 rmat<- apply(mat,2,function(x) x[sample(N)])
 rdists<- getCovFractions(rmat)
 rpos.cv<- sum(rdists$pos)/(N*(N-1)/2)
 rneg.cv<- sum(rdists$neg)/(N*(N-1)/2)
 rtot.cv<- rpos.cv + rneg.cv
 c(rtot.cv,rpos.cv,rneg.cv)
}


#####Part I - FUNCTIONS FOR SIMULATING ENVIORNMENTS AND COMMUNITIES#####
##1.1##
gen.env<-function(M,env='fractal',D=NA,err=.1){
 ##Purpose: to generate a 2-D landscape representation of a 
 ##single environmental gradient
 ##Arguments:
 ##'M' = number of rows or columns in the square grid
 ##'env = type of environment to generate, three options: 'fractal','grad',or 'peaks'
 ###'fractal' requires that the next argument 'D' the fractal dimension is also supplied
 ###'grad' is a strict linear gradient, 
 ###'peaks' is four similar peaks arranged in a square
 ##'D' = fractal dimension (e.g., 2 = strong gradient, 3 = no gradient)
 ##'err' = the amount of acceptable difference between the input D and the simulated D
 if(D==2)
   env='grad'
 if(D==3)
  env='random'
 if(env=='fractal'){
  require(FieldSim)
  if(is.na(D)) stop('must supply fractal dimension (D) when enviornment is a fractal')
  H = 3-D ##Hurst parameter
  res <- midpoint(H=H,nblevel=ceiling(log2(M)))##level sets the size of the grid
  res$Z<-res$Z[-nrow(res$Z),]
  res$Z<-res$Z[,-ncol(res$Z)]
  ndrop <- 2^(ceiling(log2(M))) - M
  if(ndrop > 0){
   res$Z <- res$Z[-(1:ndrop),]
   res$Z <- res$Z[,-(1:ndrop)]
  }
  z<-scale(res$Z)
  ##create new variable that will be drawn from
  z<-uni.trans(z,M)
  D.sim<- 3-quadvar(z)
  icount<-1
  while(abs(D.sim-D)>=err){
   if(is.na(D)) stop('must supply fractal dimension (D) when enviornment is a fractal')
   H = 3-D ##Hurst parameter
   res <- midpoint(H=H,nblevel=ceiling(log2(M)))##level sets the size of the grid
   res$Z<-res$Z[-nrow(res$Z),]
   res$Z<-res$Z[,-ncol(res$Z)]
   ndrop <- 2^(ceiling(log2(M))) - M
   if(ndrop > 0){
    res$Z <- res$Z[-(1:ndrop),]
    res$Z <- res$Z[,-(1:ndrop)]
   }
   z<-scale(res$Z)
   ##create new variable that will be drawn from
   z<-uni.trans(z,M)
   D.sim<- 3-quadvar(z)
 }}
 if(env == 'grad'){ ##smooth gradient
  z<-matrix(rep(seq(0,M,length.out=M),times=M),ncol=M,nrow=M)
 }
 if(env == 'peaks'){##4 peaks
  x<-rep(seq(0,M,length.out=M),each=M)
  y<-rep(seq(0,M,length.out=M),times=M)
  z<-peaks(x,y,M*c(.25,.75),M*c(.25,.75),M*.15)
  ##create new variable that will be drawn from
  z<-uni.trans(z,M)
 }
 if(env == 'random'){##totally random
  z<-matrix(runif(M^2),M,M)
  z<-uni.trans(z,M)
 }
 z
}

##1.2##
uni.trans<-function(z,len){
 ##Purpose: transforms a matrix into values that range from 0 to 1
 ##Called within the function 'gen.env'
 ##Arguments:
 ##'z' is a matrix
 ##'len' is the maximum value that the transformed matrix will take
 uni.val<-sort(unique(round(as.vector(z),3)))
 u<-seq(0,len,length.out=length(uni.val)) ##following Palmer 1992
 ##go to a cell see its z-value and replace it with a u-value
 z2<-z
 for(c in 1:ncol(z)){
  for(r in 1:nrow(z)){
   z2[r,c]<-u[uni.val==round(z[r,c],3)]
 }}
 z2
}

##1.3##
gauss.niche<-function(m,z,u,s){
 ##Purpose: to provide an exponential unimodal function 
 ##which is a model for a species response to the enviornment
 ##Called within the functions 'peaks' and 'sim.init.uni'
 ##from Palmer 1992
 ##Arguments:
 ##'m' is max perform
 ##'z' is enviornment at given coordinate
 ##'u' is env optim
 ##'s' is habitat breadth
 m*exp(-.5*(z-u)^2/s^2)
}

##1.4##
gauss.niche.diff<-function(m,diff,s){
 ##Purpose: to provide an exponential unimodal function 
 ##which is a model for a species response to the enviornment
 ##Called within the function 'sim.init.uni'
 ##from Palmer 1992
 ##m is max perform
 ##diff is the difference from current conditions and optima
 ##s is habitat breadth
 m*exp(-.5*(diff)^2/s^2)
}

##1.5##
peaks<-function(x,y,m1,m2,s){
 ##Purpose: generates gaussian peaks
 ##arranged in a square
 ##Called within the function 'gen.env'
 ##Arguments:
 ##'x'&'y' = spatial coordinates
 ##'m1'&'m2' = means along x-coord & y-coord respectively
 ##'s' = std dev of gaussian function
 z<-0
 for(i in 1:length(m1))
  for(j in 1:length(m2))
   z<- z+ gauss.niche(1,x,m1[i],s) * gauss.niche(1,y,m2[j],s) 
 z<-matrix(z,ncol=sqrt(length(x)),sqrt(length(x)))
 z
}



##1.6##
neut.sim.uni<-function (M = 7, K = 500, S = 50, pop = round(rep(max(K)/S, M * M * S)), 
    s = .1, u = 0.01, d = 0.5, b = 0.505, m = 0.5 * u, habitat = 1, 
    fitness = 1, ef.strength = 0, recruit = 1, fec = 1, 
    time = 100, cycles = 20, mig=0,grad.len=NULL){
    ##Purpose: to create communities through time
    ##this function was modified from Smith's function 'neut.simulate'
    ##only differences are that here 'sim.init.uni' is called instead of 'sim.init' 
    ##and we added the argument 's' which specifies niche width or std.dev of Guassian function
    ##'mig' indicates if old (=0) or new (=1) migration algo should be used
    ##see ?neut.simulate for arguments
    ptm <- proc.time()
    sim <- sim.init.uni(M = M, K = K, S = S, pop = pop, s=s, u = u, d = d, 
         b = b, m = m, habitat = habitat, fitness = fitness, ef.strength = ef.strength, 
         recruit = recruit, fec = fec, grad.len=grad.len)
    if(mig==0) ##old style migration
     sim <- sim.cycle(sim = sim, time = time, cycles = cycles)
    if(mig==1) ##new style migration
     sim <- sim.cycle.uni(sim = sim, time = time, cycles = cycles)  
    print("system.time")
    print(proc.time() - ptm)
    return(sim)
}


##1.7##
sim.init.uni<-function (M = 7, K = 500, S = 50, pop = round(rep(K/S, M * M * S)),
    s= 0.1, u = 0.01, d = 0.5, b = 0.505, m = 0.5 * u, habitat = 1, 
    fitness = 1, ef.strength = 0, recruit = 1, fec = 1,grad.len=NULL){
    ##Purpose: to set-up an object for input into function 'sim.cycle'
    ##this function was modified from Smith's function 'sim.init'
    ##differences are that here birth and death rates 
    ##are modeled as a Gauassian function rather than a breakpoint function
    ##and the coordinates are simply the euclidean coordinates
    ##and 's' was added as an argument to specify the std.dev of Guassian niche
    ##see ?sim.init for arguments
    #if (sum(pop) != M^2 * K[1]) 
    #    warning("Cell populations not at carrying capacity")
    ##Note DJM changes here: scaling u and s by the mag. of the linear axis of the grid 
    u <- u * M 
    if(is.null(grad.len))
     s <- s * M
    else
     s <- s * grad.len
    ##Note DJM changes here: fitness may have length 1 to tolerate assemblages
    ##that start with just one species
    if (length(habitat) > 1) {
        fd <- array(dim = S * M^2)
        dim(fd) <- c(S, M, M)
        ##Note DJM's changes here - dropped abs on diff calculation ##NEED TO CHECK
        for (i in 1:S) for (j in 1:M) for (k in 1:M) fd[i, j,k] <- (habitat[j, k] - fitness[i])
        fd <- as.vector(fd)
    }
    else {
        fd <- rep(0, S * M * M)
    }
    ##Note DJM's changes here
    ##birth rate 
    bv <- gauss.niche.diff(b,fd,s) 
    ##death rate
    dv <- 1-gauss.niche.diff(1-d,fd,s) 
    ef <- recruit - fd * ef.strength
    params <- list(K = K, M = M, S = S, s.rel = s, u.rel = u, dv = dv, bv = bv, 
        m = m, habitat = habitat, fitness = fitness, ef = ef, fec = fec)
    ##Note DJM's changes here
    coords <- data.frame('X' = rep(1:M, each = M), 'Y' = rep(1:M, M))
    ret <- list(p = params, coords = coords, snaps = list(`0` = pop))
    class(ret) <- "sim"
    return(ret)
}

##1.8##
sim.cycle.uni<-function (sim, time = 500, cycles = 20){
    for (i in 1:cycles) {
        cat(paste("Cycle", i,";", sep = " "))
        sim <- sim.iter.uni(sim, time)
    }
    return(sim)
}



##1.9##
sim.iter.uni<-function (sim, time){
    nb = rep(0, sim$p$M * sim$p$M * sim$p$S)
    K = sim$p$K
    if(length(K)==1) ##if only a single K value supplied
     K = rep(K,sim$p$M * sim$p$M * sim$p$S) 
    prelength <- length(sim$snaps)
    pretime <- as.numeric(names(sim$snaps)[prelength])
    tmp <- .C("smith08c", M = as.integer(sim$p$M), K = as.integer(K), 
        S = as.integer(sim$p$S), pop = as.integer(sim$snaps[[prelength]]), 
        nb = as.integer(nb), u = as.double(sim$p$u), dv = as.double(sim$p$dv), 
        bv = as.double(sim$p$bv), m = as.double(sim$p$m), t = as.integer(time), 
        ef = as.double(sim$p$ef), fec = as.integer(sim$p$fec))
    sim$snaps[[prelength + 1]] <- tmp$pop
    names(sim$snaps)[prelength + 1] <- as.character(pretime +  time)
    return(sim)
}


##1.91##
D.vario.est<-function(z){
 ##esimates the fractal dimension D from a 2-dimensional grid
 ##z is a matrix of real numbers
 n<-dim(z)[1]
 gcords<-expand.grid(1:n,1:n)
 v<-vario(as.vector(z),gcords)$vario
 mod<-lm(log(v$exp)~log(v$Dist))
 m<-coef(mod)[2]
 D<-(6-m)/2
 D
}

#####Part II - FUNCTIONS FOR Spatial Permutations #####

##2.1##
SpatPerm2D<-function(psp,shiftpos=NULL,rotate=NULL,meth='shift',sp=FALSE){
 ##Purpose: to permute an array of occurances under a given set of constraints in 2-dimensions of space
 ##Arguments:
 ##psp is the sp x row x col array, where rows and columns specify where on the spatial grid the sample was located
 ##shiftpos: two numbers that are the x and y places to shift the grid, this is generated randomly if needed
 ##rotate: a single number 1-4 that indicates how many counterclockwise rotations to perform, generated randomly
 ##meth the type of permutation to use, options include:
 ###"reflect": random reflection/rotations of species (only makes sence when sp are not fixed
 ###"shift": random torodial shifting with or with sp fixed
 ###"both": both reflection and shifting
 ###"random": random shuffle
 ##if 'sp' is FALSE then obs composition of quadrats is fixed to the observed pattern
 ##if 'sp' is TRUE then species are each shuffled independently
 n<-dim(psp)[2]
 if(length(dim(psp))==3){
  S<-dim(psp)[1]
  flag<-FALSE
 }
 else{
  S<-1
  psp<-array(psp,dim=c(S,n,n))
  flag<-TRUE
 }
 Rpsp<-psp
 if(sp){##then between sp associations nullified
  if(meth!='reflect'){
   if(meth!='random'){
    ##generate vectors of random shifts, one for the x- and one for y-coord
    for(j in 1:S){
     if(is.null(shiftpos)){
      shift.x <- sample(n,size=1) ; shift.y <- sample(n,size=1);
     }
     else{
      shift.x <- shiftpos[1] ; shift.y <- shiftpos[2]
     }
     #gen new coords
     if(shift.x==1) ncoord.x <- 1:n
     else ncoord.x<-c(shift.x:n,1:(shift.x-1))
     if(shift.y==1) ncoord.y <- 1:n
     else ncoord.y<-c(shift.y:n,1:(shift.y-1))
     if(meth=='shift'){
      ##begin rearranging the rows of matrix for jth sp
      Rpsp[j,,] <- psp[j,ncoord.x,ncoord.y]
     }
     if(meth=='both'){##reflecting/rotating and shifting
      if(is.null(rotate))
       rotate<-sample(4,size=1)##how many counterclockwise rotatations to make
      if(rotate==2){
       for(x in 1:n){
        for(y in 1:n){
         Rpsp[j,(n-y)+1,x] <- psp[j,x,y]
      }}}
      if(rotate==3){
       for(x in 1:n){
        for(y in 1:n){
         Rpsp[j,(n-x)+1,(n-y)+1] <- psp[j,x,y]
      }}}
      if(rotate==4){
       for(x in 1:n){
        for(y in 1:n){
         Rpsp[j,y,(n-x)+1] <- psp[j,x,y]
      }}}
      flips<-sample(2,replace=TRUE) ##generates two coin flips
      if(flips[1]==1){ #reflect along x-axis
       if(flips[2]==1) #reflect along y-axis
        Rpsp[j,n:1,n:1] <- Rpsp[j,ncoord.x,ncoord.y]
       else #not reflected along y-axis
        Rpsp[j,n:1,] <- Rpsp[j,ncoord.x,ncoord.y]
      }
      else{ #not reflected along x-axis
       if(flips[2]==1) #reflected along y-axis
        Rpsp[j,,n:1] <- Rpsp[j,ncoord.x,ncoord.y]
       else ##not reflected along either axis
        Rpsp[j,,] <- Rpsp[j,ncoord.x,ncoord.y]
      }
  }}}}
  if(meth=='reflect'){##if only want reflecting/rotating
   for(j in 1:S){
    if(is.null(rotate))
     rotate<-sample(4,size=1)##how many counterclockwise rotatations to make
    if(rotate==2){
     for(x in 1:n){
      for(y in 1:n){
       Rpsp[j,(n-y)+1,x] <- psp[j,x,y]
    }}}
    if(rotate==3){
     for(x in 1:n){
      for(y in 1:n){
       Rpsp[j,(n-x)+1,(n-y)+1] <- psp[j,x,y]
    }}}
    if(rotate==4){
     for(x in 1:n){
      for(y in 1:n){
       Rpsp[j,y,(n-x)+1] <- psp[j,x,y]
    }}}
    flips<-sample(2,replace=TRUE) ##generates two coin flips
    if(flips[1]==1){ #reflect along x-axis
     if(flips[2]==1) #reflect along y-axis
      Rpsp[j,n:1,n:1] <- Rpsp[j,,]
     else #not reflected along y-axis
      Rpsp[j,n:1,] <- Rpsp[j,,]
    }
    else{ #not reflected along x-axis
     if(flips[2]==1) #reflected along y-axis
      Rpsp[j,,n:1] <- Rpsp[j,,]
    }
  }}
  if(meth=='random'){
   for(j in 1:S){
    take<-sample(n^2) #sample w/o replacement
    Rpsp[j,,] <- matrix(psp[j,,][take],ncol=n,nrow=n)
 }}}
 else{##species co-occurances are fixed
  ##this only makes sense for "shift" or "both" meth
  if(meth=="reflect"){
   stop("Reflecting fixed species co-occurances w/o shifting is not meaningful")
  }
  ##generate vector of random shifts
  if(meth!="random"){
   if(is.null(shiftpos)){
    shift.x <- sample(n,size=1) ; shift.y <- sample(n,size=1);
   }
   else{
    shift.x <- shiftpos[1] ; shift.y <- shiftpos[2]
   }
   #gen new coords
   if(shift.x==1) ncoord.x <- 1:n
   else ncoord.x<-c(shift.x:n,1:(shift.x-1))
   if(shift.y==1) ncoord.y <- 1:n
   else ncoord.y<-c(shift.y:n,1:(shift.y-1))
   if(meth=='shift'){
    ##begin rearranging the rows of matrix for jth sp
    Rpsp <- psp[,ncoord.x,ncoord.y]
   }
   if(meth=='both'){##reflecting/rotating and shifting
    if(is.null(rotate))
     rotate<-sample(4,size=1)##how many counterclockwise rotatations to make
    if(rotate==2){
     for(x in 1:n){
      for(y in 1:n){
       Rpsp[,(n-y)+1,x] <- psp[,x,y]
    }}}
    if(rotate==3){
     for(x in 1:n){
      for(y in 1:n){
       Rpsp[,(n-x)+1,(n-y)+1] <- psp[,x,y]
    }}}
    if(rotate==4){
     for(x in 1:n){
      for(y in 1:n){
       Rpsp[,y,(n-x)+1] <- psp[,x,y]
    }}}
    if(sample(2,size=1)==1) ##equivalent to a coin flip, if 1 then reflect and shift
     Rpsp[,n:1,n:1] <- psp[,ncoord.x,ncoord.y]
    else ##just shift
     Rpsp <- psp[,ncoord.x,ncoord.y]
  }}
  else{#meth is random and sp columns are fixed
   take<-sample(n^2) #sample w/o replacement
   for(j in 1:S){
    Rpsp[j,,] <- matrix(psp[j,,][take],ncol=n,nrow=n)
 }}}
 if(flag)
  Rpsp<-drop(Rpsp)
 Rpsp
}

##2.2##
SpatPerm2D.str<-function(psp,shiftpos=NULL,rotate=NULL,meth='shift',sp=FALSE,nstrata=1){
 ##Purpose: to permute an array of occurances under a given set of constraints in 2-dimensions of space
 ##with defined spatial strata, see 'nstrata' argument below
 ##Arguments:
 ##psp is the sp x row x col array, where rows and columns specify where on the spatial grid the sample was located
 ##shiftpos: two numbers that are the x and y places to shift the grid, this is generated randomly if needed
 ##rotate: a single number 1-4 that indicates how many counterclockwise rotations to perform, generated randomly
 ##meth the type of permutation to use, options include:
 ###"reflect": random reflection/rotations of species (only makes sence when sp are not fixed
 ###"shift": random torodial shifting with or with sp fixed
 ###"both": both reflection and shifting
 ###"random": random shuffle
 ##if 'sp' is FALSE then obs composition of quadrats is fixed to the observed pattern
 ##if 'sp' is TRUE then species are each shuffled independently
 ##'nstrata' is the number of strata along a single spatial axis within which to randomize
 n<-dim(psp)[2]
 if(length(dim(psp))==3){
  S<-dim(psp)[1]
  flag<-FALSE
 }
 else{
  S<-1
  psp<-array(psp,dim=c(S,n,n))
  flag<-TRUE
 }
 strata.size<-n/nstrata
 if(round(strata.size)!=strata.size)
  stop('Number of strata must be evenly divisable by the linear dimension of the grid')
 ##now simply apply the function SpatPerm2D on subsets of the orginal matrix and append all the pieces together at end
 Rpsp<-psp
 brks<-seq(1,n,strata.size)
 for(i in 1:nstrata){
  sub.i<-brks[i]:(brks[i]+strata.size-1)
  for(j in 1:nstrata){
   sub.j<-brks[j]:(brks[j]+strata.size-1)
   Rpsp[,sub.i,sub.j]<-SpatPerm2D(psp[,sub.i,sub.j],meth=meth,sp=sp,shiftpos=shiftpos,rotate=rotate)
 }}
 if(flag)
  Rpsp<-drop(Rpsp)
 Rpsp
}


##2.3##
RandPat<-function(i,psp,rpsp,n,nstrata,pl,mtrials1=1e3,mtrials2=1e6,alpha=0.01){
 ##Purpose: to be called in serial or parallel by function "RandPatPar"
 ##this function evaulates the .C function 'randpatpar' which is the random patterns algo of
 ##Roxburgh and Chesson 1998. Returns species index, phi stat, number of actual swaps, and the 
 ##randomized presences as a single vector of numbers
 ##Arguments:
 ##i: the ith species index
 ##psp: multidimenstional S x (n+2) x (n+2) array
 ##rpsp: a randomized version of psp
 ##n: the size of the orginal 2-D array along one spatial axis (i.e., without extra rows and columns)
 ##pl: the places in rpsp that can be swaped
 ##mtrials1: the number of times to attempt a swap at the strata level
 ##mtrials2: the number of times to attempt a swap at the pixel level
 ##alpha: the cutoff value for the phi statistic of Roxburgh and Chesson 1998
 ###################################################
 psp<-psp[i,,]
 rpsp<-rpsp[i,,]
 n2<-n+2
 ##PART I##
 ##begin permuting the blocks defined by nstrata
 rpsp.tmp<-rpsp[-c(1,n2),-c(1,n2)]
 coords<-cbind(rep(1:nstrata,nstrata),rep(1:nstrata,each=nstrata))
 rcoords<-coords[sample(nstrata^2),]
 for(j in 1:nstrata^2){
  rows<-((coords[j,1]-1)*n/nstrata+1) : ((coords[j,1]-1)*n/nstrata+n/nstrata)
  cols<-((coords[j,2]-1)*n/nstrata+1) : ((coords[j,2]-1)*n/nstrata+n/nstrata)
  rrows<-((rcoords[j,1]-1)*n/nstrata+1) : ((rcoords[j,1]-1)*n/nstrata+n/nstrata)
  rcols<-((rcoords[j,2]-1)*n/nstrata+1) : ((rcoords[j,2]-1)*n/nstrata+n/nstrata)
  rpsp[-c(1,n2),-c(1,n2)][rows,cols]<-rpsp.tmp[rrows,rcols]
 }
 rpsp<-FixUnSamp(psp,rpsp)
 ostat<-.C("spatstat",as.double(as.vector(psp)),as.integer(n),as.double(rep(0,4)))[[3]]
 nstat<-.C("spatstat",as.double(as.vector(rpsp)),as.integer(n),as.double(rep(0,4)))[[3]]
 phi<-.C("calcphi",as.double(nstat),as.double(ostat),as.double(0))[[3]]
 ##now begin random swapping of blocks defined by strata 
 ntrials<-0 ; gtrials<-0
 rpsp.tmp1<-rpsp
 while(phi > alpha & ntrials < mtrials1){
  rpsp.tmp2<-rpsp.tmp1[-c(1,n2),-c(1,n2)]
  rcoords<-coords[sample(nstrata^2,2),]
  startrows<-((rcoords[1,1]-1)*n/nstrata+1) : ((rcoords[1,1]-1)*n/nstrata+n/nstrata)
  startcols<-((rcoords[1,2]-1)*n/nstrata+1) : ((rcoords[1,2]-1)*n/nstrata+n/nstrata)
  endrows<-((rcoords[2,1]-1)*n/nstrata+1) : ((rcoords[2,1]-1)*n/nstrata+n/nstrata)
  endcols<-((rcoords[2,2]-1)*n/nstrata+1) : ((rcoords[2,2]-1)*n/nstrata+n/nstrata)
  rpsp.tmp1[-c(1,n2),-c(1,n2)][startrows,startcols] <- rpsp.tmp2[endrows,endcols]
  rpsp.tmp1[-c(1,n2),-c(1,n2)][endrows,endcols] <- rpsp.tmp2[startrows,startcols]
  rpsp.tmp1<-FixUnSamp(psp,rpsp.tmp1)
  nstat<-.C("spatstat",as.double(as.vector(rpsp.tmp1)),as.integer(n),as.double(rep(0,4)))[[3]]
  phiTemp <- .C("calcphi",as.double(nstat),as.double(ostat),as.double(0))[[3]]
  if(phiTemp < phi){
   phi <- phiTemp
   gtrials <- gtrials +1
   ##make change permanent
   rpsp<-rpsp.tmp1
  }
  else{
   ##start back with orginal random mat
   rpsp.tmp1<-rpsp
  }
  ntrials<-ntrials+1
 } 
 if(phi > alpha & mtrials2>0){
  ##Part II##
  ##carry out individual pixel swapping
  psp<-as.vector(psp)
  rpsp<-as.vector(rpsp)
  tmp<-.C("randpatpar",psp = as.double(psp),
        rpsp = as.double(rpsp), n = as.integer(n),
         ostat = as.double(rep(0,4)), nstat = as.double(rep(0,4)), 
         phi = as.double(0), phiTemp = as.double(0),
         alpha = as.double(alpha), pl = as.integer(pl-1),
         nplaces = as.integer(length(pl)-1),ntrials = as.double(0),
         gtrials = as.double(0), mtrials = as.double(mtrials2))
  out<-c(i,phi,gtrials,tmp$phi,tmp$gtrials,tmp$rpsp)
 }
 else
  out<-c(i,phi,gtrials,NA,NA,as.vector(rpsp))
 out
}

##2.4##
RandPatPar<-function(psp,nstrata,mtrials1=1e3,mtrials2=1e6,alpha=0.01,npar=1){
 ##Purpose: convience function for working with RandPat which calls the .C function
 ##'randpatpar'. This function allows you to specify the number of processors to run on
 ##adding processsors only helps if working with many species as each species is evaulated
 ##on a different processor. Returns a (5+(n+2)^2) x S matrix, the first five rows of which 
 ##are species index, phi strata stat, number of strata swaps, phi pixel swap, and number of 
 ##pixel swaps, and then the remaining rows are the presences/abundances in the randomized occurances
 ##Arguments:
 ##psp: multidimenstional S x (n+2) x (n+2) array
 ##n: the size of the orginal 2-D array along one spatial axis (i.e., without extra rows and columns)
 ##pl: the places in rpsp that can be swaped
 ##mtrials: the numbef of times to attempt a swap
 ##alpha: the cutoff value for the phi statistic of Roxburgh and Chesson 1998
 ##npar: the number of processors to run the code on
 ###################################################
 n<-dim(psp)[2]
 if(length(dim(psp))==3)
  S<-dim(psp)[1]
 else{
  S<-1
  psp<-array(psp,dim=c(S,n,n))
 }
 ####first prepare psp for the randomization process####
 ##fill in empty pixels with -999##
 #sampled<-rep(ifelse(apply(psp,c(2,3),sum)>0,0,-999),each=S)
 #dim(sampled)<-c(S,n,n) 
 #psp <- psp + sampled
 ##add border of -999###
 n2<-n+2
 psp.temp<- array(0,dim=c(S,n2,n2)) ##create an array with boundary cells
 psp.temp[,-c(1,n2),-c(1,n2)]<-psp ##populate the array with the input data
 psp.temp[,1,-c(1,n2)]<--999 ##x of 0
 psp.temp[,n2,-c(1,n2)]<--999 ##x of n+1
 psp.temp[,-c(1,n2),1]<--999##y of 0
 psp.temp[,-c(1,n2),n2]<--999 ##y of n+1
 psp.temp[,1,1]<--999 ;  psp.temp[,1,n2]<--999 ;  psp.temp[,n2,1]<--999 ;  psp.temp[,n2,n2]<--999
 psp<-psp.temp
 pl<-1:n2^2
 c1<-NA
 c2<-NA
 skip<-pl[-999==as.vector(psp[1,,])]
 pl<-pl[is.na(match(pl,skip))]
 ##now read to begin intital randomization
 rpsp<-psp
 ##inital reflection/rotation for each species independently
 rpsp[,-c(1,n2),-c(1,n2)]<-SpatPerm2D.str(psp[,-c(1,n2),-c(1,n2)],meth='reflect',sp=TRUE,nstrata=nstrata)
 nplaces <- length(pl)
 nswaps <- (nplaces*(nplaces-1))/2
 if(npar>1){
  require(snowfall)
  sfInit(parallel=TRUE, cpus=npar, type="SOCK")
  sfClusterSetupRNG()
  sfExport("RandPat", "FixUnSamp", "psp", "rpsp","n","nstrata","pl","mtrials1","mtrials2","alpha")
  sfClusterEval(dyn.load("randompatternspar.dll"))
  out<-unlist(sfSapply(1:S,function(i)
    RandPat(i=i,psp=psp,rpsp=rpsp,n=n,nstrata=nstrata,pl=pl,mtrials1=mtrials1,mtrials2=mtrials2,alpha=alpha)))
  sfStop()
 }
 else{
  out<-NULL
  for(i in 1:S)
   out<-cbind(out,RandPat(i=i,psp=psp,rpsp=rpsp,n=n,nstrata=nstrata,pl=pl,mtrials1=mtrials1,mtrials2=mtrials2,alpha=alpha))
 }
 out
}

##2.5##
FixUnSamp<-function(oarray,rarray){
 ##purpose: to maintain the spatial locations
 ##of the unsampled pixels in rarray which is a random
 ##realization of oarray, -999 is the identifier for 
 ##unsampled cells, in this case oarray and rarray DO have a false border of -999
 rarray.tmp<-rarray
 if(length(dim(oarray))==3){ ##if multiple species then
  n2<-dim(oarray)[2]
  if(-999%in%oarray[1,-c(1,n2),-c(1,n2)]){ ##if there are unsampled pixels in the data, then
   S<-dim(oarray)[1]
   o.na<-oarray==-999
   r.na<-rarray==-999
   end.tmp<-which(o.na[1,-c(1,n2),-c(1,n2)]) ##identifies where in o.na it is -999
   for(i in 1:S){
    start.tmp<-which(r.na[i,-c(1,n2),-c(1,n2)])##identifies where in r.na it is -999
    if(sum(which(!is.na(match(end.tmp,start.tmp))))>0){
     ##drop ones in which end.tmp and start.tmp match
     end<-end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
     start<-start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
    }
    else{
     end <- end.tmp
     start <-start.tmp
    }
    rarray[i,-c(1,n2),-c(1,n2)][end]<- rarray.tmp[i,-c(1,n2),-c(1,n2)][start]
    rarray[i,-c(1,n2),-c(1,n2)][start]<- rarray.tmp[i,-c(1,n2),-c(1,n2)][end]
 }}}
 else{ ##only a single species
  n2<-dim(oarray)[1]
  if(-999%in%oarray[-c(1,n2),-c(1,n2)]){ ##if there are unsampled pixels in the data, then
   o.na<-oarray==-999
   r.na<-rarray==-999
   end.tmp<-which(o.na[-c(1,n2),-c(1,n2)]) ##identifies where in o.na it is -999
   start.tmp<-which(r.na[-c(1,n2),-c(1,n2)])##identifies where in r.na it is -999
   if(sum(which(!is.na(match(end.tmp,start.tmp))))>0){
    ##drop ones in which end.tmp and start.tmp match
    end<-end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
    start<-start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
   }
   else{
    end <- end.tmp
    start <-start.tmp
   }
   rarray[-c(1,n2),-c(1,n2)][end]<- rarray.tmp[-c(1,n2),-c(1,n2)][start]
   rarray[-c(1,n2),-c(1,n2)][start]<- rarray.tmp[-c(1,n2),-c(1,n2)][end]
 }}
 rarray
}

##2.6##
FixUnSamp2<-function(oarray,rarray){
 ##purpose: to maintain the spatial locations
 ##of the unsampled pixels in rarray which is a random
 ##realization of oarray, -999 is the identifier for 
 ##unsampled cells, in this case oarray and rarray DO NOT have a false border of -999
 if(-999%in%oarray){ ##if there are unsampled pixels in the data, then
  rarray.tmp<-rarray
  if(length(dim(oarray))==3){ ##if multiple species then
   S<-dim(oarray)[1]
   n<-dim(oarray)[2]
   o.na<-oarray==-999
   r.na<-rarray==-999
   end.tmp<-which(o.na[1,,]) ##identifies where in o.na it is -999
   for(i in 1:S){
    start.tmp<-which(r.na[i,,])##identifies where in r.na it is -999
    if(sum(which(!is.na(match(end.tmp,start.tmp))))>0){
     ##drop ones in which end.tmp and start.tmp match
     end<-end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
     start<-start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
    }
    else{
     end <- end.tmp
     start <-start.tmp
    }
    rarray[i,,][end]<- rarray.tmp[i,,][start]
    rarray[i,,][start]<- rarray.tmp[i,,][end]
  }}
  else{ ##only a single species
   n2<-dim(oarray)[1]
   o.na<-oarray==-999
   r.na<-rarray==-999
   end.tmp<-which(o.na) ##identifies where in o.na it is -999
   start.tmp<-which(r.na)##identifies where in r.na it is -999
   if(sum(which(!is.na(match(end.tmp,start.tmp))))>0){
    ##drop ones in which end.tmp and start.tmp match
    end<-end.tmp[-which(!is.na(match(end.tmp,start.tmp)))]
    start<-start.tmp[-which(!is.na(match(start.tmp,end.tmp)))]
   }
   else{
    end <- end.tmp
    start <-start.tmp
  }
  rarray[end]<- rarray.tmp[start]
  rarray[start]<- rarray.tmp[end]
 }}
 rarray
}
#####Part III - ANALYZING AND GRAPHING RESULTS#####

##3.1##
'getCovFractions' = function(x)
{
  ## Purpose: calculates the lower diagonal of a sp covariance matrix
  ## to provide the positive and negative fractions of covariance
  ## output is two lower triangular matrices, each in vector format
  ## Called within the function 'vario'
  ## Arguments: 
  ## x is a sitexsp matrix (sp as columns) of real numbers
  ## rows are the sites, columns are the species
  N = as.integer(nrow(x))
  S = as.integer(ncol(x)) 
  x = as.double(ifelse(is.na(x) | x == -999,-99999,x))
  pos = as.double(rep(0,(N*(N-1))/2))
  neg = as.double(rep(0,(N*(N-1))/2))
  result = .C('loopcovreal',x,N,S,pos,neg,PACKAGE = vario)
  out = list()
  out$pos = result[[4]]
  out$neg = result[[5]]
  return(out)
} 

##3.2##
'vario' = function(x, coord, grain=1, breaks=NA, hmin=NA, hmax=NA, round.int=FALSE,
                   pos.neg=FALSE, binary=TRUE, snap=NA, median=FALSE, 
                   quants=NA, direction = 'omnidirectional',
                   tolerance = pi/8, unit.angle = c('radians', 'degrees'),
                   distance.metric = 'euclidean')
{
  ## Purpose: calculates uni- and multi-variate variograms
  ##
  ## This code is largely modified from the 'vegan' function 'mso' by Helene
  ## Wagner, also parts of this code were derived from the 'geoR' function 'vairog'
  ## by Paulo J. Ribeiro Jr. and Peter J. Diggle.
  ## Citations:
  ## Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre,
  ## Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M.
  ## Henry H. Stevens and Helene Wagner (2011). vegan: Community Ecology
  ## Package. R package version 2.0-2.
  ## http://CRAN.R-project.org/package=vegan
  ## Paulo J. Ribeiro Jr & Peter J. Diggle geoR: a package for
  ## geostatistical analysis R-NEWS, 1(2):15-18. June, 2001
  ##
  ## Note: that if some areas are unsampled that obs != pos + neg + exp
  ## Arguments:
  ## x: an object of class 'sim' or may be a sitexsp matrix, a vector of values,
  ##   if x is the later then missing samples should be coded as -999
  ## coord: the spatial coordinates
  ## grain: interval size for distance classes, only used if 'breaks' not supplied
  ## breaks: the spatial breaks that define the spatial lags to be compared
  ## hmin: the minimum spatial lag, default value of NA is treated as a minimum
  ##   of 1
  ## hmax: the maximum spatial lag, default value of NA is treated as half of 
  ##   the maximum distance
  ## round.int: if TRUE the spatial lags are rounded to nearest integer
  ## pos.neg: if TRUE the positive and negative parts of the covariance matrix
  ##   are output
  ## binary: if TRUE and x is class sim then abundances are converted to
  ## binary 0 1 values
  ## snap: indicates which generation from an object class 'sim' to draw from
  ## median: indicates if in addition to the mean the medians of the distance
  ##   matrices are calculated
  ## direction: a numerical value for the directional (azimuth) angle. This
  ##   used to specify directional variograms. Default defines the
  ##   omnidirectional variogram. The value must be in the interval [0, pi] 
  ##   radians ([0, 180] degrees).
  ## quants: any quantiles to compute, these offer rough estimates of variability
  ##   in the empirical variogram
  ## tolerance: numerical value for the tolerance angle, when computing
  ##   directional variograms. The value must be in the interval [0, pi/2]
  ##   radians ([0, 90] degrees). Defaults to pi/8.
  ## unit.angle: defines the unit for the specification of angles in the two
  ##   previous arguments. Options are 'radians' and 'degrees', with default to
  ##   'radians'.
  ## distance.metric': can be one of the speices turnover metrics listed by the
  ##   vegan function vegdist(). This is only appropriate if pos.neg = FALSE.
  ##   Common options include, 'jaccard' and 'bray'. If computed on pres/abse
  ##   data then soreson index is computed by 'bray'.
  ##----------------------------------------------------------------------------
  if (distance.metric != 'euclidean') {
    if (pos.neg)
      stop("cannot commpute pos-neg covariance using a turnover metric")
    else
      require(vegan)
  }
  ## geoR code from function variog starts here'
  unit.angle = match.arg(unit.angle)
  if (mode(direction) == "numeric") {
    if (length(direction) > 1)
      stop("only one direction is allowed")
    if (length(tolerance) > 1)
      stop("only one tolerance value is allowed")
    if (unit.angle == "degrees") {
      ang.deg = direction
      ang.rad = (ang.deg * pi) / 180
      tol.deg = tolerance
      tol.rad = (tol.deg * pi) / 180
    }
    else {
     ang.rad = direction
     ang.deg = (ang.rad * 180) / pi
     tol.rad = tolerance
     tol.deg = (tol.rad * 180) / pi
    }
    if (ang.rad > pi | ang.rad < 0)
      stop("direction must be an angle in the interval [0,pi[ radians")
    if (tol.rad > pi/2 | tol.rad < 0)
      stop("tolerance must be an angle in the interval [0,pi/2] radians")
    if (tol.deg >= 90) {
      direction = "omnidirectional"
    }
  }  
  ## geoR code from function variog ends here'
  if (class(x) == "sim"){
    coord = x$coords
    if (is.na(snap))
      snap = length(sim$snaps)
  }
  else
    x = ifelse(x == -999, NA, x)
  Dist = dist(coord)
  maxDist = max(Dist)
  if (is.na(breaks[1])) {
    if (is.na(hmin))
      hmin = grain
    if (is.na(hmax))
      hmax = round((maxDist / 2) / grain) * grain
    H = round(Dist / grain) * grain
  }
  else {
    if (is.na(hmin))
      hmin = 1
    if (is.na(hmax))
      hmax = maxDist / 2
    H = Dist
    if (is.numeric(breaks[1])) {
      if (length(breaks) == 1)
        breaks = seq(hmin, hmax, length.out=breaks)
    }
    else {
      if (breaks[1] == 'log')
        base = 'exp'
      else if (breaks[1] == 'log2')
        base = '2^'
      else if (breaks[1] == 'log10')
        base = '10^'
      else 
        stop('Specification of breaks using a character string must be log, log2, 
             or log10')
      incre = (hmax - hmin) / as.numeric(breaks[2])
      hmax = hmax + incre
      breaks = eval(parse(text = paste(base, '(seq(', breaks[1], '(hmin),', 
                                       breaks[1], '(hmax), length.out=', 
                                       breaks[2], '))', sep='')))
    } 
    if (round.int)
      breaks = round(breaks)
    for (i in 1:(length(breaks) - 1)) {
      H[H >= breaks[i] & H < breaks[i + 1]] = breaks[i]
    }  
  }
  H[H < hmin] = NA
  H[H > hmax] = NA
  H = as.vector(H)
  if (is.vector(x)) {
    S = 1
    N = length(x)
  }
  else {
    S = ncol(x)
    N = nrow(x)
  } 
  vobject = list()
  vobject$parms = data.frame(grain, hmax, S=S, N=N, pos.neg, median, direction,
                             tolerance, unit.angle, distance.metric, 
                             quants = ifelse(is.na(quants[1]), NA, 
                                             paste(quants* 100, collapse=", ")))
  if(class(x) == "sim"){
    if(binary)
      x = apply(census(x, snap=snap), c(1, 2), as.logical) * 1
    else
      x = census(x, snap=snap)
    vobject$parms = cbind(vobject$parms, niche.wid.rel=x$p$s.rel,
                          disp.wid.rel=x$p$u.rel) 
  }
  ## geoR code from function variog with slight modifications starts here'
  if (direction != "omnidirectional") {
    ## note that the changes: 'u' has been changed for as.vector(Dist) and
    ## coords changed to coord
    u.ang = .C("tgangle", as.double(as.vector(coord[ , 1])),
               as.double(as.vector(coord[ , 2])), as.integer(dim(coord)[1]),
               res = as.double(rep(0, length(as.vector(Dist)))), 
               PACKAGE = vario)$res
    if (any(is.na(u.ang)))
      stop("NA returned in angle calculations maybe due to co-located data")
    u.ang = atan(u.ang)
    u.ang[u.ang < 0] = u.ang[u.ang < 0] + pi
    ang.lower = ang.rad - tol.rad
    ang.upper = ang.rad + tol.rad
    if (ang.lower >= 0 & ang.upper < pi)
      ang.ind = (!is.na(u.ang) & ((u.ang >= ang.lower) & (u.ang <= ang.upper)))
    if (ang.lower < 0)
      ang.ind = (!is.na(u.ang) & ((u.ang < ang.upper) | (u.ang > (pi + ang.lower))))
    if (ang.upper >= pi)
      ang.ind = (!is.na(u.ang) & ((u.ang > ang.lower) | (u.ang < (ang.upper - pi))))
    Dist[!ang.ind] = NA
    H[!ang.ind] = NA
  }
  ## geoR code from function variog ends here
  Dist = sapply(split(Dist, H), mean, na.rm=TRUE)
  vobject$vario = data.frame(H = as.numeric(names(table(H))), Dist = Dist,
                             n = as.numeric(table(H)))
  ## below 'exp.gamma' is the expected variogram if 'x' is a sitexsp
  ## pres/abse matrix. The expectation is based upon the assumption of zero
  ## sp x sp covariances
  if(distance.metric == 'euclidean')
    exp.split = split(dist(x)^2 * .5, H)
  else
    exp.split = split(vegdist(x, method=distance.metric), H)
  exp.gamma = sapply(exp.split, mean, na.rm=TRUE)
  if (!is.na(quants[1])) {
    exp.qt = sapply(exp.split, function(x) quantile(x, quants, na.rm=TRUE))
    exp.qt = t(exp.qt)
    colnames(exp.qt) = paste(quants * 100)
  }  
  vobject$vario = cbind(vobject$vario, exp.var=exp.gamma)
  if (median)
    exp.med = sapply(exp.split, median, na.rm=TRUE)
  if (!is.vector(x)) { ## i.e. x is a site x sp matrix and not simply a vector
    ## if 'x' is a sitexsp pres/abse matrix the following computes site species richness
    rich = apply(x, 1, sum)
    ## see equation 7 in Wagner, H. 2003. Spatial covariance in plant
    ## communities... Ecology 84:1045-1057 to see that observed multivariate
    ## variogram can be computed from the species richness vector
    obs.gamma = sapply(split(dist(rich)^2 * .5, H), mean, na.rm=TRUE)
    vobject$vario = cbind(vobject$vario, obs.var = obs.gamma,
                          ratio = obs.gamma/exp.gamma)
    if (pos.neg) {
      cov.mat = getCovFractions(x)
      pos.split = split(cov.mat$pos, H)
      neg.split = split(cov.mat$neg, H)
      pos = sapply(pos.split, mean)
      neg = sapply(neg.split, mean)
      vobject$vario = cbind(vobject$vario, pos = pos, neg = neg)
      if (median) {
        pos.med = sapply(pos.split, median)
        neg.med = sapply(neg.split, median)
        vobject$vario = cbind(vobject$vario, exp.med = exp.med, pos.med = pos.med,
                              neg.med = neg.med)
      }
    ##Note: obs.var = exp.var + pos.var + neg.var
    }
    else{
      if (median) {
        obs.gamma.med = sapply(split(dist(rich)^2 * .5, H),median, na.rm=TRUE)
        vobject$vario = cbind(vobject$vario, obs.med = obs.gamma.med,
                              exp.med = exp.med)
      }  
    }
  }
  if (!is.na(quants[1]))
    vobject$vario = cbind(vobject$vario, exp.qt = exp.qt)
  if (is.vector(x))
    vobject$p = sum(x) / length(x)
  else
    vobject$p = apply(x, 2, sum, na.rm=TRUE) / nrow(x)
  vobject$perm = FALSE
  ## the above line indicates to the function 'vGraph' if the variogram is the
  ## result of a randomization or not
  row.names(vobject$vario) = NULL
  return(vobject)
}

##3.3##
null.perms<-function(x,vobject,nperm,coords=NULL,meth='both',sp=TRUE,all=FALSE,snap=NULL,npar=1,linux=FALSE,RPargs=FALSE){
 ##Purpose: to generate statistical null expectations for the variograms
 ##Arguments:
 ##"x" is either an output of class 'sim' that is the output of 'sim.neut.uni' OR an site x species matrix
 ##"vobject" is the output of 'vario', this informs the function of what parameterization of vario to use
 ###specifically it indiates if the pos.neg components and median should be calculated
 ##"nperm" is the number of permutations
 ##"coords" the spatial coordinates of the sites, not needed if x is of class 'sim'
 ##"meth" the type of permutation to use, options include:
 ##"reflect": random reflection/rotations of species (only makes sence when sp are not fixed
 ###"shift": random torodial shifting with or with sp fixed
 ###"both": both reflection and shifting
 ###"random": random shuffle
 ###"randpat": random patterns algo of Roxburgh and Chesson 1998, must parameterize RPargs (See below)
 ##"sp" = FALSE then obs composition of quadrats is fixed to the observed pattern
 ##"sp" = TRUE then species are each shuffled independently
 ##"all" if TRUE then all relevant nulls calculated
 ##"snap" is the time period of the simulation to analyze,defaults to NULL val, if not set gets internally set to last time period
 ##"npar" = number of processors to run the function on
 ##"median" if TRUE then median is also calculated in addition to mean
 ##"linux" if TRUE then function assumes you are on a linux cluster and therefore exports a different compiled code
 ##"RPargs" is a vector of arguments that are needed to perform the Random Patterns spatial null model
 ###see the notes associated with the function 'null.gen' that indicate how "RPargs" should be parameterized
 ##Note: "meth" and "sp" are arguments to randomization function "SpatPerm2D"
 dists<-vobject$vario$Dist
 grain = vobject$parms$grain
 hmax = vobject$parms$hmax
 pos.neg = vobject$parms$pos.neg
 median = vobject$parms$median
 if(class(vobject$parms$direction) == "factor") 
  direction = as.character(vobject$parms$direction)
 else
  direction = as.numeric(vobject$parms$direction)
 tolerance = vobject$parms$tolerance
 unit.angle = as.character(vobject$parms$unit.angle)
 hmax<-vobject$parms$hmax
 if(class(x)=='sim'){
  coords<-x$coords
  if(is.null(snap)) snap <- length(x$snaps)
 }
 else{
  if(is.null(coords))
   stop('need to supply spatial coordinates if not a simulation product')
 }
 r.vals<-list()
 r.vals$parms<-vobject$parms
 r.vals$p<-vobject$p
 if(median&!pos.neg)
  stop("if computing medians must also compute pos.neg fractions, set pos.neg=TRUE")
 if(pos.neg){ ##pos and neg fractions
  if(all){ ##all relevant null models
   if(median){ ##will compute mean and median
    r.vals$vario<-array(0,dim=c(length(dists),6,3,nperm+1))##dists,results,methods,perms
    r.vals$vario[,,,1]<-as.matrix(vobject$vario[,c(5,7:11)])   
   }
   else{ ##only compute means
    r.vals$vario<-array(0,dim=c(length(dists),3,3,nperm+1))##dists,results,methods,perms
    r.vals$vario[,,,1]<-as.matrix(vobject$vario[,c(5,7:8)])
  }}
  else{ ##only a single null used
   if(median){
    r.vals$vario<-array(0,dim=c(length(dists),6,nperm+1))##dists,results,perms
    r.vals$vario[,,1]<-as.matrix(vobject$vario[,c(5,7:11)])
   }
   else{
    r.vals$vario<-array(0,dim=c(length(dists),3,nperm+1))##dists,results,perms
    r.vals$vario[,,1]<-as.matrix(vobject$vario[,c(5,7:8)])
 }}}
 else{##only exp and obs fractions
  if(all){
   if(median){ ##will compute mean and median
    r.vals$vario<-array(0,dim=c(length(dists),4,3,nperm+1))
    r.vals$vario[,,,1]<-as.matrix(vobject$v[,c(4:5,7:8)])
   }
   else{
    r.vals$vario<-array(0,dim=c(length(dists),2,3,nperm+1))
    r.vals$vario[,,,1]<-as.matrix(vobject$v[,4:5])
  }}
  else{ ##only a single null used
   if(median){
    r.vals$vario<-array(0,dim=c(length(dists),4,nperm+1))
    r.vals$vario[,,1]<-as.matrix(vobject$v[,c(4:5,7:8)])
   }
   else{ 
    r.vals$vario<-array(0,dim=c(length(dists),2,nperm+1))
    r.vals$vario[,,1]<-as.matrix(vobject$v[,4:5])
 }}}
 if(class(x)=='sim'){
  pop<-as.logical(x$snaps[[snap]]) ##converts it to a pres/absence vector
  dim(pop) <- c(x$p$S, x$p$M, x$p$M)
 }
 else{
  pop <- array(x,dim=c(sqrt(nrow(x)),sqrt(nrow(x)),ncol(x)))
  pop <- aperm(pop,c(3,1,2))
 }
 if(RPargs[[1]]&npar==1){
  r.vals$p.conv1<-0 ##average proportion of species that converged with strata swaps
  r.vals$p.conv2<-0 ##average proportion of species that converged with pixel swaps
 }
 if(npar==1){ ##all permutations option not yet implemented for 1 processor
  pb <- txtProgressBar(min = 0, max = nperm, style = 3)
  for(i in 1:nperm){
   if(RPargs[[1]]){ ##use the random pattern algo for the spatial null
    out<-RandPatPar(psp=pop,nstrata=RPargs[[2]],mtrials1=RPargs[[3]],mtrials2=RPargs[[4]],alpha=RPargs[[5]],npar=RPargs[[6]])
    S<-dim(pop)[1]
    n2<-dim(pop)[2]+2
    r.vals$p.conv1 <- r.vals$p.conv1 + sum(out[2,]<=RPargs[[5]])/S/nperm
    r.vals$p.conv2 <- r.vals$p.conv2 + sum(out[4,]<=RPargs[[5]],na.rm=TRUE)/S/nperm
    rpop<-array(0,dim=c(S,n2,n2))
    for(k in 1:S){
     rpop[k,,]<-array(out[-(1:5),k],dim=c(n2,n2))
    }
    rpop<-rpop[,-c(1,n2),-c(1,n2)]
   }
   else{
    rpop<-SpatPerm2D(pop,meth=meth,sp=sp)
    rpop<-FixUnSamp2(pop,rpop)
   }
   rmat<-apply(rpop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
   rv<-vario(x=rmat,coord=coords,grain=grain,hmax=hmax,pos.neg=pos.neg,median=median,
             direction=direction,tolerance=tolerance,unit.angle=unit.angle)$vario
   if(pos.neg){
    if(all){
     if(median)
      r.vals$vario[,,,i+1] <- as.matrix(rv[,,c(5,7:11)])
     else
      r.vals$vario[,,,i+1] <- as.matrix(rv[,,c(5,7:8)])
    } 
    else{
     if(median)
      r.vals$vario[,,i+1] <- as.matrix(rv[,c(5,7:11)])
     else
      r.vals$vario[,,i+1] <- as.matrix(rv[,c(5,7:8)])
   }}
   else{
    if(all){
     if(median)
      r.vals$vario[,,,i+1] <- as.matrix(rv[,,c(4:5,7:8)])
     else
       r.vals$vario[,,,i+1] <- as.matrix(rv[,,4:5])
    }
    else{
     if(median)
      r.vals$vario[,,i+1] <- as.matrix(rv[,c(4:5,7:8)])
     else
      r.vals$vario[,,i+1] <- as.matrix(rv[,4:5])
   }}
   #print(i)
   Sys.sleep(0.1)
   # update progress bar
   setTxtProgressBar(pb, i)
  }
  close(pb)
 }
 else{ ##computing in parallel
  require(snowfall)
  sfInit(parallel=TRUE, cpus=npar, type="SOCK")
  sfClusterSetupRNG()
  sfExport("pop", "vobject", "coords", "meth", "all", "sp", "RPargs", "RandPatPar", "RandPat", "FixUnSamp","FixUnSamp2", "SpatPerm2D", "SpatPerm2D.str", "vario", "getCovFractions", "null.gen")
  if(linux)
   sfClusterEval(dyn.load("danspkg.so"))
  else{
   sfLibrary(danspkg)
  }
  out<- unlist(sfLapply(1:nperm,function(...)null.gen(pop=pop,vobject=vobject,coords=coords,meth=meth,sp=sp,all=all,RPargs=RPargs)))
  sfStop()
  if(pos.neg){
   if(all){
    if(median){
     dim(out)<-c(length(dists),6,3,nperm)
     r.vals$vario[,,,-1]<-out 
    }
    else{
     dim(out)<-c(length(dists),3,3,nperm)
     r.vals$vario[,,,-1]<-out 
   }}
   else{
    if(median){
     dim(out)<-c(length(dists),6,nperm)  
     r.vals$vario[,,-1]<-out     
    }
    else{
     dim(out)<-c(length(dists),3,nperm)  
     r.vals$vario[,,-1]<-out     
  }}}
  else{
   if(all){
    if(median){
     dim(out)<-c(length(dists),4,3,nperm)  
     r.vals$vario[,,,-1]<-out 
    }
    else{
     dim(out)<-c(length(dists),2,3,nperm)  
     r.vals$vario[,,,-1]<-out 
   }}
   else{
    if(median){
     dim(out)<-c(length(dists),4,nperm)  
     r.vals$vario[,,-1]<-out 
    }
    else{
     dim(out)<-c(length(dists),2,nperm)  
     r.vals$vario[,,-1]<-out 
 }}}}
 r.vals$perm<-TRUE
 r.vals$vdists<-vobject$vario$Dist
 r.vals
} 

##3.4##
null.gen<-function(pop,vobject,coords,meth,sp,all=FALSE,RPargs=FALSE,median=FALSE){
 ##Purpose: to mediate the generation of statistical null values for the variograms 
 ##to be used in a parrallel processing loop which will generate a population of null values
 ##Arguments:
 ##"pop" the species x row x col array where row and column refer to spatial location
 ##"vobject" the output of the vario function with serves as the basis for empirical comparisons
 ##'coords' the spatial coordinates of a M^2 x S matrix
 ##"meth" the type of permutation to use, options include:
 ###"reflect": random reflection/rotations of species (only makes sence when sp are not fixed
 ###"shift": random torodial shifting with or with sp fixed
 ###"both": both reflection and shifting
 ###"random": random shuffle
 ###"randpat": random patterns algo of Roxburgh and Chesson 1998, must parameterize RPargs (See below)
 ##"sp" = FALSE then obs composition of quadrats is fixed to the observed pattern
 ##"sp" = TRUE then species are each shuffled independently
 ##'all' if TRUE then all relevant nulls are calculated
 ##'RPargs' is a list of arguments that must be supplied if the random patterns algo is desired
 ###the arguments of RPargs are used input into the function 'RandPatPar', they include:
 ###1)'allRP' if TRUE & 'all' = TRUE, then Random Patterns algo used as the spatial null
 ###2)'nstrata',3)'mtrials1',4)'mtrials2',5)'alpha',6)'npar'
 ##'median' if TRUE then means and medians are calculated
 ##Note: "meth" and "sp" are arguments to randomization function "SpatPerm2D"
 S<-dim(pop)[1]
 n<-dim(pop)[2]
 n2<-n+2
 grain = vobject$parms$grain
 hmax = vobject$parms$hmax
 pos.neg = vobject$parms$pos.neg
 median = vobject$parms$median
 if(class(vobject$parms$direction) == "factor") 
  direction = as.character(vobject$parms$direction)
 else
  direction = as.numeric(vobject$parms$direction)
 tolerance = vobject$parms$tolerance
 unit.angle = as.character(vobject$parms$unit.angle)
 if(all){
  rmat<-apply(pop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
  rv<-vobject$vario
  if(pos.neg){
   if(median) 
    r.vals<-array(NA,dim=c(length(rv$Dist),6,3)) ##dists,6 results, 3 methods
   else
    r.vals<-array(NA,dim=c(length(rv$Dist),3,3)) ##dists,3 results, 3 methods
  }
  else
   r.vals<-array(NA,dim=c(length(rv$Dist),2,3)) ##dists, 2 results, 3 methods
   ###need to add in logicals for handling median values here
  for(j in 1:3){ #cycle through permutation methods
   if(j == 1){
    rpop<-SpatPerm2D(pop,meth='random',sp=TRUE)##random model
    rpop<-FixUnSamp2(pop,rpop)
   }
   if(j == 2){
    if(RPargs[[1]]){ ##use the random pattern algo for the spatial null
     out<-RandPatPar(psp=pop,nstrata=RPargs[[2]],mtrials1=RPargs[[3]],mtrials2=RPargs[[4]],alpha=RPargs[[5]],npar=RPargs[[6]])
     rpop<-array(0,dim=c(S,n2,n2))
     for(i in 1:S){
      rpop[i,,]<-array(out[-(1:5),i],dim=c(n2,n2))
     }
     rpop<-rpop[,-c(1,n2),-c(1,n2)]
    }
    else
     rpop<-SpatPerm2D(pop,meth='both',sp=TRUE)##Random Shifts & Reflections spatial null
   }
   if(j == 3){
    rpop<-SpatPerm2D(pop,meth='random',sp=FALSE) ##species model
    rpop<-FixUnSamp2(pop,rpop)
   }
   rmat<-apply(rpop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
   rv<-vario(x=rmat,coord=coords,grain=grain,hmax=hmax,pos.neg=pos.neg,median=median,
             direction=direction,tolerance=tolerance,unit.angle=unit.angle)$vario
   if(pos.neg){
    if(median)
     r.vals[,,j] <- as.matrix(rv[,c(5,7:11)])
    else
     r.vals[,,j] <- as.matrix(rv[,c(5,7:8)])
   }
   else
    r.vals[,,j] <- as.matrix(rv[,4:5])
 }}
 else{ ##not using every null model, only one null
  if(RPargs[[1]]){ ##use the random pattern algo for the spatial null
   out<-RandPatPar(psp=pop,nstrata=RPargs[[2]],mtrials1=RPargs[[3]],mtrials2=RPargs[[4]],alpha=RPargs[[5]],npar=RPargs[[6]])
   rpop<-array(0,dim=c(S,n2,n2))
   for(i in 1:S){
    rpop[i,,]<-array(out[-(1:5),i],dim=c(n+2,n+2))
   }
   rpop<-rpop[,-c(1,n+2),-c(1,n+2)]
  } 
  else{
   rpop<-SpatPerm2D(pop,meth=meth,sp=sp)
   rpop<-FixUnSamp2(pop,rpop)
  }
  rmat<-apply(rpop,1,as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
  rv<-vario(x=rmat,coord=coords,grain=grain,hmax=hmax,pos.neg=pos.neg,median=median,
            direction=direction,tolerance=tolerance,unit.angle=unit.angle)$vario
  if(pos.neg){
   if(median)
    r.vals <- as.matrix(rv[,c(5,7:11)])
   else
    r.vals <- as.matrix(rv[,c(5,7:8)])
  }
  else
   r.vals <- as.matrix(rv[,4:5])
 }
 r.vals
}

##3.5##
v.graph<-function(vobject,optim=NA,exp.only=FALSE,flip.neg=FALSE,ylim=NULL,xlab=NULL,ylab=NULL,cls=NULL){
 ##Purpose: to graph community variograms, the results of function 'vario'
 ##Arguments
 ##'vobject' the output of the function 'vario' or the function 'null.perms'
 ##'optim' is the location of the species niches
 ##'exp.only' if TRUE then only the two expected (spatial and nonspatial) components of variance displayed
 ##'flip.neg' if TRUE then negative fraction is expressed as a positive value
 n<-sqrt(vobject$parms$N)
 N<-n^2
 S<-vobject$parms$S
 coords<-cbind(rep(1:n,each=n),rep(1:n,times=n))
 v<-vobject$vario
 p<-vobject$p
 if(is.null(cls)){
  ##green, purple, orange with transparency
  cls<-c(rgb(127, 201, 127,alpha=255*.5,maxColorValue = 255),
        rgb(190, 174, 212,alpha=255*.5,maxColorValue = 255), 
        rgb(253, 192, 134,alpha=255*.5,maxColorValue = 255))
 }
 if(is.null(xlab))
  xlab<-'Spatial Lag'
 if(is.null(ylab))
  ylab<-'Variance'
 if(!is.na(optim[1])){
  y<-1:S
  plot(sort(optim),y,axes=F,frame.plot=F,xlab='grad',ylab='',xlim=c(0,max(optim)))
  axis(side=1)
  axis(side=2,at=c(1,S))
  arrows(sort(optim)-vobject$parms$niche.wid.rel/2,y,sort(optim)+vobject$parms$niche.wid.rel/2,y,angle=90,len=.01,code=3)
 }
 if(vobject$perm){ ##if vobject is result of permutations
  d<-vobject$vdists ##assign distances
  ##for each distance class I want to calculate the upper 95 and lower 95
  ##for both the neg and pos curves
  quants<-apply(v,c(1,2),quantile,c(.025,.975))
  if(is.null(ylim))
   ylim=range(list(quants,v[,,1]))
  if(vobject$parms$pos.neg){
   plot(d,v[,1,1],col=1,pch=19,type='o',ylab=ylab,xlab=xlab,ylim=ylim)
   points(d,v[,1,1]-v[,2,1]-v[,3,1],col='green3',type='o',pch=19)
   lines(d,v[,2,1],col='purple',lwd=2)
   lines(d,v[,3,1],col='orange',lwd=2)
  }
  else{
   plot(d,v[,1,1],pch=19,type='o',col='green3',ylab=ylab,xlab=xlab,ylim=ylim) #expected values
   lines(d,v[,2,1],col=1,lwd=2,type='o',pch=19) ##observed variance
  }
  lines(d,rep(sum(p*(1-p)),length(d)),col='blue',lwd=2)
  for(i in 1:dim(quants)[3])
   polygon(c(d,d[length(d):1]),c(quants[1,,i],quants[2,length(d):1,i]),col=cls[i],border = NA)
 }
 else{
  if(is.null(ylim)){
   if(vobject$parms$pos.neg){
    if(flip.neg) v[,8]<-v[,8]*-1
    ylim=range(list(v[,c(4:5,7:8)],sum(p*(1-p))))
   }
   else{
    if(exp.only)
     ylim=range(list(v[,4],sum(p*(1-p))))
    else
     ylim=range(list(v[,4:5],sum(p*(1-p))))
  }}
  plot(v$Dist,v$obs.var,col=1,pch=19,type='n',ylab=ylab,xlab=xlab,ylim=ylim)
  points(v$Dist,v$exp.var,col='green3',type='o',pch=19)
  if(!exp.only)
   points(v$Dist,v$obs.var,col=1,pch=19,type='o')
  if(vobject$parms$pos.neg){
    lines(v$Dist,v$pos,col='purple',lwd=2)
    lines(v$Dist,v$neg,col='orange',lwd=2)
  }
  lines(v$Dist,rep(sum(p*(1-p)),length(v$Dist)),col='blue',lwd=2)
 }
}



##3.6##
v.graph.all<-function(vrand=NULL,vspat=NULL,obs.var=FALSE,ylims=NA,xlims=NA,ylab='variance',xlab='lag',cls=NA,lwd=1,plot.new=TRUE){
 ##Purpose: to graph community variograms, the results of function 'vario'
 ##when all permutations have been run
 ##Arguments
 ##'vobject' the output of the function 'vario' or the function 'null.perms'
 ##'draw' may be 'exp','total', 'both' or 'pos-neg'
 rflag<-!is.null(vrand)
 sflag<-!is.null(vspat)
 if( rflag ){
  dr<-vrand$vdists ##assign distances
  vr<-vrand$vario[,1,] ##just the expected variance component
 }
 if( sflag){
  ds<-vspat$vdists
  vs<-vspat$vario
 }
 if(is.na(xlims[1]))
  xlims<-range(dr,ds)
 if(is.na(cls)){
  ##purple, blue, green
  cls<-c(rgb(190, 174, 212,alpha=255*.5,maxColorValue = 255),
         '#99CCFF',
         rgb(127, 201, 127,alpha=255*.5,maxColorValue = 255))
 }
 if(rflag & sflag & plot.new){
  if(obs.var) 
   par(mfrow=c(1,3))
  else
   par(mfrow=c(1,2))
 }
 if(rflag){
  q.rand<-apply(vr,1,quantile,c(.025,.975))
  ylims=range(q.rand,vr[,1])
  plot(dr,vr[,1],type='n',ylab=ylab,xlab=xlab,ylim=ylims,xlim=xlims,main='Within-species Agg.')
  polygon(c(dr,dr[length(dr):1]),c(q.rand[1,],q.rand[2,length(dr):1]),border=NA,col=cls[1])
  lines(dr,vr[,1],col=1,lwd=lwd)
 }
 if(sflag){
  if(obs.var){
   ovar<-vs[,2,] + vs[,3,] ##pos + neg fractions
   q.obs<-apply(ovar,1,quantile,c(.025,.975))
   ylims=range(q.obs,ovar)
   plot(ds,ovar[,1],type='n',ylab=ylab,xlab=xlab,ylim=ylims,xlim=xlims,main='Total Between-species Agg.')
   polygon(c(ds,ds[length(ds):1]),c(q.obs[1,],q.obs[2,length(ds):1]),border=NA,col=cls[1])
   lines(ds,ovar[,1],col=1,lwd=lwd)
  }
  q.spat<-apply(vs,1:2,quantile,c(.025,.975))
  ylims=range(q.spat,vs[,,1])
  plot(ds,vs[,2,1],type='n',ylab=ylab,xlab=xlab,ylim=ylims,xlim=xlims,main='Pos/Neg Between-species Agg.')
  polygon(c(ds,ds[length(ds):1]),c(q.spat[1,,2],q.spat[2,length(ds):1,2]),border=NA,col=cls[2])
  polygon(c(ds,ds[length(ds):1]),c(q.spat[1,,3],q.spat[2,length(ds):1,3]),border=NA,col=cls[3])
  lines(ds,vs[,2,1],col='dodgerblue',lwd=lwd)
  lines(ds,vs[,3,1],col='green3',lwd=lwd)
}}


##3.7## This version has additional options for setting ylim and xlim that do not fully work
v.graph.all2<-function(vrand=NULL,vspat=NULL,obs.var=FALSE,flip.neg=FALSE,
                      ylim1=NA,ylim2=NA,xlim=NA,ylab='variance',xlab='lag',
                      cls=NA,lwd=1,cex.axis=1,box=TRUE,plot.new=TRUE){
 ##Purpose: to graph community variograms, the results of function 'vario'
 ##when all permutations have been run
 ##Arguments
 ##'vrand': the output from null.perms for the random permutations
 ##'vspat' the output from null.perms for the spatial permutations
 ##obs.var: if TRUE then this fraction is ploted
 ##flip.neg: if TRUE then neg fraction of covariance is plotted as a positive value
 rflag<-!is.null(vrand)
 sflag<-!is.null(vspat)
 dr<-NULL;ds<-NULL
 if( rflag ){
  dr<-vrand$vdists ##assign distances
  vr<-vrand$vario[,1,] ##just the expected variance component
 }
 if( sflag){
  ds<-vspat$vdists
  vs<-vspat$vario
  if(flip.neg) vs[,3,]<-vs[,3,]*-1
 }
 if(is.na(xlim[1]))
  xlim<-range(dr,ds)
 if(is.na(cls)){
  ##purple, blue, green
  cls<-c(rgb(190, 174, 212,alpha=255*.5,maxColorValue = 255),
         '#99CCFF',
         rgb(127, 201, 127,alpha=255*.5,maxColorValue = 255))
 }
 if(rflag & sflag & plot.new){
  if(obs.var) 
   par(mfrow=c(1,3))
  else
   par(mfrow=c(1,2))
 }
 if(rflag){
  q.rand<-apply(vr,1,quantile,c(.025,.975))
  if(is.na(ylim1[1])) ylim1=range(q.rand,vr[,1])
  plot(dr,vr[,1],type='n',ylab=ylab,xlab=xlab,ylim=ylim1,xlim=xlim,main='Within-species Agg.',axes=F)
  axis(side=1,cex.axis=cex.axis)
  axis(side=2,cex.axis=cex.axis)
  polygon(c(dr,dr[length(dr):1]),c(q.rand[1,],q.rand[2,length(dr):1]),border=NA,col=cls[1])
  lines(dr,vr[,1],col=1,lwd=lwd)
  if(box)box()
 }
 if(sflag){
  if(obs.var){
   if(flip.neg)  ovar<-vs[,2,] - vs[,3,] ##pos + neg fractions
   else ovar<-vs[,2,] + vs[,3,] ##pos + neg fractions
   q.obs<-apply(ovar,1,quantile,c(.025,.975))
   if(is.na(ylim1[1])) ylim1=range(q.obs,ovar)
   plot(ds,ovar[,1],type='n',ylab=ylab,xlab=xlab,ylim=ylim1,xlim=xlim,main='Total Between-species Agg.',axes=F)
   axis(side=1,cex.axis=cex.axis)
   axis(side=2,cex.axis=cex.axis)
   polygon(c(ds,ds[length(ds):1]),c(q.obs[1,],q.obs[2,length(ds):1]),border=NA,col=cls[1])
   lines(ds,ovar[,1],col=1,lwd=lwd)
   if(box)box()
  }
  q.spat<-apply(vs,1:2,quantile,c(.025,.975))
  if(is.na(ylim2[1])) ylim2=range(q.spat,vs[,,1])
  plot(ds,vs[,2,1],type='n',ylab=ylab,xlab=xlab,ylim=ylim2,xlim=xlim,main='Pos/Neg Between-species Agg.',axes=F)
  axis(side=1,cex.axis=cex.axis)
  axis(side=2,cex.axis=cex.axis)
  polygon(c(ds,ds[length(ds):1]),c(q.spat[1,,2],q.spat[2,length(ds):1,2]),border=NA,col=cls[2])
  polygon(c(ds,ds[length(ds):1]),c(q.spat[1,,3],q.spat[2,length(ds):1,3]),border=NA,col=cls[3])
  lines(ds,vs[,2,1],col='dodgerblue',lwd=lwd)
  lines(ds,vs[,3,1],col='green3',lwd=lwd)
  if(box)box()
}}

'mat2psp' = function(x,N=NULL,M=NULL)
{
  ##place site by species matrix (x) into an S x N x M array where N >= M
  ##a multidimensional array that eases computing 
  ##replaces the old function 'grid.pres'
  S = ncol(x)
  if(is.null(N))
    N = sqrt(nrow(x))
  if(is.null(M))
    M = N
  psp = array(x,dim=c(N,M,S))
  psp = aperm(psp,c(3,1,2))
  spSums = apply(psp,1,sum)
  ## drop species that never occur
  if(any(spSums %in% 0))
    psp = psp[spSums > 0,,]
  return(psp)
}

##3.8##
'getSAR' = function(psp, grains, mv_window=FALSE)
{
  ## Purpose: to construct spatially explict SAR based upon a
  ## mapped grid of occurances
  ## this function replaces the older function 'grid.SAR'
  ## Arguments:
  ## psp: community array (i.e., S x N x M pres/absen array where N >= M)
  ## grains: the areas in pixels for which to compute the SAR
  ##         only grains that have integer log base 2 are considered
  ## mv_window: FALSE indicates that a non-moving window SAR will be calculated
  ## Note: This implementation may require that the grain of 
  if (class(psp) != 'array')
    stop('psp must be a community array (S X N X M)')
  grains = grains[log2(grains) == round(log2(grains))]
  grainsSqr = grains[sqrt(grains) == round(sqrt(grains))]
  ## define the size of sampling units on each side
  lenN = 2^ceiling(log2(grains) / 2)
  lenM = 2^floor(log2(grains) / 2)
  sr = rep(0, length(grains))
  ind = rep(0, length(grains))
  cs = rep(0, length(grains))
  N = dim(psp)[2]
  M = dim(psp)[3]
  if (M > N) {
    stop('The first spatial dimension of psp must be larger than or equal to the second 
         (i.e. psp[S,N,M] where N >= M)')
  }       
  for (l in seq_along(grains)) {
    if (grains[l] == 1) {  # if area=1
      sr[l] = sum(psp > 0)
      ind[l] = sum(psp)
      cs[l] = N * M
    }
    else{
      if (mv_window) {
        brksN = 1:(N - lenN[l] + 1)
        brksM = 1:(M - lenM[l] + 1)
      }
      else{
        brksN = seq(1, N, lenN[l])
        brksM = seq(1, M, lenM[l])
      }  
      for (n in brksN) {
        for (m in brksM) {
          sr[l] = sr[l] + sum(apply(psp[, n:(n + (lenN[l] - 1)),
                                          m:(m + (lenM[l] - 1))] > 0, 1, sum) > 0)
          ind[l] = ind[l] + sum(apply(psp[, n:(n + (lenN[l] - 1)),
                                            m:(m + (lenM[l] - 1))], 1, sum))
          cs[l] = cs[l] + 1
        }
      }
    }
  }
  out = cbind(grains, sr / cs, ind / cs, cs)  
  colnames(out) = c('grains', 'richness', 'indiv', 'count')
  return(out)
}  


##3.9##
'quadAggregator' = function(mat, coords, grains, binary=FALSE){
  ## Purpose: This function generates aggregated community matrices for each
  ## spatial grain that is specified. This function is only approrpriate for 
  ## data from a regular square spatial grid.  If coordinates are not supplied
  ## they will be generated as will spatial grains. 
  ## Inputs:
  ## mat: site x species matrix assumed to be sampled from a square grid
  ## coords: two column spatial x and y coordinates
  ## grains: spatial grains to aggregate community at, thse  must represent 
  ##         quadrasections (i.e. aggregations of groups of 4 quadrats)
  ## binary: boolean, if TRUE binary pres/abse matrix returned
  ## Note: assumes that spatially defined square was sampled.
  ## ToDo: 1) Allow rectrangular quadrats, 2) allow bisection aggregation
  side_length = sqrt(nrow(mat)) 
  if (missing(coords)) {
    coords = expand.grid(1:side_length, 1:side_length)
  }
  if (missing(grains)) {
    grains = (2^(0:(log2(side_length) - 1)))^2
  }
  grains = grains[sqrt(grains) == round(sqrt(grains))]
  lens = sqrt(grains) 
  nRows = sum((side_length / lens)^2)
  irow = 1
  out = matrix(NA, ncol=ncol(mat) + 3, nrow=nRows)
  for (i in seq_along(grains)) {
    if (lens[i] == 1) {  # if area == 1
      out[1:nrow(mat),] = as.matrix(cbind(grains[i], coords, mat))
      irow = nrow(mat) + 1
    }
    else{
      breaks = seq(1, side_length, lens[i])
      for(x in breaks){
        for(y in breaks){
          true = coords[,1] >= x & coords[,1] < (x + lens[i]) &
                 coords[,2] >= y & coords[,2] < (y + lens[i])
          out[irow, ] = c(grains[i], 
                          apply(coords[true,], 2, mean), 
                          apply(mat[true,], 2, sum))
          irow = irow +1
        }
      }
    }
  }  
  colnames(out) = c('grain', 'x', 'y', colnames(mat))
  if (binary) {
    out[,-(1:3)] = as.numeric(out[,-(1:3)] > 0)
  }
  return(out)
}  


##3.10##
varExp = function(mat){
  #first evaulate if pres/abs or abundance matrix
  spmeans = apply(mat,2,mean)
  if(sum(mat > 1)) #  for abundance
   sum(spmeans) 
  else             #  for occupancy
   sum(spmeans * (1 - spmeans)) 
}

##3.11##
spCommExpPoi = function(a,b,n){
   sum(apply(sapply(n,function(x) sapply(c(a,b),function(y) 1 - exp(-x*y))),2,prod))
}

##3.12##
spAvgExpPoi = function(a,b,n){
   .5 * sum(sapply(n,function(x) sapply(c(a,b),function(y) 1 - exp(-x*y))))
}


##3.13##
spCommExpBin = function(a,b,areaTot,numOcc){
  sum(apply(sapply(numOcc,function(x) sapply(c(a,b),function(y) 1 - (1 - (x/areaTot))^(y*areaTot))),2,prod))
}

##3.14##
spAvgExpBin= function(a,b,areaTot,numOcc){
   .5 * sum(sapply(numOcc,function(x) sapply(c(a,b),function(y) 1 - (1 - (x/areaTot))^(y*areaTot))))
}

##3.15##
sorExp = function(mat,areaSampA,areaSampB=NULL){
  if(is.null(areaSampB)) 
    areaSampB = areaSampA
  areaTot = nrow(mat)  
  a = areaSampA / areaTot
  b = areaSampB / areaTot
  #first evaulate if pres/abs or abundance matrix
  if(sum(mat > 1)){ #  for abundance use Poisson expectation
    n = apply(mat,2,sum)
    sor = spCommExpPoi(a,b,n) / spAvgExpPoi(a,b,n)
  }
  else{
    numOcc = apply(mat,2,sum)
    sor = spCommExpBin(a,b,areaTot,numOcc) / spAvgExpBin(a,b,areaTot,numOcc)
  }
 return(sor)
}

##3.16##
jacExp = function(mat,areaSampA,areaSampB=NULL){
  if(is.null(areaSampB)) 
    areaSampB = areaSampA
  areaTot = nrow(mat)  
  a = areaSampA / areaTot
  b = areaSampB / areaTot
  #first evaulate if pres/abs or abundance matrix
  if(sum(mat > 1)){ #  for abundance use Poisson expectation
    n = apply(mat,2,sum)
    jac = spCommExpPoi(a,b,n) / (2*spAvgExpPoi(a,b,n) - spCommExpPoi(a,b,n))
  }
  else{
    numOcc = apply(mat,2,sum)
    areaTot = nrow(mat)
    jac = spCommExpBin(a,b,areaTot,numOcc) / (2*spAvgExpBin(a,b,areaTot,numOcc) - spCommExpBin(a,b,areaTot,numOcc))

  }
 return(jac)
}

##3.17##
calcMetrics = function(comms,metricsToCalc,dataType,grain=1,breaks=NA,hmin=NA,
                       hmax=NA,quants=NA,
                       direction='omnidirectional',tolerance=NA,nperm=NULL,npar,
                       RPargs=NULL,writeToFile=FALSE,fileSuffix=NULL){
  ## Purpose: to compuate spatial distance decay metrics for community data.
  ## Metrics to choose from are varWithin,varBetween, jaccard, and sorensen
  ## indices.  
  ## Arguments: 
  ## comms: a site x sp community matrix the first column specifies the
  ##        community id, the 2nd and 3rd columns specify the spatial
  ##        coordinates of the quadrat
  ## metricsToCalc: a string specifing which metrics to calc, can be 'all'
  ## dataType: if == 'binary' then comms is converted to a pres/absence matrix
  ##           prior to analysis. If == 'abu' then matrix is not transformed
  ##           and an additional analytical null expectation is calculated
  ## grain: interval size for distance classes, only used if 'breaks' not supplied
  ## breaks: the spatial breaks that define the spatial lags to be compared
  ## hmin: the minimum spatial lag
  ## hmax: the maximum spatial lag
  ## quants: the quantiles to compute
  ## nperm: number of permutations to carry out for null models
  ## npar: number of processesors to run null models on
  ## RPargs: arguments to parameterize the Random Patterns null model
  ## writeToFile: if True an .Rdata file is written for each metric calculated
  ## fileSuffix: add a file identifying string here
  if(metricsToCalc == 'all')
    metricsToCalc = c('varWithin','varBetween','jaccard','sorensen')
  if(writeToFile){
    if(direction != 'omnidirectional')
      fileSuffix = paste(fileSuffix,'_',direction,'deg',sep='') 
  }
  grains = unique(comms[,1])
  out = vector('list',length(grains))
  names(out) = paste('comm',grains,sep='')
  for(i in seq_along(grains)){
    coords = as.matrix(comms[comms[,1] == grains[i], 2:3]) * grains[i]
    mat = as.matrix(comms[comms[,1] == grains[i],-c(1:3)])
    if (!is.na(breaks)) {
      if (is.list(breaks))
        brks = breaks[[i]]
      else
        brks = breaks
    }  
    if(dataType == 'binary')
      mat = (mat > 0) * 1
    if(any('varWithin' %in% metricsToCalc)){
      if(i == 1){
        varWithin = vector('list', length(grains))      
        names(varWithin) = grains
      }  
      varWithinObs = vario(mat,coords,grain,brks,hmin,hmax,pos.neg=FALSE,
                           quants=quants,direction=direction,tolerance=tolerance,
                           unit.angle='degrees')
      if(!is.null(nperm)){ 
        varWithinNull = null.perms(mat,varWithinObs,nperm,coords=coords,
                                   meth='random',npar=npar)
      }
      else{
        varWithinNull = NULL
      }       
      varWithinExp = varExp(mat)
      varWithin[[i]] = list(varWithinObs=varWithinObs,varWithinNull=varWithinNull,
                            varWithinExp=varWithinExp)
    }
    if(any('varBetween' %in% metricsToCalc)){
      if(i == 1){
        varBetween = vector('list', length(grains))      
        names(varBetween) = grains
      }  
      varBetweenObs = vario(mat,coords,grain,brks,hmin,hmax,pos.neg=TRUE,
                            quants=quants,direction=direction,tolerance=tolerance,
                            unit.angle='degrees') 
      if(!is.null(nperm)){ 
        varBetweenNull = null.perms(mat,varBetweenObs,nperm,coords=coords,
                                    meth='randpat',RPargs=RPargs,npar=npar)
      }
      else{
        varBetweenNull = NULL
      }   
      varBetween[[i]] = list(varBetweenObs=varBetweenObs,
                             varBetweenNull=varBetweenNull)
    }
    if(any('jaccard' %in% metricsToCalc)){
      if(i == 1){
        jaccard = vector('list', length(grains))  
        names(jaccard) = grains
      }  
      jaccardObs  = vario(mat,coords,grain,brks,hmin,hmax,distance.metric='jaccard',
                          quants=quants, direction=direction,tolerance=tolerance,
                           unit.angle='degrees') 
      jaccardNull = NULL
      jaccardExp = 1 - jacExp(mat,1) #  to convert into a dissimiarlity
      if(dataType == 'abu'){
        jaccardExpAbuToBinary = 1 - jacExp(mat > 0,1)
      }
      else{
        jaccardExpAbuToBinary = NULL
      } 
      jaccard[[i]] = list(jaccardObs=jaccardObs,jaccardNull=jaccardNull,
                          jaccardExp=jaccardExp,
                          jaccardExpAbuToBinary=jaccardExpAbuToBinary)
    }
    if(any('sorensen' %in% metricsToCalc)){
      if(i == 1){
        sorensen = vector('list', length(grains))
        names(sorensen) = grains
      }  
      ## bray-curtis is equiv to sorensen        
      sorensenObs  = vario(mat,coords,grain,brks,hmin,hmax,distance.metric='bray',
                           quants=quants,direction=direction,tolerance=tolerance,
                           unit.angle='degrees') 
      sorensenNull = NULL
      sorensenExp = 1 - sorExp(mat,1) #  to convert into a dissimiarlity
      if(dataType == 'abu'){
        sorensenExpAbuToBinary = 1 - sorExp(mat > 0,1)
      }
      else{
        sorensenExpAbuToBinary = NULL
      } 
      sorensen[[i]] = list(sorensenObs=sorensenObs, sorensenNull=sorensenNull,
                           sorensenExp=sorensenExp,
                           sorensenExpAbuToBinary=sorensenExpAbuToBinary)
    }
    out[[i]] = list()
    for(j in metricsToCalc){ 
      out[[i]][[j]] = eval(parse(text=paste(j,'[[',i,']]')))
    }
    if(writeToFile){      
      ## update result files as loop proceeds
      if(any('varWithin' %in% metricsToCalc)){
        save(varWithin,file=paste(getwd(),'/varWithin/varWithin_',
             fileSuffix,'_',dataType,'.Rdata',sep=''))  
      }
      if(any('varBetween' %in% metricsToCalc)){
        save(varBetween,file=paste(getwd(),'/varBetween/varBetween_',
             fileSuffix,'_',dataType,'.Rdata',sep=''))
      }
      if(any('jaccard' %in% metricsToCalc)){
        save(jaccard,file=paste(getwd(),'/jaccard/jaccard_',
             fileSuffix,'_',dataType,'.Rdata',sep=''))
      }
      if(any('sorensen' %in% metricsToCalc)){
        save(sorensen,file=paste(getwd(),'/sorensen/sorensen_',
             fileSuffix,'_',dataType,'.Rdata',sep=''))
      }             
    }
  }
  return(out)
}


calcMetricsPar = function(comms,metricsToCalc,dataType,npar,grain=1,breaks=NA,
                          hmax=NA,quants=NA,direction='omidirectional',
                          tolerance=NA,nperm=NULL,RPargs=NULL,writeToFile=FALSE,
                          fileSuffix=NULL){
  ## Purpose: to compuate spatial distance decay metrics for community data.
  ## Metrics to choose from are varWithin,varBetween, jaccard, and sorensen
  ## indices.  
  ## Arguments: 
  ## comms: a site x sp community matrix the first column specifies the
  ##        community id, the 2nd and 3rd columns specify the spatial
  ##        coordinates of the quadrat
  ## metricsToCalc: a string specifing which metrics to calc, can be 'all'
  ## dataType: if == 'binary' then comms is converted to a pres/absence matrix
  ##           prior to analysis. If == 'abu' then matrix is not transformed
  ##           and an additional analytical null expectation is calculated
  ## npar: number of processers to run the community for loop across  
  ## grain: interval size for distance classes, only used if 'breaks' not supplied
  ## breaks: the spatial breaks that define the spatial lags to be compared
  ## hmax: the maximum spatial lag
  ## nperm: number of permutations to carry out for null models
  ## writeToFile: if True an .Rdata file is written for each metric calculated
  ## fileSuffix: add a file identifying string here
  require(foreach)
  require(doSNOW) ##needed to initialize the cluster with foreach
  require(snowfall) ##needed to create cluster
  sfInit(parallel=TRUE, cpus=npar, type="SOCK")
  registerDoSNOW(sfGetCluster())
  if(metricsToCalc == 'all')
    metricsToCalc = c('varWithin','varBetween','jaccard','sorensen')
  if(writeToFile){
    if(direction != 'omnidirectional')
      fileSuffix = paste(fileSuffix,'_',direction,'deg',sep='') 
  }
  commNames = unique(comms[,1])
  out = vector('list',length(commNames))
  names(out) = paste('comm',commNames,sep='')
  foreach(i = 1:length(commNames), .inorder = TRUE) %dopar% {
    coords = as.matrix(comms[comms[,1] == commNames[i],2:3])
    mat = as.matrix(comms[comms[,1] == commNames[i],-c(1:3)])
    if(dataType == 'binary')
      mat = (mat > 0) * 1
    if(any('varWithin' %in% metricsToCalc)){
      if(i == 1){
        varWithin = vector('list', length(commNames))      
        names(varWithin) = commNames
      }     
      varWithinObs = vario(mat,coords,grain,breaks,hmax,pos.neg=FALSE,
                           quants=quants,direction=direction,tolerance=tolerance,
                           unit.angle='degrees')
      if(!is.null(nperm)){ 
        varWithinNull = null.perms(mat,varWithinObs,nperm,coords=coords,
                                   meth='random')
      }
      else{
        varWithinNull = NULL
      }       
      varWithinExp = varExp(mat)
      varWithin[[i]] = list(varWithinObs=varWithinObs,varWithinNull=varWithinNull,
                            varWithinExp=varWithinExp)
    }
    if(any('varBetween' %in% metricsToCalc)){
      if(i == 1){
        varBetween = vector('list', length(commNames))      
        names(varBetween) = commNames
      }       
      varBetweenObs = vario(mat,coords,grain,breaks,hmax,pos.neg=TRUE,quants=quants,
                           direction=direction,tolerance=tolerance,
                           unit.angle='degrees') 
      if(!is.null(nperm)){ 
        varBetweenNull = null.perms(mat,varBetweenObs,nperm,coords=coords,
                                    meth='randpat',RPargs=RPargs)
      }
      else{
        varBetweenNull = NULL
      }   
      varBetween[[i]] = list(varBetweenObs=varBetweenObs,
                             varBetweenNull=varBetweenNull)
    }
    if(any('jaccard' %in% metricsToCalc)){
      if(i == 1){
        jaccard = vector('list', length(commNames))  
        names(jaccard) = commNames
      }  
      jaccardObs  = vario(mat,coords,grain,breaks,hmax,distance.metric='jaccard',
                           quants=quants, direction=direction,tolerance=tolerance,
                           unit.angle='degrees') 
      jaccardNull = NULL
      jaccardExp = 1 - jacExp(mat,1) #  to convert into a dissimiarlity
      if(dataType == 'abu'){
        jaccardExpAbuToBinary = 1 - jacExp(mat > 0,1)
      }
      else{
        jaccardExpAbuToBinary = NULL
      } 
      jaccard[[i]] = list(jaccardObs=jaccardObs,jaccardNull=jaccardNull,
                          jaccardExp=jaccardExp,
                          jaccardExpAbuToBinary=jaccardExpAbuToBinary)
    }
    if(any('sorensen' %in% metricsToCalc)){
      if(i == 1){
        sorensen = vector('list', length(commNames))
        names(sorensen) = commNames
      }   
      ## bray-curtis is equiv to sorensen        
      sorensenObs  = vario(mat,coords,grain,breaks,hmax,distance.metric='bray',
                           quants=quants,direction=direction,tolerance=tolerance,
                           unit.angle='degrees') 
      sorensenNull = NULL
      sorensenExp = 1 - sorExp(mat,1) #  to convert into a dissimiarlity
      if(dataType == 'abu'){
        sorensenExpAbuToBinary = 1 - sorExp(mat > 0,1)
      }
      else{
        sorensenExpAbuToBinary = NULL
      } 
      sorensen[[i]] = list(sorensenObs=sorensenObs, sorensenNull=sorensenNull,
                           sorensenExp=sorensenExp,
                           sorensenExpAbuToBinary=sorensenExpAbuToBinary)
    }
    out[[i]] = list()
    for(j in metricsToCalc){ 
      out[[i]][[j]] = eval(parse(text=paste(j,'[[',i,']]')))
    }
    if(writeToFile){      
      ## update result files as loop proceeds
      if(any('varWithin' %in% metricsToCalc)){
        save(varWithin,file=paste(getwd(),'/varWithin/varWithin',
             fileSuffix,'_',dataType,'.Rdata',sep=''))  
      }
      if(any('varBetween' %in% metricsToCalc)){
        save(varBetween,file=paste(getwd(),'/varBetween/varBetween',
             fileSuffix,'_',dataType,'.Rdata',sep=''))
      }
      if(any('jaccard' %in% metricsToCalc)){
        save(jaccard,file=paste(getwd(),'/jaccard/jaccard',
             fileSuffix,'_',dataType,'.Rdata',sep=''))
      }
      if(any('sorensen' %in% metricsToCalc)){
        save(sorensen,file=paste(getwd(),'/sorensen/sorensen',
             fileSuffix,'_',dataType,'.Rdata',sep=''))
      }             
    }
  }
  return(out)
}


n_pixels_long = function(i_bisections){
  ## returns the number of pixels on the side of a grid with more or equal pixels
  ## after i bisection events
  ## old function name: len
  2^floor((i_bisections + 1) / 2) 
}

n_pixels_wide = function(i_bisections){
  ## returns the number of pixels on the side of a grid with less or equal pixels
  ## after i bisection events
  ## old function name: wid
  2^floor(i_bisections / 2) 
} 

##3.19 ##
'make_comm_matrix' = function(spnum, S, coords, n_quadrats, domain, abu = NULL,
                              grainSuffix=NULL)
{ 
  ## Output: 
  ## A community matrix where each row is a differnet pixel on a grid.  
  ## Arguments:
  ## spnum : an integer specifying species identities
  ## S : the size of the species pool may be larger than the number of unique 
  ##     spnum
  ## coords : two column matrix (x,y) specifying the spatial coordinates
  ## n_quadrats : the number of quadrats at each spatial grain
  ## domain : specifies the spatial domain of the area:  (xmin, xmax, ymin, ymax)
  ## abu: abundance associated with each record, if NULL then it is set to 1
  ##      individual per record
  ## grainSuffix : if supplied the grain column will have this appended to it
  ##               so that it is clear what community this corresponds with
  xdiff = abs(domain[2] - domain[1])
  ydiff = abs(domain[4] - domain[3])
  if (xdiff > ydiff) {
    xlengths = xdiff / n_pixels_long(log2(n_quadrats))
    ylengths = ydiff / n_pixels_wide(log2(n_quadrats))
  }  
  else if (xdiff < ydiff) {
    xlengths = xdiff / n_pixels_wide(log2(n_quadrats))
    ylengths = ydiff / n_pixels_long(log2(n_quadrats))
  }
  else if (xdiff == ydiff) {
    xlengths = xdiff / sqrt(n_quadrats)
    ylengths = ydiff / sqrt(n_quadrats)
  }
  else
    stop('Function cannot figure out how to split up the area')
  comms = matrix(NA, nrow=sum(n_quadrats), ncol=S + 3)
  colnames(comms) = c('grain', 'x', 'y', paste0('sp', 1:S))
  irow = 1
  for (i in seq_along(n_quadrats)) {
    xbreaks = seq(domain[1], domain[2], xlengths[i])
    ybreaks = seq(domain[3], domain[4], ylengths[i]) 
    for (x in 1:(length(xbreaks) - 1)) {
      for (y in 1:(length(ybreaks) - 1)) {
        inQuad =  xbreaks[x] <= coords[,1] & coords[,1] < xbreaks[x + 1] & 
                  ybreaks[y] <= coords[,2] & coords[,2] < ybreaks[y + 1]
        if (is.null(grainSuffix)) {
          comms[irow, c(1:3)] = c(paste0(round(xlengths[i] * ylengths[i], 2)),
                                  x, y)
        }
        else {
          comms[irow, c(1:3)] = c(paste0(round(xlengths[i] * ylengths[i], 2),
                                  grainSuffix), x, y)
        }
        if (is.null(abu) ){
          comms[irow, -c(1:3)] = as.integer(table(c(spnum[inQuad],1:S)) - 1)
        }
        else {
          comms[irow, -c(1:3)] =  as.integer(table(c(unlist(mapply(
                                     rep, spnum[inQuad], abu[inQuad])), 1:S)) - 1)
        }  
        irow = irow + 1 
      }
    }
  }
  return(comms)
}

'get_bisection_history' = function(grains, abu) {
  ## retrives the record of how a species' abunance was bisected in a particular
  ## landscape the input into this function can be generated by the function 
  ## 'make_comm_matrix'
  ## arguments:
  ## grains: a vector of spatial grains that is associated with each abundance
  ## abu: a vector of abundances across the spatial scales of interest
  abu = as.numeric(abu)
  grains = as.numeric(grains)
  grains_table = table(grains)
  grains_unique = as.integer(names(grains_table))
  i_bisections = log2(as.integer(grains_table))
  out = data.frame(n1 = NULL, n2 = NULL)
  for (i in seq_along(grains_unique)) {
    abu_tmp = abu[grains == grains_unique[i]]
    if (i_bisections[i] == 1) {
      n1 = abu_tmp[1] 
      n2 = abu_tmp[2]
    }  
    else {
      if (i_bisections[i]%%2 == 0) {
        true = 1:length(abu_tmp)%%2 == 1
        n1 = abu_tmp[true]
        n2 = abu_tmp[!true]          
      }
      else {
        block_size = 2^(i_bisections[i] - ceiling(i_bisections[i] / 2))
        true = rep(rep(c(T,F), each=block_size), times=block_size)
        n1 = abu_tmp[true]
        n2 = abu_tmp[!true]
      }
    }
    out = rbind(out,data.frame(n1 = n1, n2 = n2))
  }              
  return(out)
}

'getResults' = function(names,metric,dataType)
{
  ## Purpose: to import the results of the 'calcMetrics' function
  ## and to load them into a list
  ## Arguments:
  ## names: the short names that were used in the naming of the output files
  ## metric: the metric calculated, see 'calcMetrics' for options
  ## dataType: the data type of interest, 'binary' or 'abu'
  results = vector('list',length=length(names))
  names(results) = names
  for(i in seq_along(results)){
    load(paste('./',metric,'/',metric,'_',names[i],'_',dataType,'.Rdata',
               sep=''))
    results[[i]] = eval(parse(text=metric))
  }
  return(results)
}

'reshapeResults' = function(results, metric){
  ## Purpose: to reshape the results from a nested list to a matrix
  ## of the most important information
  out = vector('list', length=length(results))
  names(out) = names(results)
  for (i in seq_along(results)) { ## the datset
    if (is.null(results[[i]]))
      next
    vExp = vector('list',length=length(results[[i]]))
    Dist = vExp
    n = vExp
    rpExp = vExp 
    for(j in seq_along(results[[i]])){ ## the grain/community
      if(is.null(results[[i]][[j]]))
        next
      vobject = results[[i]][[j]][[1]]$vario
      quants = !is.na(results[[i]][[j]][[1]]$parms$quants)
      Dist[[j]] = vobject$Dist 
      n[[j]] = vobject$n
      rpExp[[j]] = rep(results[[i]][[j]][[3]], length(Dist[[j]]))
      if (quants) {
        qt.names = names(vobject)[grep('exp.qt',names(vobject))]
        vExp[[j]] = vobject[ , qt.names]
      }
      else 
        vExp[[j]] = vobject$exp.var
      if (any(metric %in% c('sorensen','jaccard'))) {
        vExp[[j]] = 1 - vExp[[j]]
      }  
      if (j == 1) 
        vExp_dat = vExp[[j]]
      else
        vExp_dat = rbind(vExp_dat, vExp[[j]])
    }
    if(is.null(names(results[[i]])))
      commNames = 1:length(results[[i]])
    else
      commNames = names(results[[i]])
    names(vExp_dat) = sub('exp.qt.','',names(vExp_dat))
    names(vExp_dat) = rev(names(vExp_dat))
    out[[i]] = data.frame(Dist = unlist(Dist), Metric = vExp_dat,
                          Exp = unlist(rpExp), N = unlist(n))
    out[[i]] = data.frame(out[[i]],
               Comm = as.numeric(unlist(mapply(rep,commNames,each=sapply(n, length),
                                       SIMPLIFY=FALSE))))
  }
  return(out)
}  

avgResults = function(results,combine = NULL){
  ## Purpose: to compute the averages of the results given a vector 'combine'
  ## Arguments
  ## combine:  a matrix, each row specifies an index, each column specifies a 
  ##           different one to average over. If combine is NA then no elements
  ##           are averaged over.
  out = vector('list',length(results))
  names(out) = names(results)
  for(i in seq_along(results)){
    if(is.null(results[[i]]))
      next
    if(is.na(combine[[i]][1]))
      combine[[i]] = results[[i]]$Comm
    unicombine = unique(combine[[i]])
    for(j in seq_along(unicombine)){
      true = combine[[i]] == unicombine[j]
      Dist = tapply(results[[i]]$Dist[true],round(results[[i]]$Dist[true],3),mean)
      Metric = tapply(results[[i]]$Metric[true],round(results[[i]]$Dist[true],3),mean)
      MetricLo = tapply(results[[i]]$Metric[true],round(results[[i]]$Dist[true],3),quantile,.025)
      MetricHi = tapply(results[[i]]$Metric[true],round(results[[i]]$Dist[true],3),quantile,.975)    
      Exp = tapply(results[[i]]$Exp[true],round(results[[i]]$Dist[true],3),mean)
      ExpLo = tapply(results[[i]]$Exp[true],round(results[[i]]$Dist[true],3),quantile,.025)
      ExpHi = tapply(results[[i]]$Exp[true],round(results[[i]]$Dist[true],3),quantile,.975)    
      N = tapply(results[[i]]$N[true],round(results[[i]]$Dist[true],3),sum)
      Comm = rep(unicombine[j],length(Metric))
      if(j == 1){
        out[[i]] = data.frame(Dist,Metric,MetricLo,MetricHi,Exp,ExpLo,ExpHi,N,Comm)
      }  
      else{
        out[[i]] = rbind(out[[i]],data.frame(Dist,Metric,MetricLo,MetricHi,Exp,
                                             ExpLo,ExpHi,N,Comm))
      }  
    }  
  }
  return(out)
}

plotEmpir = function(results,log="",quants=FALSE, alpha=1/3,
                     col=NULL,lwd=NULL, ...){
  ## Purpose: to plot the results, expects that the graphical window has been 
  ## setup appropriately 
  ## Arguments
  ## results: a list from which the results should be ploted
  ## log: which axes are to be log transformed
  ## quants: if TRUE will plot quantiles as well
  ## col: what colors to use
  ## lwd: the line width of the lines
  if (is.null(lwd))
    lwd = 2
  for(i in seq_along(results)){
    #grains = as.numeric(as.character(results[[i]]$Comm))
    grains = results[[i]]$Comm
    unigrains = unique(grains)
    results[[i]]$Dist = results[[i]]$Dist * sqrt(grains)
    resp_var = names(results[[i]])[grep('Metric', names(results[[i]]))]
    if (length(resp_var) > 1) {
      names(results[[i]]) = sub('.50','',names(results[[i]]))
    }   
    else {
      if (quants)
        stop("This results file does not contain data on quantiles")
    } 
    if(is.null(col))
      col = palette()[-1]
    if (quants) {
      col_poly = apply(col2rgb(col), 2, 
                  function(x) rgb(x[1],x[2],x[3], 
                                  alpha=alpha * max(x), 
                                  maxColorValue = max(x)))
    }  
    plot(Metric ~ Dist, data = results[[i]],
         xlim = range(Dist), ylim = range(Metric), type='n', log=log,
         main=names(results)[i])
    for (j in seq_along(unigrains)) {
      dat = subset(results[[i]], Comm == unigrains[j])
      if (quants) {
        xvals = dat$Dist
        xvals = c(xvals, rev(xvals))
        yvals = dat$Metric.25
        yvals = c(yvals, rev(dat$Metric.75))
        polygon(xvals, yvals, border=NA, col=col_poly[j])
        lines(Metric ~ Dist, data = dat, col=col[j], lwd=lwd, ...)
      }  
      else { 
        lines(Metric ~ Dist, data = dat, col=col[j], lwd=lwd, ...)
      }  
    }  
  }
}

#####Part IV - BATCH FUNCTIONS FOR GENERATING LARGE SETS OF RESULTS#####

##4.1##
sim.vario <- function(PREFIX, DESTDIR = paste(getwd(),"/results_main",sep=""), REPS = 10, REP.ID = NULL,
                     M = 64, ENV = 'fractal', K = 250, S = 10, NICH.WD = c(.1,.5), DISP.WD = c(.1,.5), 
                     DEATH = 0.5, BIRTH = 0.505, IMMIGRATE = 0.005, FEC = 1, FRACTAL = c(2.01, 2.25, 2.5, 2.75, 2.99),
                     CYCLES = 1, TIME = 10000, MIG = 1, NPERM = 100, ALL = TRUE, NPAR = 5,
                     RPARGS = c(TRUE,M/3,1e3,1e6,0.01,1), DESTDIRRAW = paste(getwd(),"/results_raw",sep=""),
                     SIMINPUT=NULL) { 
 ## simulation parameters:
 ## PREFIX : prefix to use to identify simulation runs
 ## DESTDIR : directory where only variogram results are returned
 ## REPS : number of reps per each combination of parameters
 ## REP.ID : the current replicate id to identify the output file as (this is only really for large serial batches)
 ## M : size of landscape along one spatial axis
 ## ENV : type of envrionment - default is fractal
 ## K : carrying capacity
 ## S : number of species
 ## NICH.WD : niche width, parameter 's' in "neut.sim.uni"
 ## DISP.WD : dispersal width, parameter 'u' in "neut.sim.uni"
 ## D : base probability of death
 ## B : base probability of birth
 ## MIGRATE : migration rate, parameter 'm' in "neut.sim.uni"
 ## IMMIGRATE : immigration rate
 ## FEC : fecunity - number of reproductive attempts per generation per individual
 ## FRACTAL : fractal dimensions to use
 ## CYCLES : number of cycles
 ## TIME : length of each cycle
 ## MIG : indicates if old (=0) or new (=1) migration algo should be used
 ## NPERM : number of permutations for randomizations of community file
 ## METH : method of permutations
 ## SP : if TRUE then species are each shuffled independently
 ## ALL : if TRUE then all relevant combinations of permutations are performed
 ## NPAR : the number of processors to run on
 ##'RPARGS: is a list of arguments that must be supplied if the random patterns algo is desired
 ###the arguments of RPARGS are used input into the function 'null.perms', they include:
 ###1)'allRP' if TRUE & 'all' = TRUE, then Random Patterns algo used as the spatial null
 ###2)'nstrata',3)'mtrials1',4)'mtrials2',5)'alpha',6)'npar'
 ## DESTDIRDRAW : directory where all raw results are returned (i.e. community sim file and vario files)
 ## SIMINPUT : this is a list of empirically derived simulation parameters
 #####################
 if(!is.null(SIMINPUT)){
  M = SIMINPUT$M
  S = SIMINPUT$S
  K = SIMINPUT$K
  no.bbs = SIMINPUT$no.bbs
  ENV = SIMINPUT$envir
  grad.len = dist(range(ENV,na.rm=TRUE))[1]
  FRACTAL = round(D.vario.est(ENV),2)
 }
 else
  grad.len = M
 coord <- cbind(rep(1:M,each=M),rep(1:M,times=M))
 Dist <- dist(coord)
 H <- round(Dist)
 hmax <- round(max(Dist)/2) 
 H[H > hmax] <- NA
 H <- as.vector(H)
 Dist <- sapply(split(Dist, H), mean)
 OUT <- NULL
 OBS <- array(NA,dim=c(REPS,length(Dist),4)) ##observed variogram value dim(replicate,distance,vals: exp, obs, pos, & neg 
 MEAN <- array(NA,dim=c(REPS,length(Dist),4)) ##avg of the permutation output
 SD <- array(NA,dim=c(REPS,length(Dist),4)) ##sd of the permutation output
 Z <- array(NA,dim=c(REPS,length(Dist),4)) ##Z-score of observed value relative to null permutations
 P.val <- array(NA,dim=c(REPS,length(Dist),4)) ##p.val of significance test
 Reject.025 <- array(NA,dim=c(REPS,length(Dist),4)) ##TRUE if significantly smaller than minimum confidence bound
 Reject.975 <- array(NA,dim=c(REPS,length(Dist),4)) ##TRUE if significantly larger than maximum confidence bound
 if(is.null(REP.ID)) REP.ID <- 1:REPS
 for(D in 1:length(FRACTAL)) {
  for(s in 1:length(NICH.WD)) {
   for(u in 1:length(DISP.WD)) { 
    for(m in 1:REPS) {
     if(is.character(ENV)){
      if(ENV == 'neutral'){
       HAB = 1
       OPTIM = 1
      }
      else{
       HAB <- gen.env(M = M, env = ENV, D = FRACTAL[D])
       OPTIM <- runif(S,.1,.9)*M #species ennviornmental optima
      }
     }
     else{
      HAB <- ifelse(is.na(ENV),0,ENV)
      OPTIM <- runif(S,range(ENV,na.rm=TRUE)) #species ennviornmental optima
     }
     if(is.null(SIMINPUT)){
      if(ENV == 'neutral')
       D.est <- NA
      else
       D.est <- D.vario.est(HAB)
     }
     else{
      D.est <- FRACTAL
     }
     SIM <- neut.sim.uni(M = M, K = K, S = S, s = NICH.WD[s], u = DISP.WD[u],b = BIRTH, d = DEATH, m = IMMIGRATE,habitat = HAB,
                         fitness = OPTIM,fec = FEC, time = TIME, cycles = CYCLES, mig = MIG, grad.len = grad.len)
     ##make SIM into a site x species pres/absence matrix to be input into vario
     mat <- apply(census(SIM,snap=length(SIM$snaps)),c(1,2),as.logical)*1
     mat <- mat[,apply(mat,2,sum)>0] ##drop species that never occur
     if(exists('no.bbs'))
       mat[as.vector(no.bbs),] <- -999 ##the code that identifies cells not sampled
     ###
     V.rand <- vario(mat,SIM$coords,pos.neg=FALSE) ##because we only need the expected variance
     VP.rand <- null.perms(mat,V.rand,NPERM,coords=SIM$coords,meth='random',npar=NPAR)
     V.spat <- vario(mat,SIM$coords,pos.neg=TRUE) 
     VP.spat <- null.perms(mat,V.spat,NPERM,coords=SIM$coords,meth='randpat',RPargs=RPARGS,npar=NPAR)
     ##create a comlomerate of the randomization results
     VP<- array(NA,dim=dim(VP.spat$vario)+c(0,1,0))
     VP[,1,]<-VP.rand$vario[,1,] ##expected val
     VP[,-1,]<-VP.spat$vario ## obs and pos/neg 
     VP[,2,] <- VP.spat$vario[,2,]+VP.spat$vario[,3,] ##add pos and neg fractions together to get the 'true' observed var
     if(!is.null(DESTDIRRAW)){
      ##save time intensive results 
      OBJS <- paste(PREFIX, 'M', M, 'D', FRACTAL[D], 'K', max(K), 'S', S, 's', NICH.WD[s], 'u', DISP.WD[u],'rep', REP.ID[m],'.Rdata', sep = '')
      save(SIM,VP.rand,VP.spat,file=paste(DESTDIRRAW,"/",OBJS,sep=''))
     }    
     ##add summarized results to output tables
     OBS[m,,] <- VP[,,1] 
     MEAN[m,,] <- apply(VP,1:2,mean)
     SD[m,,] <- apply(VP,1:2,sd)
     Z[m,,] <- (OBS[m,,] - MEAN[m,,]) / SD[m,,]
     for(k in 1:length(Dist)){ 
      for(i in 1:4){ 
       if(VP[k,i,1] > MEAN[m,k,i])
        P.val[m,k,i] <- sum(VP[k,i,1]<=VP[k,i,])/(NPERM+1)
       else
        P.val[m,k,i] <- sum(VP[k,i,1]>=VP[k,i,])/(NPERM+1)
     }}
     Reject.025[m,,] <- apply(VP,1:2,quantile,.025) ##lower confidence interval
     Reject.975[m,,] <- apply(VP,1:2,quantile,.975) ##upper confidence interval
     if(is.null(OUT)){##very first run through
      len<-length(Dist)*4
      len2<-length(Dist)
      OUT <- data.frame(M = rep(M,len), K = rep(max(K),len), S = rep(S,len), S.rel = rep(ncol(mat),len), D = rep(FRACTAL[D],len),
                        D.est = rep(D.est,len), s = rep(NICH.WD[s],len), s.rel = rep(round(NICH.WD[s]*grad.len,2),len), 
                        u = rep(DISP.WD[u],len), u.rel = rep(round(DISP.WD[u]*M,2),len), rep = rep(REP.ID[m],len),
                        mod = rep(c('rand',rep('spat',3)),each=len2), type = rep(c('exp','obs','pos','neg'),each=len2),
                        dist = rep(Dist,4),var = as.vector(OBS[m,,]), mean = as.vector(MEAN[m,,]), z = as.vector(Z[m,,]),
                        p = as.vector(P.val[m,,]), avg.rand = as.vector(MEAN[m,,]), CI.low = as.vector(Reject.025[m,,]),CI.up = as.vector(Reject.975[m,,]))
      
     }
     else{
      OUT <- rbind(OUT,
             data.frame(M = rep(M,len), K = rep(max(K),len), S = rep(S,len), S.rel = rep(ncol(mat),len), D = rep(FRACTAL[D],len),
                        D.est = rep(D.est,len), s = rep(NICH.WD[s],len), s.rel = rep(NICH.WD[s]*grad.len,len), 
                        u = rep(DISP.WD[u],len), u.rel = rep(DISP.WD[u]*M,len), rep = rep(REP.ID[m],len),
                        mod = rep(c('rand',rep('spat',3)),each=len2), type = rep(c('exp','obs','pos','neg'),each=len2),
                        dist = rep(Dist,4),var = as.vector(OBS[m,,]), mean = as.vector(MEAN[m,,]), z = as.vector(Z[m,,]),
                        p = as.vector(P.val[m,,]), avg.rand = as.vector(MEAN[m,,]), CI.low = as.vector(Reject.025[m,,]),CI.up = as.vector(Reject.975[m,,]))
             )
    }}
    if(length(REP.ID)==1){
     OBJS.sum1 <- paste(PREFIX, 'M', M, 'D', FRACTAL[D], 'K', max(K), 'S', S, 's', NICH.WD[s], 'u', DISP.WD[u],'rep', REP.ID,'.Rdata', sep = '') 
     OBJS.sum2 <- paste(PREFIX, 'M', M, 'D', FRACTAL[D], 'K', max(K), 'S', S, 's', NICH.WD[s], 'u', DISP.WD[u],'rep', REP.ID,'.csv', sep = '')
    }
    else{
     OBJS.sum1 <- paste(PREFIX, 'M', M, 'D', FRACTAL[D], 'K', max(K), 'S', S, 's', NICH.WD[s], 'u', DISP.WD[u],'.Rdata', sep = '') 
     OBJS.sum2 <- paste(PREFIX, 'M', M, 'D', FRACTAL[D], 'K', max(K), 'S', S, 's', NICH.WD[s], 'u', DISP.WD[u],'.csv', sep = '')
    }
    ##save summary data
    save(OBS,MEAN,Z,SD,P.val,Reject.025,Reject.975,file=paste(DESTDIR,"/",OBJS.sum1,sep=''))
    ##output a temporary csv file in case full job doesn't complete
    write.csv(OUT, file = paste(DESTDIR,"/",OBJS.sum2,sep=''))
 }}}
 ##save output dataframe if reps more that 1
 if(length(REP.ID)>1){
  OBJS.final <- paste(PREFIX, 'M', M, 'D', FRACTAL[D], 'K', max(K), 'S', S,'.csv', sep = '')
  write.csv(OUT, file = paste(DESTDIR,"/",OBJS.final,sep=''))
 }
 OUT
}

##4.2##
sim.dat.prep<-function(stcrd,psp,envir,perc.land,res,n,K,K.cont=FALSE,export=FALSE,DESTDIR=getwd()){
 ##################################################################
 ##PURPOSE:take a starting coordinate and generate the empirically constrained
 ##input files for the simulation framework. This script generates and Rdata file if export is TRUE
 ##ARGUMENTS:
 ##stcrd: the starting x and y coordinate
 ##psp: the sp x row x col array
 ##envir: the enviornmental grid to simulated over
 ##perc.land: the prop of a pixel that must fall on land to be considered good
 ##res: the spatial resolution of the grid
 ##n: the spatial extend in number of pixels of one side of the grid
 ##K: the baseline carrying capacity for a given pixel
 ##K.cont: TRUE if continously varying Ks are desired which are linear related to envir
 ##export: TRUE if you wish to export the .Rdata file
 ##DESTDIR: the directory where the .Rdata file will be written to
 ##################################################################
 i<-(stcrd[1]/res)+1
 j<-(stcrd[2]/res)+1
 envir.grd<-envir[1,i:(i+n-1),j:(j+n-1)]
 psp.grd<-psp[,i:(i+n-1),j:(j+n-1)]
 ##drop species that do not occur in this grid##
 psp.grd<-psp.grd[apply(psp.grd,1,sum,na.rm=TRUE)>0,,]
 S.grd<-dim(psp.grd)[1]
 if(K.cont)
  K<-rep(as.vector(K*envir[2,i:(i+n-1),j:(j+n-1)]),each=S.grd)
 else
  K<-rep(as.vector(ifelse(envir[2,i:(i+n-1),j:(j+n-1)]>=perc.land,K,0)),each=S.grd) ##specifies which pixels will be allowed indviduals
 no.bbs<-ifelse(apply(psp.grd,c(2,3),sum)==0,TRUE,FALSE) ##specifies which pixels will be included in randomizations
 sim.input<-list(M=n,S=S.grd,K=K,no.bbs=no.bbs,envir=envir.grd)
 if(export)
  save(sim.input,file=paste(DESTDIR,"/SimEmpirInputs.M",n,".x",stcrd[1],".y",stcrd[2],".Rdata",sep=""))
 else
  sim.input
}

##4.3##
##big.graph<-function(out,...){
 ##PURPOSE: to graph the results of a call to var.test
 ##which is a mechanism for generating batch output of the sim and 


#####Part V - FUNCTIONS FOR EMPIRICAL BBS PROCESSING#####
##5.1##
GridCount<-function(samps,grd,land=NULL,perc.land=NULL){
 ##count number of samples in each cell 
 ##purpose: to count the number of samples that fall in a given grid section
 ##samps is the x,y coords of each sample
 ##grd is the x,y coords of the grid
 ##land is the land array which tells us where there is land 
 ##perc.land the prop of the pixel that must be on land to be a good pixel
 dat.x<-samps[,1]; dat.y<-samps[,2]
 grd.x<-grd[,1]; grd.y<-grd[,2]
 uni.x<-sort(unique(grd.x)); uni.y<-sort(unique(grd.y))
 luni.x<-length(uni.x); luni.y<-length(uni.y)
 if(is.null(land)){
  land<-matrix(0,nrow=luni.x-1,ncol=luni.y-1)
  perc.land<-0
 }
 n<-matrix(0,nrow=luni.x-1,ncol=luni.y-1)
 for(i in 1:(luni.x-1)){
  for(j in 1:(luni.y-1)){
   if(land[i,j]>=perc.land){
    n[i,j]<-sum(dat.x >= uni.x[i] & dat.x < uni.x[i+1] & 
                dat.y >= uni.y[j] & dat.y < uni.y[j+1] )
 }}}     
 out<-list(x = uni.x, y = uni.y, n = n)
 out
}

##5.2##
StartCoord<-function(samps,grd,n,nsamp,psamp,land=NULL,perc.land=NULL){
 ##purpose: to generate SW starting coordinates for a grid
 ##samps is the x,y coords of each sample
 ##grd is the x,y coords of the grid
 ##n the size of the grid
 ##nsamp the minimum number of atlas points
 ##psamp the prop of filled grid cells to be considered complete
 ##land is the land array which tells us where there is land 
 ##perc.land the prop of the pixel that must be on land to be a good pixel
 if(is.null(land))
  cnt.map<-GridCount(samps,grd)
 else
  cnt.map<-GridCount(samps,grd,land,perc.land)
 startcoord<-NULL
 for(i in 1:(length(cnt.map$x)-n)){
  for(j in 1:(length(cnt.map$y)-n)){ 
   prop.fill<-sum(cnt.map$n[i:(i+n-1),j:(j+n-1)] >= nsamp)/n^2
   if(prop.fill >= psamp)
    startcoord<-rbind(startcoord,c(cnt.map$x[i],cnt.map$y[j],prop.fill))
 }}
 startcoord
}  

##5.3##
Empir<-function(stcrds,xmn,ymn,len,n,pos.neg=TRUE,meth='randpat',nstrata=3,rperm=1e3,nperm=299,npar=8,all.dat=psp.bbs,linux=TRUE){
 ##purpose: to generate empirical variogram analyses for a given set of starting coordinates
 ##arguments
 ##'stcrds' :the spatial staring coordinates
 ##'xmn': the minimum x value of the spatial grid, for place keeping on the species data
 ##'ymn': the minimum y value of the spatial grid, for place keeping on the species data
 ##'len' : the spatial distance between grid cells
 ##'n' : the size of the grid along one spatial dimension
 ##"meth" the type of permutation to use, options include:
 ##"reflect": random reflection/rotations of species (only makes sence when sp are not fixed
 ###"shift": random torodial shifting with or with sp fixed
 ###"both": both reflection and shifting
 ###"random": random shuffle
 ###"randpat": random patterns algo of Roxburgh and Chesson 1998, must parameterize RPargs (See below)
 ##'nstrata : the number of strata to devide the grid into along one axis for randomization
 ##rperm: the number of permutations to use for strata swapping
 ##nperm: the number of permutations to use for calculating variogram envalope
 ##npar : the number of processors to use for calculating all of the many variograms
 ##all.dat: the global empirical S x n x n array of occupancies
 ##linux: whether or not this is to be run on a linux machine
 ###############################################################
 ####for each starting coordinate
 for(icord in 1:length(stcrds[,1])){
  xmin<-stcrds[icord,1] ; xmax<-xmin+n*len
  ymin<-stcrds[icord,2] ; ymax<-ymin+n*len
  xs<-seq(xmin,xmax,len)
  ys<-seq(ymin,ymax,len)
  psp<-all.dat[,((xmin-xmn+len)/len):((xmax-xmn)/len),((ymin-ymn+len)/len):((ymax-ymn)/len)]
  ##drop species that do not occur in this grid
  psp<-psp[apply(psp,1,sum)>0,,]
  S<-dim(psp)[1]
  mat<-matrix(as.vector(psp),ncol=S,nrow=n^2,byrow=TRUE)
  coords<-expand.grid(xs[-(n+1)],ys[-(n+1)])
  v<-vario(mat,coords,pos.neg=pos.neg,grain=len) 
  if(meth=='randpat')
   rpargs<-list(TRUE,nstrata,rperm,0,0.01,1)
  else
   rpargs<-list(FALSE,nstrata,rperm,0,0.01,1)
  out<-null.perms(x=psp,vobject=v,coords=coords,nperm=nperm,meth=meth,npar=npar,RPargs=rpargs,linux=linux)
  save(out,file=paste("empir.n",n,".",meth,".x",xmin,".y",ymin,".Rdata",sep=""))
 }
 out
}


##5.4##
Empir.pdf<-function(stcrds,n,nsamp=2,both=TRUE){
 ###purpose: to generate .pdfs of the empirical variograms####
 if(both)
  pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"/empir_n",n,".pdf",sep=""),width=7*3,height=7*2)
 else
  pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"/empir_n",n,".pdf",sep=""),width=7*2.5,height=7)
 icount<-0
 for(icord in 1:length(stcrds[,1])){
  out<-NULL
  try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"/empir.n",n,".x",
      stcrds[icord,1],".y",stcrds[icord,2],".Rdata",sep="")),silent=TRUE)
  if(!is.null(out)){
   if(icount==0){
    quants<-apply(out$vario,1:2,quantile,c(.025,.975))
    ylims<-range(list(quants,out$vario[,,1]))
   }
   if(both)   
    par(mfrow=c(2,3))
   else
    par(mfrow=c(1,3))
   brks<-seq(nsamp,85,5)
   cls<-terrain.colors(length(brks)-1)
   image(Bn100$x,Bn100$y,Bn100$n,breaks=brks,col=cls,ylim=c(-100,3500),xlim=c(-100,5500),main=paste('min samp',nsamp),xlab='',ylab='')
   points(stcrds[icord,1],stcrds[icord,2],col='purple',pch=19,cex=3)
   polygon(c(rep(stcrds[icord,1],2),rep(stcrds[icord,1]+n*len1,2)),c(stcrds[icord,2],rep(stcrds[icord,2]+n*len1,2),stcrds[icord,2])) 
   lines(wrld.pts,col='red',lwd=2)
   v.graph.all(out,plot=FALSE,ylims=ylims,main=paste('geo coords: x=',stcrds[icord,1],' y=',stcrds[icord,2],sep=''))
   if(both){ ##also generate delta plots
    plot(1:10,1:10,type='n',axes=F,xlab='',ylab='',frame.plot=F)
    ##calculate delta
    diff<-array(rep(out$vario[,,1],dim(out$vario)[3]),dim=dim(out$vario)) - out$vario
    DIFF<- apply(diff,1:2,mean)
    sds<-apply(diff,1:2,sd)
    DELTA <- DIFF/sds
    if(icount==0)
     ylims2<-c(0,max(DELTA))
    plot(out$vdist,DELTA[,1],ylim=ylims2,type='o',main=paste('geo coords: x=',stcrds[icord,1],' y=',stcrds[icord,2],sep=''))
    plot(out$vdist,DELTA[,2],ylim=ylims2,type='o',col='red',main=paste('geo coords: x=',stcrds[icord,1],' y=',stcrds[icord,2],sep=''))
    lines(out$vdist,DELTA[,3],type='o',col='blue')
    legend('bottomleft',c('pos','neg'),col=c('red','blue'),pch=19,lty=1,lwd=2,bty='n')
   }
   icount<-1
 }}
 dev.off()
}

##5.4.b##
Empir.pdf.b<-function(stcrds,n,nsamp=2,both=FALSE){
 ###purpose: to generate .pdfs of the empirical variograms####
 if(both)
  pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"_random/empir_n",n,"_random.pdf",sep=""),width=7*3,height=7*2)
 else
  pdf(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"_random/empir_n",n,"_random.pdf",sep=""),width=7*2.5,height=7)
 icount<-0
 for(icord in 1:length(stcrds[,1])){
  for(i in 1:2){
   out<-NULL
   if(i == 1){
    try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"_random/empir.n",n,".random.x",
        stcrds[icord,1],".y",stcrds[icord,2],".Rdata",sep="")),silent=TRUE)
    ##plot the geographical map##
    if(both)   
     par(mfrow=c(2,3))
    else
     par(mfrow=c(1,3))
    brks<-seq(nsamp,85,5)
    cls<-terrain.colors(length(brks)-1)
    image(Bn100$x,Bn100$y,Bn100$n,breaks=brks,col=cls,ylim=c(-100,3500),xlim=c(-100,5500),main=paste('min samp',nsamp),xlab='',ylab='')
    points(stcrds[icord,1],stcrds[icord,2],col='purple',pch=19,cex=3)
    polygon(c(rep(stcrds[icord,1],2),rep(stcrds[icord,1]+n*len1,2)),c(stcrds[icord,2],rep(stcrds[icord,2]+n*len1,2),stcrds[icord,2])) 
    lines(wrld.pts,col='red',lwd=2)
    v.graph.all(out,plot=FALSE,draw='exp',main=paste('geo coords: x=',stcrds[icord,1],' y=',stcrds[icord,2],sep=''))
   }
   if(i == 2){
    try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"/empir.n",n,".x",
        stcrds[icord,1],".y",stcrds[icord,2],".Rdata",sep="")),silent=TRUE)
    v.graph.all(out,plot=FALSE,draw='pos-neg',main=paste('geo coords: x=',stcrds[icord,1],' y=',stcrds[icord,2],sep=''))
 }}}
 dev.off()
}

##5.5##
Empir.fit<-function(stcrds,n,sim.array,meth='avg',on='desktop'){
 ##to compare the fit of the simulated and empirical results
 output<-NULL
 for(icord in 1:length(stcrds[,1])){
  if(on=='desktop'){
   try(load(paste("C:/Users/dmcglinn/Documents/Lab data/abundance surface/BBS/empir_study_n",n,"/empir.n",n,".x",
       stcrds[icord,1],".y",stcrds[icord,2],".Rdata",sep="")),silent=TRUE)
  }
  if(on=='laptop'){
   try(load(paste("C:/Users/Daniel McGlinn/Desktop/Lab data/abundance surface/BBS/empir_study_n",n,"/empir.n",n,".x",
       stcrds[icord,1],".y",stcrds[icord,2],".Rdata",sep="")),silent=TRUE)
  }
  if(!is.null(out)){
   null.means<-apply(out$vario,c(1,2),mean)
   null.sds<-apply(out$vario,c(1,2),sd)
   d.pos<-(out$vario[,2,1]-null.means[,2])/null.sds[,2]
   d.neg<-(out$vario[,3,1]-null.means[,3])/null.sds[,3]
   ##drop unneed distances from sim.array
   sim.array<-sim.array[,1:length(d.pos),,,]
   obs.array<-array(rbind(d.neg,d.pos),dim=c(2,length(d.pos)))
   obs.array<-array(rep(obs.array,prod(dim(sim.array[1,1,,,]))),dim=dim(sim.array))
   resid<-obs.array-sim.array
   ##now need to find the parameters that are closest to zero
   ##after summing across all the distances and two types of indices
   tot.resid<-apply(abs(resid),3:5,sum)
   ##now unroll tot.resid into 4 columns
   dat<-cbind(rep(unique(D),25),rep(rep(unique(s),5),each=3),rep(unique(u),each=15),tot.resid)
   colnames(dat)<-c('D','s','u','tot.resid')
   best<-dat[order(dat[,4])[1:5],] ##best 5
   worst<-dat[order(dat[,4],decreasing=TRUE)[1:5],]##worst 5
   ##now just need to output what we want
   x<-stcrds[icord,1]+(n*100)/2
   y<-stcrds[icord,2]+(n*100)/2
   if(meth=='avg')
    output<-rbind(output,c(1,x,y,apply(best,2,mean),apply(worst,2,mean)))
   if(meth=='top')
    output<-rbind(output,c(1,x,y,best[1,],worst[1,]))
   if(meth=='all')
    output<-rbind(output,cbind(1:5,rep(x,5),rep(y,5),best,worst))
 }}
 colnames(output)<-c('rank','x','y','Dbest','sbest','ubest','resid.best',
                    'Dwor','swor','uwor','resid.wor')
 output
}




