##Part II - FUNCTIONS FOR Spatial Permutations-----------------------------------

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
  n<-dim(psp)[2]
  if(length(dim(psp))==3)
    S<-dim(psp)[1]
  else{
    S<-1
    psp<-array(psp,dim=c(S,n,n))
  }
  ##first prepare psp for the randomization process
  ##fill in empty pixels with -999
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

## ANALYZING AND GRAPHING RESULTS----------------------------

##3.1##
getCovFractions = function(x)
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

get_breaks = function(breaks, hmin, hmax, maxDist, log=FALSE) {
  ## compute breaks of distance bins used by the function vario()
  ## Arguments:
  ## breaks: either a vector of breaks or an integer number of 
  ##  breaks to compute.  If a vector of breaks is supplied then
  ##  the function returns that exact same vector back. 
  ## hmin: minimum distance of interest 
  ## hmax: maximum distance of interest
  ## maxDist: maximum possible distance (greater or equal to hmax)
  ## log: boolean, if true the breaks are equidistance on a log scale
  if (length(breaks) == 1) {
    if (log) {
      if (round(hmax, 2) == round(maxDist / 2, 2)) {
        incre = (hmax - hmin) / breaks
        hmax = hmax + incre
      }
      breaks = exp(seq(log(hmin), log(hmax), length.out=breaks))
    }
    else 
      breaks = seq(hmin, hmax, length.out=breaks)    
  }
  return(breaks)
}

check_vario_direction_args = function(direction = 'omnidirectional',
                                      tolerance = pi/8,
                                      unit.angle = c('radians', 'degrees'))
{
  ## This function carries out sanity checks on the directional
  ## arguments that are supplied to the function vario(), if 
  ## these checks are failed then vario() will stop with an error message
  ## Note: this code was copied from geoR in the function variog
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
}
##3.2##
vario = function(x, coord, grain=1, breaks=NA, log=FALSE, hmin=NA,
                 hmax=NA, round.int=FALSE, pos.neg=FALSE, binary=TRUE,
                 snap=NA, median=FALSE, quants=NA, direction = 'omnidirectional',
                 tolerance = pi/8, unit.angle = c('radians', 'degrees'),
                 distance.metric = 'euclidean', univariate=FALSE)
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
  ## log: boolean, if true then the breaks are equidistance  
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
  ## univariate: if TRUE then results are computed on a per species basis
  if (class(x) == "sim"){
    coord = x$coords
    if (is.na(snap))
      snap = length(sim$snaps)
  }
  else if (is.vector(x))
    x = matrix(x)
  else if (!is.matrix(x))
    stop('x must be either of class sim, a vector, or a matrix')
  #x = ifelse(x == -999, NA, x)  ## best to fix these before entering into the function
  if (univariate) {
    vobject = vario_uni(x, bisect=FALSE, coord, grain, breaks, log, hmin, hmax,
                        round.int, pos.neg, binary, snap, median, quants, direction,
                        tolerance, unit.angle, distance.metric)
  }
  else {
    if (distance.metric != 'euclidean') {
      if (pos.neg)
        stop("cannot commpute pos-neg covariance using a turnover metric")
      else
        require(vegan)
    }
    unit.angle = match.arg(unit.angle)
    check_vario_direction_args(direction, tolerance, unit.angle)
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
        hmin = min(Dist)  ## potentially this should be set when breaks is NA as well
      if (is.na(hmax))
        hmax = maxDist / 2
      H = Dist
      breaks = get_breaks(breaks, hmin, hmax, maxDist, log)
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
    class(vobject) = 'vario'
    vobject$parms = data.frame(grain, hmin, hmax, S=S, N=N, pos.neg, median, direction,
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
                 res = as.double(rep(0, length(as.vector(Dist)))))$res
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
  }
  return(vobject)
}
  
vario_uni = function(x, bisect=FALSE, ...)
{
  ## Purpose: to compute the multivariate variogram as well as individual univariate 
  ## variograms for evey column of x
  ## Arguments:
  ## x : site x sp matrix
  ## bisect : if the bisection style vairogram should be computed
  ## ... : arguments supplied to the function vario()
  ## Note: speed gains would be significant if partitioning of computation between
  ## species was carried out within the vario function after computing the
  ## distance matrix because that is a time intensive step
  require(snowfall)
  S = ncol(x)
  if (bisect)
    v = vario_bisect(x, ...)
  else
    v = vario(x, ...)
  n_cpus = length(suppressMessages(sfGetCluster()))
  if (n_cpus > 0) {
    sfSource('./scripts/spat_functions.R')
    sfLibrary(vegan)
    if (bisect) 
      exp_var = sfSapply(1:S, function(sp) vario_bisect(x[ , sp], ...)$vario$exp.var)
    else
      exp_var = sfSapply(1:S, function(sp) vario(x[ , sp], ...)$vario$exp.var)
  }
  else {
    if (bisect)
      exp_var = sapply(1:S, function(sp) vario_bisect(x[ , sp], ...)$vario$exp.var) 
    else
      exp_var = sapply(1:S, function(sp) vario(x[ , sp], ...)$vario$exp.var) 
  }
  colnames(exp_var) = paste('sp', 1:S, sep='')
  v$exp.var = exp_var  
  return(v)
}


##3.3##
null.perms = function(x, vobject, nperm, coords=NULL, meth='both',
                      sp=TRUE, all=FALSE, snap=NULL, npar=1, linux=FALSE,
                      RPargs=FALSE, breaks=NA) {
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
  ##"snap" is the time period of the simulation to analyze, defaults to NULL val, if not set gets internally set to last time period
  ##"npar" = number of processors to run the function on
  ##"median" if TRUE then median is also calculated in addition to mean
  ##"linux" if TRUE then function assumes you are on a linux cluster and therefore exports a different compiled code
  ##"RPargs" is a vector of arguments that are needed to perform the Random Patterns spatial null model
  ###see the notes associated with the function 'null.gen' that indicate how "RPargs" should be parameterized
  ##"breaks" gives either the number or the position of the breaks for the function vario
  ##Note: "meth" and "sp" are arguments to randomization function "SpatPerm2D"
  dists = vobject$vario$Dist
  grain = vobject$parms$grain
  hmin = vobject$parms$hmin
  hmax = vobject$parms$hmax
  pos.neg = vobject$parms$pos.neg
  median = vobject$parms$median
  if (class(vobject$parms$direction) == "factor") 
    direction = as.character(vobject$parms$direction)
  else
    direction = as.numeric(vobject$parms$direction)
  tolerance = vobject$parms$tolerance
  unit.angle = as.character(vobject$parms$unit.angle)
  if (is.na(unit.angle))
    unit.angle = 'degrees'
  distance.metric = as.character(vobject$parms$distance.metric)
  if (class(x) == 'sim') {
    coords = x$coords
    if (is.null(snap))
      snap = length(x$snaps)
  }
  else if (is.null(coords))
    stop('need to supply spatial coordinates if not a simulation product')
  r.vals = list()
  r.vals$parms = vobject$parms
  r.vals$p = vobject$p
  if (median & !pos.neg)
    stop("if computing medians must also compute pos.neg fractions, set pos.neg=TRUE")
  if (pos.neg) { ##pos and neg fractions
    if (all) { ##all relevant null models
      if (median) { ##will compute mean and median
        r.vals$vario = array(0, dim=c(length(dists), 6, 3, nperm + 1))##dists, results, methods, perms
        r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos',
                                                              'neg', 'exp.med',
                                                              'pos.med', 'neg.med')]) 
      }
      else { ##only compute means
        r.vals$vario = array(0, dim=c(length(dists), 3, 3, nperm+1))##dists, results, methods, perms
        r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos', 'neg')])
      }
    }  
    else { ##only a single null used
      if (median) {
        r.vals$vario = array(0, dim=c(length(dists), 6, nperm + 1))##dists, results, perms
        r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos',
                                                            'neg', 'exp.med',
                                                            'pos.med', 'neg.med')]) 
      }
      else {
        r.vals$vario = array(0, dim=c(length(dists), 3, nperm + 1))##dists, results, perms
        r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('obs.var', 'pos', 'neg')])
      }
    }
  }
  else {##only exp and obs fractions
    if (all) {
      if (median) { ##will compute mean and median
        r.vals$vario = array(0, dim=c(length(dists), 4, 3, nperm + 1))
        r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var',
                                                              'obs.med', 'exp.med')])
      }
      else {
        r.vals$vario = array(0, dim=c(length(dists), 2, 3, nperm + 1))
        r.vals$vario[ , , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var')])
      }
    }
    else { ##only a single null used
      if (median) {
        r.vals$vario = array(0, dim=c(length(dists), 4, nperm + 1))
        r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var',
                                                            'obs.med', 'exp.med')])
      }
      else { 
        r.vals$vario = array(0, dim=c(length(dists), 2, nperm + 1))
        r.vals$vario[ , , 1] = as.matrix(vobject$vario[, c('exp.var', 'obs.var')])
      }
    }
  }
  if (class(x) == 'sim') {
    pop = as.logical(x$snaps[[snap]]) ##converts it to a pres/absence vector
    dim(pop) = c(x$p$S, x$p$M, x$p$M)
  }
  else {
    pop = array(x, dim=c(sqrt(nrow(x)), sqrt(nrow(x)), ncol(x)))
    pop = aperm(pop, c(3, 1, 2))
  }
  if (RPargs[[1]] & npar == 1) {
    r.vals$p.conv1 = 0 ##average proportion of species that converged with strata swaps
    r.vals$p.conv2 = 0 ##average proportion of species that converged with pixel swaps
  }
  if (npar == 1) { ##all permutations option not yet implemented for 1 processor
    pb = txtProgressBar(min = 0, max = nperm, style = 3)
    for (i in 1:nperm) {
      if (RPargs[[1]]) { ##use the random pattern algo for the spatial null
        out = RandPatPar(psp=pop, nstrata=RPargs[[2]], mtrials1=RPargs[[3]],
                         mtrials2=RPargs[[4]], alpha=RPargs[[5]], npar=RPargs[[6]])
        S = dim(pop)[1]
        n2 = dim(pop)[2] + 2
        r.vals$p.conv1 = r.vals$p.conv1 + sum(out[2, ] <= RPargs[[5]]) / S / nperm
        r.vals$p.conv2 = r.vals$p.conv2 + sum(out[4, ] <= RPargs[[5]], na.rm=TRUE) / S / nperm
        rpop = array(0, dim=c(S, n2, n2))
        for (k in 1:S) {
          rpop[k, , ] = array(out[-(1:5), k], dim=c(n2, n2))
        }
        rpop = rpop[ , -c(1, n2), -c(1, n2)]
      }
      else {
        rpop = SpatPerm2D(pop, meth=meth, sp=sp)
        rpop = FixUnSamp2(pop, rpop)
      }
      rmat = apply(rpop, 1, as.vector) ##converts to a M^2 x S matrix - same effect as loop in 'census' function 
      rv = vario(x=rmat, coord=coords, grain=grain, breaks=breaks, hmin=hmin,
                 hmax=hmax, pos.neg=pos.neg, median=median, direction=direction,
                 tolerance=tolerance, unit.angle=unit.angle,
                 distance.metric=distance.metric)$vario
      if (pos.neg) {
        if (all) {
          if (median) {
            r.vals$vario[ , , , i + 1] = as.matrix(rv[ , ,
                                                      c('obs.var','pos','neg',
                                                        'exp.med','pos.med',
                                                        'neg.med')])
          }
          else {
            r.vals$vario[ , , , i + 1] = as.matrix(rv[ , , c('obs.var', 'pos', 'neg')])
          }
        } 
        else {
          if (median) {
            r.vals$vario[ , , i + 1] = as.matrix(rv[ , ,
                                                    c('obs.var','pos','neg',
                                                      'exp.med','pos.med',
                                                      'neg.med')])
          }  
          else {
            r.vals$vario[ , , i + 1] = as.matrix(rv[ , c('obs.var', 'pos', 'neg')])
          }  
        }
      }
      else {
        if (all) {
          if (median) {
            r.vals$vario[ , , , i + 1] = as.matrix(rv[ , , c('exp.var', 'obs.var',
                                                             'obs.med', 'exp.med')])
          }
          else {
            r.vals$vario[ , , , i + 1] = as.matrix(rv[ , , c('exp.var', 'obs.var')])
          }
        }
        else {
          if (median) {
            r.vals$vario[ , , i + 1] = as.matrix(rv[ , c('exp.var', 'obs.var',
                                                         'obs.med', 'exp.med')])
          }
          else {
            r.vals$vario[ , , i + 1] = as.matrix(rv[ , c('exp.var', 'obs.var')])
          }
        }
      }
      #print(i)
      Sys.sleep(0.1)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  else { ##computing in parallel
    require(snowfall)
    sfInit(parallel=TRUE, cpus=npar, type="SOCK")
    sfClusterSetupRNG()
    sfExport("pop", "vobject", "coords", "meth", "all", "sp", "RPargs",
             "breaks", "RandPatPar", "RandPat", "FixUnSamp", "FixUnSamp2",
             "SpatPerm2D", "SpatPerm2D.str", "vario", "getCovFractions",
             "null.gen", "check_vario_direction_args")
    if (linux)
      sfClusterEval(dyn.load("danspkg.so"))
    else {
      sfLibrary(danspkg)
    }
    out = unlist(sfLapply(1:nperm, function(...) null.gen(pop=pop, vobject=vobject,
                                                          coords=coords, meth=meth,
                                                          sp=sp, all=all, RPargs=RPargs,
                                                          breaks=breaks)))
    sfStop()
    if (pos.neg) {
      if (all) {
        if (median) {
          dim(out) = c(length(dists), 6, 3, nperm)
          r.vals$vario[ , , , -1] = out 
        }
        else {
          dim(out) = c(length(dists), 3, 3, nperm)
          r.vals$vario[ , , , -1] = out 
        }
      }
      else {
        if (median) {
          dim(out) = c(length(dists), 6, nperm)  
          r.vals$vario[ , , -1] = out     
        }
        else {
          dim(out)=c(length(dists), 3, nperm)  
          r.vals$vario[ , , -1] = out     
        }
      }
    }
    else {
      if (all) {
        if (median) {
          dim(out) = c(length(dists), 4, 3, nperm)  
          r.vals$vario[ , , , -1] = out 
        }
        else {
          dim(out) = c(length(dists), 2, 3, nperm)  
          r.vals$vario[ , , , -1] = out 
        }
      }
      else {
        if (median) {
          dim(out) = c(length(dists), 4, nperm)  
          r.vals$vario[ , , -1] = out 
        }
        else {
          dim(out) = c(length(dists), 2, nperm)  
          r.vals$vario[ , , -1] = out 
        }
      }
    }
  }
  r.vals$perm = TRUE
  r.vals$vdists = vobject$vario$Dist
  return(r.vals)
}

shuffle_comm = function(comm, swap) {
  ## Purpose: to returned a shuffled community site x species matrix
  ## Arguments:
  ## comm: site x species matrix with abundance or pres/absen data
  ## swap: two options: 'indiv' or 'sample' for individual or sample-based shuffling
  if (swap == 'indiv') {
    nquad = nrow(comm)
    comm_shuffled = comm
    for (j in 1:ncol(comm)) {
      rand_samp = sample(nquad, size=sum(comm[ , j]), replace=T)
      comm_shuffled[ , j] = as.numeric(table(c(rand_samp, 1:nquad)) - 1)
    }
  }
  else if (swap == 'sample')
    comm_shuffled = comm[sample(nrow(comm)), ]
  else
    stop('swap must be "indiv" or "sample"')
  return(comm_shuffled)
}

random_shuffle = function(x, vobject, swap, nperm, coords=NULL, breaks=NA) {
  ## Purpose: to generate individual or sample-based random shuffle statistical
  ## null expectations for variograms. 
  ## Note: With the sample-based shuffling the mean of the null distribution is
  ## the same as the average variance across all distance classes. So shuffling
  ## is only useful in deriving the variance around that expectation.
  ## Note: if x is a presence-absence matrix then the individual and sample-based
  ## approaches will generate identical within-species variograms
  ##Arguments:
  ##"x" is either an output of class 'sim' that is the output of 'sim.neut.uni' OR an site x species matrix
  ##"vobject" is the output of 'vario', this informs the function of what parameterization of vario to use
  ###specifically it indiates if the pos.neg components and median should be calculated
  ##"swap" two options: 'indiv' or 'sample' for individual or sample-based shuffling
  ### repsectively. 
  ##"nperm" is the number of permutations
  ##"coords" the spatial coordinates of the sites, not needed if x is of class 'sim'
  ##"breaks" what spatial breaks to use
  require(snowfall)
  npar = length(suppressMessages(sfGetCluster()))
  dists = vobject$vario$Dist
  grain = vobject$parms$grain
  hmin = vobject$parms$hmin
  hmax = vobject$parms$hmax
  pos.neg = vobject$parms$pos.neg
  median = vobject$parms$median
  if (class(vobject$parms$direction) == "factor") 
    direction = as.character(vobject$parms$direction)
  else
    direction = as.numeric(vobject$parms$direction)
  tolerance = vobject$parms$tolerance
  unit.angle = as.character(vobject$parms$unit.angle)
  if (is.na(unit.angle))
    unit.angle = 'degrees'
  distance.metric = as.character(vobject$parms$distance.metric)
  if (is.null(coords))
    stop('need to supply spatial coordinates if not a simulation product')
  r.vals = list()
  r.vals$parms = vobject$parms
  r.vals$p = vobject$p
  if (npar > 0) { ##all permutations option not yet implemented for 1 processor
    sfClusterSetupRNG()
    sfSource('./scripts/spat_functions.R')
    sfExport('x', 'swap', 'coords', 'distance.metric')
    if (direction == 'bisection') {
      rv = sfSapply(1:nperm, function(...)
                    vario_bisect(shuffle_comm(x, swap), coords,
                                 distance.metric=distance.metric)$vario[ , 'exp.var'])
    }  
    else {
      sfExport('grain','breaks','hmin','hmax','pos.neg','median','direction',
               'tolerance','unit.angle')
      rv = sfLapply(1:nperm, function(...)
        vario(shuffle_comm(x, swap), coords, grain=grain, breaks=breaks, hmin=hmin,
              hmax=hmax, pos.neg=pos.neg, median=median, direction=direction,
              tolerance=tolerance, unit.angle=unit.angle,
              distance.metric=distance.metric)$vario[ , c('exp.var', 'obs.var')])
    }  
  }
  else {
    if (direction == 'bisection')
      rv = sapply(1:nperm, function(...)
        vario_bisect(shuffle_comm(x, swap), coords,
                     distance.metric=distance.metric)$vario[ , 'exp.var'])
    else
      rv = lapply(1:nperm, function(...)
        vario(shuffle_comm(x, swap), coords, grain=grain, breaks=breaks, hmin=hmin,
              hmax=hmax, pos.neg=pos.neg, median=median, direction=direction,
              tolerance=tolerance, unit.angle=unit.angle,
              distance.metric=distance.metric)$vario[ , c('exp.var', 'obs.var')])
  }
  if (is.matrix(rv)) 
    rv_dim = dim(rv)
  else if (is.list(rv))
    rv_dim = dim(rv[[1]])
  r.vals$vario = array(NA, dim= rv_dim + c(0, 1))
  if (direction == 'bisection') {
    r.vals$vario[ , 1] = as.matrix(vobject$vario[ , 'exp.var'])
    r.vals$vario[ , -1] = rv
  }
  else {
    r.vals$vario[ , , 1] = as.matrix(vobject$vario[ , c('exp.var', 'obs.var')])
    r.vals$vario[ , , -1] = array(unlist(rv), dim=c(rv_dim, nperm))  
  }
  r.vals$perm = TRUE
  r.vals$vdists = vobject$vario$Dist
  return(r.vals)
}

    
##3.4##
null.gen<-function(pop,vobject,coords,meth,sp,all=FALSE,RPargs=FALSE,median=FALSE,breaks=NA){
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
 ##"breaks" gives either the number or the position of the breaks for the function vario
 ##Note: "meth" and "sp" are arguments to randomization function "SpatPerm2D"
 S<-dim(pop)[1]
 n<-dim(pop)[2]
 n2<-n+2
 grain = vobject$parms$grain
 hmin = vobject$parms$hmin
 hmax = vobject$parms$hmax
 pos.neg = vobject$parms$pos.neg
 median = vobject$parms$median
 if(class(vobject$parms$direction) == "factor") 
  direction = as.character(vobject$parms$direction)
 else
  direction = as.numeric(vobject$parms$direction)
 tolerance = vobject$parms$tolerance
 unit.angle = as.character(vobject$parms$unit.angle)
 distance.metric = as.character(vobject$parms$distance.metric)
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
   rv<-vario(x=rmat,coord=coords,grain=grain,breaks=breaks,hmin=hmin,hmax=hmax,
             pos.neg=pos.neg,median=median,direction=direction,tolerance=tolerance,
             unit.angle=unit.angle,distance.metric=distance.metric)$vario
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
  rv<-vario(x=rmat,coord=coords,grain=grain,breaks=breaks,hmin=hmin,hmax=hmax,
            pos.neg=pos.neg,median=median,direction=direction,tolerance=tolerance,
            unit.angle=unit.angle,distance.metric=distance.metric)$vario
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

mat2psp = function(sp_mat, xy_coord, N=NULL, M=NULL)
{
  ##place site by species matrix (x) into an S x N x M array where N >= M
  ##a multidimensional array that eases computing 
  ##replaces the old function 'grid.pres'
  ##Note: if area rectangular the first dimension of xy_coord does NOT
  ## have to be larger than the second dimension
  if (is.null(N) )
    N = sqrt(nrow(sp_mat))
  if (is.null(M) )
    M = N
  if (nrow(sp_mat) != nrow(xy_coord))
    stop('Number of samples in species matrix must match number of samples in xy-coordinates')
  if (N < M)
    stop('N should be >= M')
  if (N * M != nrow(sp_mat))
    stop('Number of specified samples (N X M) must be equal to the number of samples in sp_mat')
  S = ncol(sp_mat)
  ## order rows of sp_mat based on xy_coords so that they are placed
  ## in the multidimensional array in the correct arrangement
  proper_order = order(xy_coord[ , 2], xy_coord[ , 1])
  sp_mat = sp_mat[proper_order, ]
  psp = array(sp_mat, dim=c(N, M, S))
  psp = aperm(psp, c(3, 1, 2))
  spSums = apply(psp, 1, sum)
  ## drop species that never occur
  if(any(spSums %in% 0))
    psp = psp[spSums > 0, , ]
  return(psp)
}

##3.8##
getSAR = function(psp, grains, mv_window=FALSE)
{
  ## Purpose: to construct spatially explict SAR based upon a
  ## mapped grid of occurances
  ## this function replaces the older function 'grid.SAR'
  ## Arguments:
  ## psp: community array (i.e., S x N x M abundance array where N >= M)
  ## grains: the areas in pixels for which to compute the SAR
  ##         only grains that have integer log base 2 are considered
  ## mv_window: FALSE indicates that a non-moving window SAR will be calculated
  if (class(psp) != 'array')
    stop('psp must be a community array (S X N X M)')
  grains = grains[log2(grains) == round(log2(grains))]
  ## define the size of sampling units on each side
  lenN = n_pixels_long(log2(grains))
  lenM = n_pixels_wide(log2(grains))
  sr = rep(0, length(grains))
  ind = rep(0, length(grains))
  std = rep(0, length(grains))
  cs = rep(0, length(grains))
  S = dim(psp)[1]
  N = dim(psp)[2]
  M = dim(psp)[3]
  if (M > N) {
    stop('The first spatial dimension of psp must be larger than or equal to the second 
         (i.e. psp[S,N,M] where N >= M)')
  }       
  for (l in seq_along(grains)) {
    if (grains[l] == 1) {  # if area=1
      sr[l] = sum(psp > 0)
      std[l] = sd(as.vector(psp > 0))
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
      sr_vec = NULL
      for (n in brksN) {
        for (m in brksM) {
          psp_tmp = psp[ , n:(n + (lenN[l] - 1)),
                             m:(m + (lenM[l] - 1))]
          if (S == 1) {
            sr_vec = c(sr_vec, any(psp_tmp > 0) * 1)
            ind[l] = ind[l] + sum(psp_tmp)
          }
          else {
            sr_vec = c(sr_vec, sum(apply(psp_tmp > 0, 1, sum) > 0))
            ind[l] = ind[l] + sum(apply(psp_tmp, 1, sum))
          }
          cs[l] = cs[l] + 1
        }
      }
      sr[l] = sum(sr_vec)
      std[l] = sd(sr_vec)
    }
  }
  out = cbind(grains, sr / cs, ind / cs, cs, std)  
  colnames(out) = c('grains', 'richness', 'indiv', 'count', 'sr_std')
  return(out)
}  

getEAR = function(psp, grains, mv_window=FALSE)
{
  ## Purpose: to construct spatially explict Endemics Area Relationship
  ## based upon a mapped grid of occurances
  ## Arguments:
  ## psp: community array (i.e., S x N x M abundance array where N >= M)
  ## grains: the areas in pixels for which to compute the SAR
  ##         only grains that have integer log base 2 are considered
  ## mv_window: FALSE indicates that a non-moving window SAR will be calculated
  stop('Error function has not been tested!')
  if (class(psp) != 'array')
    stop('psp must be a community array (S X N X M)')
  grains = grains[log2(grains) == round(log2(grains))]
  ## define the size of sampling units on each side
  lenN = n_pixels_long(log2(grains))
  lenM = n_pixels_wide(log2(grains))
  sr = rep(0, length(grains))
  ind = rep(0, length(grains))
  std = rep(0, length(grains))
  cs = rep(0, length(grains))
  S = dim(psp)[1]
  N = dim(psp)[2]
  M = dim(psp)[3]
  if (M > N) {
    stop('The first spatial dimension of psp must be larger than or equal to the second 
         (i.e. psp[S,N,M] where N >= M)')
  }       
  for (l in seq_along(grains)) {
    if (mv_window) {
      brksN = 1:(N - lenN[l] + 1)
      brksM = 1:(M - lenM[l] + 1)
    }
    else{
      brksN = seq(1, N, lenN[l])
      brksM = seq(1, M, lenM[l])
    }  
    sr_vec = NULL
    for (n in brksN) {
      for (m in brksM) {
        n_indices = n:(n + (lenN[l] - 1))
        m_indices = m:(m + (lenM[l] - 1))
        psp_in = psp[ , n_indices, m_indices]
        psp_out = psp[ , -n_indices, -m_indices]
        
        if (S == 1) {
          sr_vec = c(sr_vec, any(psp_tmp > 0) * 1)
          ind[l] = ind[l] + sum(psp_tmp)
        }
        else {
          
          sr_vec = c(sr_vec, sum(apply(psp_tmp > 0, 1, sum) > 0))
          ind[l] = ind[l] + sum(apply(psp_tmp, 1, sum))
        }
        cs[l] = cs[l] + 1
      }
    }
    sr[l] = sum(sr_vec)
    std[l] = sd(sr_vec)
    }
  out = cbind(grains, sr / cs, ind / cs, cs, std)  
  colnames(out) = c('grains', 'richness', 'indiv', 'count', 'sr_std')
  return(out)
}  


##3.9##
aggr_comm_matrix = function(mat, coords, bisect, grains, binary=FALSE){
  ## Purpose: This function generates aggregated community matrices for each
  ## spatial grain that is specified. This function is only approrpriate for 
  ## data from a regular square spatial grid.  Coordinates, bisect, and 
  ## grains will be generated if not supplied
  ## Inputs:
  ## mat: site x species matrix assumed to be sampled from a square grid
  ## coords: two column spatial x and y coordinates
  ## bisect: specifies the levels of bisection that the community should be 
  ##   aggregated at
  ## grains: 
  ## binary: boolean, if TRUE binary pres/abse matrix returned
  S = ncol(mat)
  N = nrow(mat)
  bisect_start = log2(N)
  if (missing(bisect)) {
    bisect = log2(bisect_start) : 2
  } 
  if (any(bisect > bisect_start)) {
    stop(paste('The number of bisections must be less than ',
               bisect_start, sep=''))  
  }
  if (missing(coords)) {
    coords = expand.grid(1 : n_pixels_long(bisect_start),
                         1 : n_pixels_wide(bisect_start))
  }
  D = dist(coords)
  minD = min(D)
  domain = c(min(coords[ , 1]), max(coords[ , 1]) + minD,
             min(coords[ , 2]), max(coords[ , 2]) + minD)
  xdiff = abs(domain[1] - domain[2]) 
  ydiff = abs(domain[3] - domain[4]) 

  if (xdiff > ydiff) {
    xlengths = xdiff / n_pixels_long(bisect)
    ylengths = ydiff / n_pixels_wide(bisect)
  }
  else if (xdiff < ydiff) {
    xlengths = xdiff / n_pixels_wide(bisect)
    ylengths = ydiff / n_pixels_long(bisect)
  }
  else if (xdiff == ydiff) { ## extent is a sqr.
    xlengths = xdiff / sqrt(2^bisect)
    ylengths = ydiff / sqrt(2^bisect)
  }
  else
    stop('Function cannot figure out how to split up the area')
  n_quadrats = sum(2^bisect)
  if (missing(grains)) {
    grains = round(xlengths * ylengths, 2)
  }
  comms = matrix(NA, nrow=sum(n_quadrats), ncol=S + 3)
  if (is.null(colnames(mat))) {
    colnames(mat) = paste('sp', 1:ncol(mat), sep='')
  }  
  colnames(comms) = c('grain', 'x', 'y', colnames(mat))
  irow = 1
  for (i in seq_along(bisect)) {
    if (bisect[i] == bisect_start) {
      indices = irow : (irow + N - 1)
      comms[indices, 1:3] = cbind(grains[i], coords[ , 1], coords[ , 2])    
      comms[indices, -(1:3)] = as.matrix(mat)
      irow = max(indices) + 1
    }
    else {
      xbreaks = seq(domain[1], domain[2], xlengths[i])
      ybreaks = seq(domain[3], domain[4], ylengths[i]) 
      for (x in 1:(length(xbreaks) - 1)) {
        for (y in 1:(length(ybreaks) - 1)) {
          inQuad =  xbreaks[x] <= coords[ , 1] & coords[ , 1] < xbreaks[x + 1] & 
                    ybreaks[y] <= coords[ , 2] & coords[ , 2] < ybreaks[y + 1]
          comms[irow, 1:3] = c(grains[i], x, y)
          comms[irow, -(1:3)] = apply(mat[inQuad, ], 2, sum)
          irow = irow + 1 
        }
      }
    }  
  }
  if (binary)
    comms[ , -(1:3)] = as.numeric(comms[ , -(1:3)] > 0)
  return(comms)
}


##3.10##
varExp = function(mat) {
  # first evaulate if pres/abs or abundance matrix
  spmeans = apply(mat, 2, mean)
  if (sum(mat > 1) > 0)  # for abundance 
    sum(spmeans) 
  else                   # for occupancy
    sum(spmeans * (1 - spmeans)) 
}

##3.11##
spCommExpPoi = function(a, b, n) {
  ## returns the expected species commonality under the Poisson expectation for
  ## two samples of area a and b and species abundances n
  ## Plotkin and Muller-Lanadu 2002
  ## Eq. 10
  ## arguments:
  ## a: proportion of area of a randomly selected sample from a larger total area
  ## b: proportion of area of a randomly selected sample from a larger total area
  ## n: the species total abundances across the study site  
   sum(apply(sapply(n, function(x) 
                    sapply(c(a, b), function(y) 
                           1 - exp(-x * y))), 2, prod))
}

##3.12##
exp_S_poi = function(a, b, n) {
  ## old name: spAvgExpPoi()
  ## Plotkin and Muller-Lanadu 2002, Eq. 8
  ## Returns:
  ## the average expected number of species under the Poisson expectation for
  ## two samples of proportional area a and b and species abundances n
  ## Note this is an approximation of what Plotkin refers to as the binomial 
  ## distribution which is also known as the coleman model. These two are 
  ## very similar when a << 1
  ## Arguments:
  ## a: proportion of area of a randomly selected sample from a larger total area
  ## b: proportion of area of a randomly selected sample from a larger total area
  ## n: the species total abundances across the study site
  .5 * sum(sapply(n, function(x) 
                  sapply(c(a, b), function(y) 
                         1 - exp(-x * y))))
}

exp_S_binom = function(A, A0, n) {
  ## old name: spAvgExpColeman()
  ## Expected number of species from Coleman (1981), Eq. 1.12
  ## See Eq. 3.11 for an alternative formultaion using the SAD
  ## Arguments:
  ## A: area of a randomly selected sample from a larger total area
  ## A0: the total area of the site  
  ## n: the species total abundances across the study site
  ## Returns:
  ## the average expected number of species under the binomial distr for
  ## a sample of area A out of A0 and species abundances n
  S = sum(1 - (1 - (A/A0)) ^ n)
  return(S)
}

exp_Svar_binom = function(A, A0, n) {
  ## old name: spVarExpColeman()
  ## Coleman (1981), Eq. 1.13
  ## See Eq. 3.12 for an alternative formulation using the SAD
  ## Arguments:
  ## a: proportion of area of a randomly selected sample from a larger total area
  ## n: the species total abundances across the study site
  ## Returns:
  ## the variance in the expected number of species under the binomial 
  ## distr for a sample of proportional area a and species abundances n
  a = A / A0
  Svar = sum((1 - a) ^ n) - sum((1 - a) ^ (2 * n))
  return(Svar)
}

exp_S_logser_binom = function(A, A0, S0, N0, version='harte') {
  ## version == 'coleman' returns Coleman 1981 equ. 3.11
  ## version == 'harte' returns expectation following Harte 2011 equ 3.13
  joint_prob = rep(NA, N0)
  beta_sad = as.numeric(get_beta_sad_mle(S0, N0)[2])
  for (n in 1:N0) {
    prob_not_occur = (1 - (A/A0))^n
    prob_occur = 1 - prob_not_occur
    prob_of_n = exp(-beta_sad * n) / (n * log(beta_sad^-1))
    if (version == 'harte')
      joint_prob[n] = prob_occur * prob_of_n 
    if (version == 'coleman')
      joint_prob[n] = prob_not_occur * prob_of_n 
  }  
  if (version == 'harte')
    S = S0 * sum(joint_prob)
  if (version == 'coleman')
    S = S0 * (1 -  sum(joint_prob))
  return(S)
}

sar_logser_binom = function(Avals, A0, S0, N0, both=FALSE) {
  if (both) {
    Sexp = matrix(NA, ncol=2, nrow= length(Avals))
    Sexp[ , 1] = sapply(Avals, function(A) exp_S_logser_binom(A, A0, S0, N0, 'harte'))
    Sexp[ , 2] = sapply(Avals, function(A) exp_S_logser_binom(A, A0, S0, N0, 'coleman'))
  }
  else
    Sexp = sapply(Avals, function(A) exp_S_logser_binom(A, A0, S0, N0, 'harte'))
  return(Sexp)
}

##3.13##
spCommExpBin = function(a, b, areaTot, numOcc) {
  ## expected species commonality when pres/abse data is only available
  ## arguments:
  ## a: proportion of area of a randomly selected sample from a larger total area
  ## b: proportion of area of a randomly selected sample from a larger total area
  ## areaTot: the total area from which a and b are drawn
  ## numOcc: the number of occurances in the dataset of each species
  sum(apply(sapply(numOcc, function(x) 
                   sapply(c(a, b), function(y) 
                          1 - (1 - (x / areaTot))^(y * areaTot))), 2, prod))
}

##3.14##
spAvgExpBin = function(a, b, areaTot, numOcc) {
  ## expected species richness when pres/abse data is only available
  ## arguments:
  ## a: proportion of area of a randomly selected sample from a larger total area
  ## b: proportion of area of a randomly selected sample from a larger total area
  ## areaTot: the total area from which a and b are drawn
  ## numOcc: the number of occurances in the dataset of each species
  .5 * sum(sapply(numOcc, function(x) 
                  sapply(c(a, b), function(y) 
                         1 - (1 - (x / areaTot))^(y * areaTot))))
}

##3.15##
sorExp = function(mat, areaSampA, areaSampB=NULL){
  ## the expected binary sorensen similarity index
  if(is.null(areaSampB)) 
    areaSampB = areaSampA
  areaTot = nrow(mat)  
  a = areaSampA / areaTot
  b = areaSampB / areaTot
  #first evaulate if pres/abs or abundance matrix
  if (sum(mat > 1) > 0) { #  for abundance use Poisson expectation
    n = apply(mat, 2, sum)
    sor = spCommExpPoi(a, b, n) / exp_S_poi(a, b, n)
  }
  else{
    numOcc = apply(mat, 2, sum)
    sor = spCommExpBin(a, b, areaTot, numOcc) / spAvgExpBin(a, b, areaTot, numOcc)
  }
  return(sor)
}

##3.16##
jacExp = function(mat, areaSampA, areaSampB=NULL){
  ## the expected binary jaccard similarity index
  if(is.null(areaSampB)) 
    areaSampB = areaSampA
  areaTot = nrow(mat)  
  a = areaSampA / areaTot
  b = areaSampB / areaTot
  #first evaulate if pres/abs or abundance matrix
  if (sum(mat > 1) > 0) { #  for abundance use Poisson expectation
    n = apply(mat, 2, sum)
    jac = spCommExpPoi(a, b ,n) / (2 * exp_S_poi(a, b, n) - spCommExpPoi(a, b, n))
  }
  else{
    numOcc = apply(mat, 2, sum)
    areaTot = nrow(mat)
    jac = spCommExpBin(a, b, areaTot, numOcc) / 
          (2 * spAvgExpBin(a, b, areaTot, numOcc) - spCommExpBin(a, b, areaTot, numOcc))

  }
  return(jac)
}

##3.17##
calc_metrics = function(comms, metricsToCalc, dataType, grain=1, breaks=NA, log=FALSE,
                        hmin=NA, hmax=NA, quants=NA, direction='omnidirectional',
                        tolerance=NA, nperm=NULL, npar=1, RPargs=NULL, univariate=FALSE,
                        writeToFile=FALSE,fileSuffix=NULL) {
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
  ## univariate: if TRUE then results are computed on a per species basis
  ## writeToFile: if True an .Rdata file is written for each metric calculated
  ## fileSuffix: add a file identifying string here
  if(metricsToCalc == 'all')
    metricsToCalc = c('varWithin', 'varBetween', 'jaccard', 'sorensen')
  if(writeToFile) {
    if(direction != 'omnidirectional')
      fileSuffix = paste(fileSuffix,'_', direction, 'deg', sep='') 
  }
  grains = unique(comms[,1])
  out = vector('list', length(grains))
  names(out) = paste('comm', grains,sep='')
  for (i in seq_along(grains)) {
    true = comms[ , 1] == grains[i]
    ## because we only consider grains that are perfect squares 
    coords = as.matrix(comms[true, 2:3]) * sqrt(grains[i])
    mat = as.matrix(comms[true, -c(1:3)])
    if (!is.na(breaks[i]))
      brks = breaks[i]
    else
      brks = NA
    if(dataType == 'binary')
      mat = (mat > 0) * 1
    if(any('varWithin' %in% metricsToCalc)){
      if(i == 1){
        varWithin = vector('list', length(grains))      
        names(varWithin) = grains
      }  
      varWithinObs = vario(mat,coords,grain,brks,log,hmin_val,hmax,pos.neg=FALSE,
                           quants=quants,direction=direction,tolerance=tolerance,
                           unit.angle='degrees',univariate=univariate)
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
    if (any('varBetween' %in% metricsToCalc)) {
      if (i == 1) {
        varBetween = vector('list', length(grains))      
        names(varBetween) = grains
      }  
      varBetweenObs = vario(mat,coords,grain,brks,log,hmin,hmax,pos.neg=TRUE,
                            quants=quants,direction=direction,tolerance=tolerance,
                            unit.angle='degrees',univariate=univariate) 
      if (!is.null(nperm)) { 
        varBetweenNull = null.perms(mat,varBetweenObs,nperm,coords=coords,
                                    meth='randpat',RPargs=RPargs,npar=npar)
      }
      else {
        varBetweenNull = NULL
      }   
      varBetween[[i]] = list(varBetweenObs=varBetweenObs,
                             varBetweenNull=varBetweenNull)
    }
    if (any('jaccard' %in% metricsToCalc)) {
      if (i == 1) {
        jaccard = vector('list', length(grains))  
        names(jaccard) = grains
      }  
      jaccardObs  = vario(mat,coords,grain,brks,log,hmin,hmax,distance.metric='jaccard',
                          quants=quants, direction=direction,tolerance=tolerance,
                          unit.angle='degrees',univariate=univariate) 
      jaccardNull = NULL
      if (!is.null(nperm)) {
        jaccardNull = null.perms(mat, jaccardObs, nperm, coords=coords,
                                 meth='random', npar=npar, breaks=brks)
      }        
      jaccardExp = NULL
      if(dataType == 'binary') {
        jaccardExp = 1 - jacExp(mat, 1) #  to convert into a dissimiarlity
      } 
      jaccard[[i]] = list(jaccardObs=jaccardObs,jaccardNull=jaccardNull,
                          jaccardExp=jaccardExp)
    }
    if (any('sorensen' %in% metricsToCalc)) {
      if (i == 1) {
        sorensen = vector('list', length(grains))
        names(sorensen) = grains
      }  
      ## bray-curtis is equiv to sorensen        
      sorensenObs  = vario(mat,coords,grain,brks,log,hmin,hmax,distance.metric='bray',
                           quants=quants,direction=direction,tolerance=tolerance,
                           unit.angle='degrees',univariate=univariate) 
      sorensenNull = NULL
      if (!is.null(nperm)) {
          sorensenNull = null.perms(mat, sorensenObs, nperm, coords=coords,
                                    meth='random', npar=npar, breaks=brks)
      }    
      sorensenExp = NULL
      if (dataType == 'binary') {
        sorensenExp = 1 - sorExp(mat,1) #  to convert into a dissimiarlity
      } 
      sorensen[[i]] = list(sorensenObs=sorensenObs, sorensenNull=sorensenNull,
                           sorensenExp=sorensenExp)
    }
    out[[i]] = list()
    for (j in metricsToCalc) { 
      out[[i]][[j]] = eval(parse(text=paste(j,'[[',i,']]')))
    }
    if (writeToFile) {      
      ## update result files as loop proceeds
      if (any('varWithin' %in% metricsToCalc)) {
        save(varWithin, file=paste('./varWithin/varWithin_',
             fileSuffix, '_', dataType, '.Rdata', sep=''))  
      }
      if (any('varBetween' %in% metricsToCalc)) {
        save(varBetween, file=paste('./varBetween/varBetween_',
             fileSuffix, '_', dataType, '.Rdata', sep=''))
      }
      if (any('jaccard' %in% metricsToCalc)) {
        save(jaccard, file=paste('./jaccard/jaccard_',
             fileSuffix, '_', dataType, '.Rdata', sep=''))
      }
      if (any('sorensen' %in% metricsToCalc)) {
        save(sorensen, file=paste('./sorensen/sorensen_',
             fileSuffix, '_', dataType, '.Rdata', sep=''))
      }             
    }
  }
  return(out)
}


calc_metrics_par = function(comms,metricsToCalc,dataType,npar,grain=1,breaks=NA,
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

calc_metrics_bisect = function(comms, metricsToCalc, dataType, swap, 
                               quants=NA, nperm=NA, univariate=FALSE, 
                               writeToFile=FALSE, fileSuffix=NULL)
{
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
  ## swap: two options: indiv or sample for individual or sample-based shuffling
  ## quants: the quantiles to compute
  ## nperm: number of permutations to carry out for null models
  ## univariate: if TRUE then results are computed on a per species basis
  ## writeToFile: if True an .Rdata file is written for each metric calculated
  ## fileSuffix: add a file identifying string here
  if(metricsToCalc == 'all')
    metricsToCalc = c('varWithin', 'varBetween', 'jaccard', 'sorensen')
  if(writeToFile) {
    if (univariate)
      fileSuffix = paste(fileSuffix,'_uni_bisect', sep='') 
    else
      fileSuffix = paste(fileSuffix,'_bisect', sep='') 
  }
  grains = unique(comms[,1])
  out = vector('list', length(grains))
  names(out) = paste('comm', grains,sep='')
  for (i in seq_along(grains)) {
    true = comms[ , 1] == grains[i]
    ## because we only consider grains that are perfect squares 
    coords = as.matrix(comms[true, 2:3]) * sqrt(grains[i])
    mat = as.matrix(comms[true, -c(1:3)])
    if(dataType == 'binary')
      mat = (mat > 0) * 1
    if(any('varWithin' %in% metricsToCalc)){
      if(i == 1){
        varWithin = vector('list', length(grains))      
        names(varWithin) = grains
      }  
      varWithinObs = vario_bisect(mat,coords,quants=quants,univariate=univariate)
      if(!is.na(nperm)){ 
        varWithinNull = random_shuffle(mat,varWithinObs, swap, nperm,coords)
      }
      else{
        varWithinNull = NULL
      }       
      varWithinExp = varExp(mat)
      varWithin[[i]] = list(varWithinObs=varWithinObs,varWithinNull=varWithinNull,
                            varWithinExp=varWithinExp)
    }
    if (any('varBetween' %in% metricsToCalc)) {
      if (i == 1) {
        varBetween = vector('list', length(grains))      
        names(varBetween) = grains
      }  
      varBetweenObs = vario_bisect(mat,coords,quants=quants,univariate=univariate) 
      if (!is.na(nperm)) { 
        varBetweenNull = random_shuffle(mat,varBetweenObs,swap,nperm,coords)
      }
      else {
        varBetweenNull = NULL
      }   
      varBetween[[i]] = list(varBetweenObs=varBetweenObs,
                             varBetweenNull=varBetweenNull)
    }
    if (any('jaccard' %in% metricsToCalc)) {
      if (i == 1) {
        jaccard = vector('list', length(grains))  
        names(jaccard) = grains
      }  
      jaccardObs  = vario_bisect(mat,coords, distance.metric='jaccard',
                          quants=quants, univariate=univariate) 
      jaccardNull = NULL
      if (!is.na(nperm)) {
        jaccardNull = random_shuffle(mat, jaccardObs, swap, nperm, coords)
      }        
      jaccardExp = NULL
      if(dataType == 'binary') {
        jaccardExp = 1 - jacExp(mat, 1) #  to convert into a dissimiarlity
      } 
      jaccard[[i]] = list(jaccardObs=jaccardObs,jaccardNull=jaccardNull,
                          jaccardExp=jaccardExp)
    }
    if (any('sorensen' %in% metricsToCalc)) {
      if (i == 1) {
        sorensen = vector('list', length(grains))
        names(sorensen) = grains
      }  
      ## bray-curtis is equiv to sorensen        
      sorensenObs  = vario_bisect(mat,coords, distance.metric='bray',
                                  quants=quants, univariate=univariate) 
      sorensenNull = NULL
      if (!is.na(nperm)) {
        sorensenNull = random_shuffle(mat, sorensenObs, swap, nperm, coords)
      }    
      sorensenExp = NULL
      if (dataType == 'binary') {
        sorensenExp = 1 - sorExp(mat,1) #  to convert into a dissimiarlity
      } 
      sorensen[[i]] = list(sorensenObs=sorensenObs, sorensenNull=sorensenNull,
                           sorensenExp=sorensenExp)
    }
    out[[i]] = list()
    for (j in metricsToCalc) { 
      out[[i]][[j]] = eval(parse(text=paste(j,'[[',i,']]')))
    }
    if (writeToFile) {      
      ## update result files as loop proceeds
      if (any('varWithin' %in% metricsToCalc)) {
        save(varWithin, file=paste('./varWithin/varWithin_',
                                   fileSuffix, '_', dataType, '.Rdata', sep=''))  
      }
      if (any('varBetween' %in% metricsToCalc)) {
        save(varBetween, file=paste('./varBetween/varBetween_',
                                    fileSuffix, '_', dataType, '.Rdata', sep=''))
      }
      if (any('jaccard' %in% metricsToCalc)) {
        save(jaccard, file=paste('./jaccard/jaccard_',
                                 fileSuffix, '_', dataType, '.Rdata', sep=''))
      }
      if (any('sorensen' %in% metricsToCalc)) {
        save(sorensen, file=paste('./sorensen/sorensen_',
                                  fileSuffix, '_', dataType, '.Rdata', sep=''))
      }             
    }
  }
  return(out)
}



n_pixels_long = function(i_bisect){
  ## returns the number of pixels on the side of a grid with more or equal pixels
  ## after i bisection events
  ## old function name: len
  2^floor((i_bisect + 1) / 2) 
}

n_pixels_wide = function(i_bisect){
  ## returns the number of pixels on the side of a grid with less or equal pixels
  ## after i bisecttion events
  ## old function name: wid
  2^floor(i_bisect / 2) 
} 

##3.19 ##
make_comm_matrix = function(spnum, S, coords, n_quadrats, domain, abu = NULL,
                              grainSuffix=NULL)
{ 
  ## Output: 
  ## A community matrix where each row is a differnet pixel on a grid.  
  ## Arguments:
  ## spnum : an integer specifying species identities
  ## S : the size of the species pool may be larger than the number of unique 
  ##     spnum
  ## coords : two column matrix (x,y) specifying the spatial coordinates of each stem
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
    xlengths = ylengths = rep(NA, length(n_quadrats))
    for (i in seq_along(n_quadrats)) {
      if (log2(n_quadrats[i]) %% 2 == 1) {
        ## if # of bisections is odd then arbitrarily 
        ## make x dimension have more cells
        xlengths[i] = xdiff / n_pixels_long(log2(n_quadrats[i]))
        ylengths[i] = ydiff / n_pixels_wide(log2(n_quadrats[i]))
      }  
      else {
        xlengths[i] = xdiff / sqrt(n_quadrats[i])
        ylengths[i] = ydiff / sqrt(n_quadrats[i])
      } 
    }  
  }
  else
    stop('Function cannot figure out how to split up the area')
  comms = matrix(NA, nrow=sum(n_quadrats), ncol=S + 3)
  colnames(comms) = c('grain', 'x', 'y', paste('sp', 1:S, sep=''))
  irow = 1
  for (i in seq_along(n_quadrats)) {
    xbreaks = seq(domain[1], domain[2], xlengths[i])
    ybreaks = seq(domain[3], domain[4], ylengths[i]) 
    for (x in 1:(length(xbreaks) - 1)) {
      for (y in 1:(length(ybreaks) - 1)) {
        inQuad =  xbreaks[x] <= coords[ , 1] & coords[ , 1] < xbreaks[x + 1] & 
                  ybreaks[y] <= coords[ , 2] & coords[ , 2] < ybreaks[y + 1]
        if (is.null(grainSuffix)) {
          comms[irow, c(1:3)] = c(paste(round(xlengths[i] * ylengths[i], 2), sep=''),
                                  x, y)
        }
        else {
          comms[irow, c(1:3)] = c(paste(round(xlengths[i] * ylengths[i], 2),
                                  grainSuffix, sep=''), x, y)
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

get_bisection_history = function(grains, abu) {
  ## retrives the record of how a species' abunance was bisected in a particular
  ## landscape the input into this function can be generated by the function 
  ## 'make_comm_matrix'
  ## Note: this function is necessary for mle fitting of Conlisk et al. models
  ## Note: does not seem to work properly when shape is rectangular
  ## arguments:
  ## grains: a vector of spatial grains that is associated with each abundance
  ## abu: a vector of abundances across the spatial scales of interest
  abu = as.numeric(abu)
  grains = as.numeric(grains)
  grains_table = table(grains)
  grains_unique = as.integer(names(grains_table))
  i_bisect = log2(as.integer(grains_table))
  out = data.frame(n1 = NULL, n2 = NULL)
  for (i in seq_along(grains_unique)) {
    abu_tmp = abu[grains == grains_unique[i]]
    if (i_bisect[i] == 1) {
      n1 = abu_tmp[1] 
      n2 = abu_tmp[2]
    }  
    else {
      if (i_bisect[i] %% 2 == 0) { 
        true = 1:length(abu_tmp) %% 2 == 1
        n1 = abu_tmp[true]
        n2 = abu_tmp[!true]          
      }
      else {
        block_size = 2^(i_bisect[i] - ceiling(i_bisect[i] / 2))
        true = rep(rep(c(T,F), each=block_size), times=block_size)
        n1 = abu_tmp[true]
        n2 = abu_tmp[!true]
      }
    }
    out = rbind(out,data.frame(n1 = n1, n2 = n2))
  }              
  return(out)
}

paste_list = function(a_list) {
  ## this function takes each list element and pastes it to the 
  ## other list elements sequentially
  if (length(a_list) > 1) # if not at the end of the list 
    out = paste(a_list[[1]], paste_list(a_list[-1]), sep='')
  else  # you are at the end of the list
    out = paste(a_list[[1]], sep='')
  return(out)
}

get_bisect_matrix = function(bisect_string) {
  ## splits a vector of 0/1 strings into a matrix where
  ## each digit is in its own column. This effectively converts
  ## the string representation of a bisection coordinate into a 
  ## numerical matrix.
  ## E.g. '010' becomes 0 1 0 
  ## Arguments
  ## bisect_coords: a string representation of the bisection coordinates
  ## in which 0 represents a left bisection and 1 reprents a right bisection
  if (class(bisect_string) != 'character')
    stop('bisect_string must be a character string')
  i_bisect = nchar(bisect_string[1])
  bisect_matrix = matrix(NA, ncol=i_bisect, nrow=length(bisect_string))
  colnames(bisect_matrix) = paste('i=', 1:i_bisect, sep='')
  for (i in 1:i_bisect) {
    bisect_matrix[ , i] = as.numeric(substr(bisect_string, i, i))
  }
  return(bisect_matrix)
}

get_bisect_coords = function(i_bisect) {
  ## generates a data.frame of coordinates that indictates what pattern 
  ## of bisections (0 = left bisection or 1 = right bisect)
  ## that resulted in a given quadrat. The function also
  ## returns a set of traditional xy-coordinates
  ## Arguments:
  ## i_bisect: the number of bisections
  ## Returns:
  ## a matrix with i_bisect + 2 columns 
  ## x, y, and the bisection coordinates each as their own column
  ## for the bisection coordinates, 0 indicates the 
  ## quadrat was on the left of the bisection, and 1 indicates
  ## the quadrat was on the right of the bisection.
  ## Example: 0 1 1 indicates that 1st bisection was left, and 
  ## 2nd and 3rd bisections where right
  ## Note:
  ## this function assumes that the long side of rectangular
  ## shaped A0 is along the x-axis (i.e., the 1st bisection is 
  ## along the x-axis)
  if (i_bisect < 0 )
    stop('i_bisect must be greater than zero')
  if (i_bisect != round(i_bisect))
    stop('i_bisect must be an integer')
  mat_list = vector('list', length=i_bisect)
  n_cols = n_pixels_long(i_bisect)
  n_rows = n_pixels_wide(i_bisect)
  flag = TRUE
  for (j in i_bisect:1) {
    fill = rep(c(0,1), each= n_pixels_wide(j - 1))
    mat_list[[i_bisect - j + 1]] = matrix(fill, byrow=flag, 
                                          ncol=n_cols, nrow=n_rows)
    flag = !flag
  }
  bisect_string = paste_list(mat_list)
  bisect_matrix = get_bisect_matrix(bisect_string)
  xy_coords = expand.grid(y = n_pixels_wide(i_bisect):1, 
                          x = 1:n_pixels_long(i_bisect))[ , 2:1]
  bisect_coords = cbind(xy_coords, bisect_matrix)
  return(bisect_coords)
}

get_pairs_matrix = function(bisect_coords, j) {
  ## Returns a binary matrix that indicates which samples are to be 
  ## compared for a given seperation order (j). 
  ## Arguments
  ## bisect_coords: a data.frame that contains the x, y, and 
  ##   bisection coordinates, the result of the function get_bisect_coords()
  ## j: seperation order, 1 <= j <= i_bisect
  if (j < 1)
    stop('j must be greater than 1')
  if (j > (ncol(bisect_coords) - 2))
    stop('j must be less than or equal to the number of bisections')
  bisect_mat = bisect_coords[ , -(1:2)]
  if (j == 1) {
    pairs_mat = as.matrix(dist(bisect_mat[ , j])) == 1 
  }
  else {
    same_bisect = as.matrix(dist(bisect_mat[ , 1:(j - 1)])) == 0
    diff_bisect = as.matrix(dist(bisect_mat[ , j])) == 1
    pairs_mat  = same_bisect * diff_bisect
  }
  pairs_lower_tri = pairs_mat  * lower.tri(pairs_mat)
  return(pairs_lower_tri)
}

load_get_sep_dist = function(c_path='./scripts'){
  if (!is.loaded('get_sep_dist')) {
    OS = Sys.info()['sysname']
    if (OS == 'Linux')
      dyn.load(file.path(c_path, 'get_sep_dist.so'))
    else
      dyn.load(file.path(c_path, 'get_sep_dist.dll'))
  }
}

dist_bisect = function(i_bisect, use_c=FALSE, c_path='./scripts') {
  ## Returns a list with two elements
  ## dist: a lower-triangular distance matrix (class dist)
  ## crd: the spatial euclidean coordinates of 
  ## that is populated with the bisection seperation order of
  ## each unique pairwise comparison
  N = 2^i_bisect
  coords = get_bisect_coords(i_bisect)
  coords = coords[order(coords[ , 2], coords[ , 1]), ]
  vec_coords = as.vector(apply(coords[ , -(1:2)], 1, as.vector))
  sep_dist = rep(0, (N * (N - 1)) / 2)
  if (use_c) {
    load_get_sep_dist(c_path)
    sep_dist = .C("get_sep_dist", vec_coords = as.integer(vec_coords),
                  i_bisect = as.integer(i_bisect), N = as.integer(N), 
                  sep_dist = as.integer(sep_dist))$sep_dist
  }
  else {
    icount = 1
    for(i in 1:(N-1)) {
      for(j in (i + 1):N) {
        diff = 0
        sep = 0
        while (diff == 0) {
          sep = sep + 1
          diff = vec_coords[(i - 1) * i_bisect + sep] -
                 vec_coords[(j - 1) * i_bisect + sep]
        }
        sep_dist[icount] = sep
        icount = icount + 1
      }
    }
  }  
  place_matrix = matrix(1:N^2, N, N)
  places = as.vector(as.dist(place_matrix))
  sep_dist_mat = matrix(NA, N, N)
  sep_dist_mat[places] = sep_dist
  sep_dist = as.dist(sep_dist_mat)
  crd = coords[ , 1:2]
  rownames(crd) = 1:nrow(crd)
  out = list(dist = sep_dist, crd = crd)
  return(out)
}

mete_sor_transform = function(dat) {
  ## this effectively creates a new distance matrix populated with 
  ## sorensen predictions 
  bisections = unique(dat$i)
  out = vector('list', length(bisections))
  names(out) = paste('sor, i=', bisections, sep='')
  icount = 1
  for (i_bisect in bisections) {
    dat_tmp = dat[dat$i == i_bisect, ]
    sep_dist = as.matrix(dist_bisect(i_bisect, use_c=TRUE)$dist)
    n = nrow(sep_dist)
    sor_vect = dat_tmp$sor[match(as.vector(sep_dist), dat_tmp$j)]
    out[[icount]] = as.dist(matrix(sor_vect, n, n))
    icount = icount + 1
  }
  return(out)
}

dimensional_match = function(coord1, coord2) {
  ## Returns boolean indicating whether or not the spatial dimensions of 
  ## two sets of coordinates match each other. 
  ## Note: this does not check that the exact values match, but simply that
  ## the number of unique x and y values in the two coordinate sets match
  ## Note: coordinates of an arbitrary dimension can be checked
  ## Arguments
  ## coord1: the first set of coordinates
  ## coord2: the second set of coordinates to be checked against.
  flag = TRUE
  icol = 1
  while (flag & icol <= ncol(coord1)) {
    uni_coord1 = length(unique(coord1[ , icol]))  
    uni_coord2 = length(unique(coord2[ , icol]))
    if (uni_coord1 != uni_coord2)
      flag = FALSE
    icol = icol + 1
  }
  return(flag)
}

vario_bisect = function(x, coord, sep_orders=NULL, distance.metric='euclidean',
                        quants=NA, univariate=FALSE){
  ## Computes a variogram based upon the possible seperation orders
  ## Only seperation orders that compare square samples are considered.
  ## Arguments:
  ## x: a sitexsp matrix
  ## coord: the spatial coordinates
  ## sep_orders: default is NULL, the serpation orders of interest, accepts
  ##   'all' or a numeric vector of specific seperation orders to examine.
  ## distance.metric': can be one of the speices turnover metrics listed by the
  ##   vegan function vegdist(). Common options include, 'jaccard' and 'bray'.
  ##   If computed on pres/abse data then soreson index is computed by 'bray'.
  if (distance.metric != 'euclidean') {
    require(vegan)
  }
  if (is.vector(x)) 
    x = matrix(x)
  else if (!is.matrix(x))
    stop('x must be a vector or a matrix')
  ## make sure the x-coordinate is the long side of the plot
  if (length(unique(coord[ , 1])) < length(unique(coord[ , 2]))) {
    ## assuming x and y coordinates are measured in the same units
    ## transform x -> y and y -> x so that x is long end of plot
    tmp_coord = coord
    coord[ , 1] = tmp_coord[ , 2]
    coord[ , 2] = tmp_coord[ , 1]
    rm(tmp_coord)
  }
  if (univariate) {
    vobject = vario_uni(x, bisect=TRUE, coord, sep_orders, distance.metric, quants)
    class(vobject) = 'vario'
  } 
  else {
    S = ncol(x)
    N = nrow(x)
    i_bisect = log2(N)
    if (i_bisect != round(i_bisect))
      stop('Number of samples must be consistent with a bisection procedure')
    if (is.null(sep_orders)) {
      if (i_bisect %% 2 == 0)
        sep_orders = seq(i_bisect, 2, -2)
      else
        sep_orders = seq(i_bisect, 1, -2)
    }
    else if (sep_orders[1] == 'all')
      sep_orders = floor(i_bisect) : 1
    else if (is.numeric(sep_orders[1]))
      sep_orders = sort(sep_orders, dec=T)
    else 
      stop('sep_orders incorrectly specified.
            Leave blank for geometry preserving comparisons,
            set to "all" to compute all possible bisections, or
            set as a numeric vector')
    bisect_coords = get_bisect_coords(i_bisect)
    ## check that spatial dimensions of coord and bisect_coords are identical
    if (!dimensional_match(coord, bisect_coords[ , 1:2]))
      stop('The spatial dimension of the community matrix does not match that expected by the function get_bisect_coords()')
    ## arrange x, coord and bisect_mat so that they are in the same order
    sort_order = order(coord[ , 2], coord[ , 1])
    x = x[sort_order, ]
    if (is.vector(x))
      x = matrix(x)
    coord = coord[sort_order, ]
    sort_order = order(bisect_coords[ , 2], bisect_coords[ , 1])
    bisect_coords = bisect_coords[sort_order, ]
    geo_dist = dist(coord)
    geo_dist = as.matrix(geo_dist)
    if (distance.metric == 'euclidean')
      sp_dist = as.matrix(dist(x))
    else
      sp_dist = as.matrix(vegdist(x, method=distance.metric))
    geo_dist_avg = rep(NA, length(sep_orders))
    sp_dist_avg = rep(NA, length(sep_orders))
    n_pairs = rep(NA, length(sep_orders)) 
    if (!is.na(quants[1])) {
      sp_dist_qt = matrix(NA, ncol = length(quants), nrow=length(sep_orders))
      colnames(sp_dist_qt) = paste(quants * 100)
    }
    for (j in seq_along(sep_orders)) {
      pairs_mat = get_pairs_matrix(bisect_coords, sep_orders[j])
      n_pairs[j] = sum(pairs_mat)
      geo_dist_mat = geo_dist * pairs_mat
      geo_dist_avg[j] = mean(geo_dist_mat[geo_dist_mat > 0])
      sp_dist_mat = sp_dist * pairs_mat
      sp_dist_avg[j] = mean(sp_dist_mat[geo_dist_mat > 0], na.rm=T)
      if (!is.na(quants[1])) {
        sp_dist_qt[j, ] = quantile(sp_dist_mat[geo_dist_mat > 0], quants, na.rm=TRUE)
      }
    }
    vobject = list() 
    class(vobject) = 'vario'
    direction = 'bisection'
    vobject$parms = data.frame(grain=1, hmin=NA, hmax=NA, S=S, N=N,
                               pos.neg=FALSE, median=FALSE, direction,
                               tolerance=NA, unit.angle=NA, distance.metric, 
                               quants = ifelse(is.na(quants[1]), NA, 
                                               paste(quants* 100, collapse=", ")))
    vobject$vario = data.frame(j = sep_orders, Dist = geo_dist_avg, n = n_pairs,
                               exp.var = sp_dist_avg, obs.var = NA)
    if (!is.na(quants[1]))
      vobject$vario = cbind(vobject$vario, var.qt = sp_dist_qt)
    if (is.vector(x))
      vobject$p = sum(x) / length(x)
    else
      vobject$p = apply(x, 2, sum, na.rm=TRUE) / nrow(x)  
    vobject$perm = FALSE
  }
  return(vobject)
}

getResults = function(names, metric, dataType, bisect=FALSE, sim_result=FALSE)
{
  ## Purpose: to import the results of the 'calc_metrics' function
  ## and to load them into a list
  ## Arguments:
  ## names: the short names that were used in the naming of the output files
  ## metric: the metric calculated, see 'calc_metrics' for options
  ## dataType: the data type of interest, 'binary' or 'abu'
  ## sim_result: specifies if the object is a simulation result
  results = vector('list', length=length(names))
  names(results) = names
  for (i in seq_along(results)) {
    if (bisect) {
      load(paste('./', metric, '/' ,metric, '_', names[i], '_bisect_', dataType, '.Rdata',
                 sep=''))
    }
    else {
      load(paste('./', metric, '/' ,metric, '_', names[i], '_', dataType, '.Rdata',
                 sep=''))
    }
    if (sim_result)
      results[[i]] = metrics
    else
      results[[i]] = eval(parse(text=metric))
  }
  return(results)
}

reshapeResults = function(results, metric, null=TRUE, perm_null=FALSE, bisect=FALSE,
                          sim_result=FALSE){
  ## Purpose: to reshape the results from a nested list to a matrix
  ## of the most important information
  ## Arguements:
  ## results: a list of results loaded by getResults
  ## metric: which metric were the results computed for
  ## perm_null: boolean, TRUE means that permutation null models were 
  ##            carried out
  ## sim_result: boolean, TRUE means that the results are from simulated
  ##            communities
  out = vector('list', length=length(results))
  names(out) = names(results)
  for (i in seq_along(results)) { ## the datset
    if (is.null(results[[i]]))
      next
    vExp = vector('list', length=length(results[[i]]))
    Dist = vExp
    n = vExp
    if (null)
      rpExp = vExp 
    for (j in seq_along(results[[i]])) { ## the grain/community
      if (is.null(results[[i]][[j]]))
        next
      if (sim_result) {
        vobject = results[[i]][[j]][[1]][[1]]$vario
        quants = !is.na(results[[i]][[j]][[1]][[1]]$parms$quants)
        if (null) {
          if (perm_null) {
            rpExp[[j]] = t(rbind(apply(results[[i]][[j]][[1]][[2]]$vario[ , 1, -1], 1, mean),
                                 apply(results[[i]][[j]][[1]][[2]]$vario[ , 1, -1], 1, quantile,
                                       c(.25, .5, .75))))
            colnames(rpExp[[j]]) = c('avg', '25', '50', '75')
          }  
          else 
            rpExp[[j]] =  data.frame(Exp.avg = rep(results[[i]][[j]][[1]][[3]], nrow(vobject)))
        }  
      }
      else {
        vobject = results[[i]][[j]][[1]]$vario
        quants = !is.na(results[[i]][[j]][[1]]$parms$quants)
        if (null) {
          if (perm_null) {
            if (bisect)
              rpExp[[j]] = t(rbind(apply(results[[i]][[j]][[2]]$vario[ , -1], 1, mean),
                                   apply(results[[i]][[j]][[2]]$vario[ , -1], 1, quantile,
                                         c(.25, .5, .75))))
            else
              rpExp[[j]] = t(rbind(apply(results[[i]][[j]][[2]]$vario[ , 1, -1], 1, mean),
                                   apply(results[[i]][[j]][[2]]$vario[ , 1, -1], 1, quantile,
                                         c(.25, .5, .75))))
            colnames(rpExp[[j]]) = c('avg', '25', '50', '75')
          }
          else
            rpExp[[j]] = data.frame(Exp.avg = rep(results[[i]][[j]][[3]], nrow(vobject)))
        }  
      }
      Dist[[j]] = vobject$Dist 
      n[[j]] = vobject$n
      vExp[[j]] = data.frame(avg = vobject$exp.var)
      if (quants) {
        if (bisect)
          qt.names = names(vobject)[grep('var.qt',names(vobject))]
        else
          qt.names = names(vobject)[grep('exp.qt',names(vobject))]
        vExp[[j]] = cbind(vExp[[j]], vobject[ , qt.names])
      }
      if (any(metric %in% c('sorensen','jaccard'))) {
        ## convert the values into similarity measures
        vExp[[j]] = 1 - vExp[[j]]
        if (null)
          rpExp[[j]] = 1 - rpExp[[j]]
      }  
      if (j == 1) {
        vExp_dat = vExp[[j]]
        if (null)
          rpExp_dat = rpExp[[j]]
      }  
      else {
        vExp_dat = rbind(vExp_dat, vExp[[j]])
        if (null)
          rpExp_dat = rbind(rpExp_dat, rpExp[[j]])
      }  
    }
    if (is.null(names(results[[i]])))
      commNames = 1:length(results[[i]])
    else
      commNames = names(results[[i]])
    if (bisect)
      names(vExp_dat) = sub('var.qt.','', names(vExp_dat))
    else
      names(vExp_dat) = sub('exp.qt.','', names(vExp_dat))
    if (any(metric %in% c('sorensen','jaccard'))) {
      names(vExp_dat[ , -1]) = rev(names(vExp_dat[ , -1]))
      if (null)
        names(rpExp_dat[ , -1]) = rev(names(rpExp_dat[ , -1]))
    }  
    if (null)
      out[[i]] = data.frame(Dist = unlist(Dist), Metric = vExp_dat,
                            Exp = rpExp_dat, N = unlist(n))
    else
      out[[i]] = data.frame(Dist = unlist(Dist), Metric = vExp_dat,
                            N = unlist(n))
    out[[i]] = data.frame(out[[i]],
               Comm = unlist(mapply(rep,commNames,each=sapply(n, length),
                                       SIMPLIFY=FALSE)))
  }
  return(out)
}  

avgResults = function(results, null=TRUE) {
  ## convert to flat file for averaging
  N = length(results)
  for (i in seq_along(results)) {
    if (i == 1){
      flat = results[[i]] 
    }
    else
      flat = rbind(flat, results[[i]])
  }
  uniComms = unique(flat$Comm)
  for (i in seq_along(uniComms)) {
    true = flat$Comm == uniComms[i]
    if (null) {
      avgs = apply(flat[true, c('Dist', 'N', 'Metric.50', 'Metric.avg', 'Exp.avg')],
                   2, function(x) tapply(x, round(flat$Dist[true], 4), mean))
      qts_lo = apply(flat[true, c('Metric.50', 'Metric.avg', 'Exp.avg')], 2, function(x)
                     tapply(x, round(flat$Dist[true], 4), quantile, 0.025))
      qts_hi = apply(flat[true, c('Metric.50', 'Metric.avg', 'Exp.avg')], 2, function(x)
                     tapply(x, round(flat$Dist[true], 4), quantile, 0.975))
      dat = data.frame(grain = as.numeric(sub('comm','',uniComms[i])),
                       Dist = avgs[,1], N = avgs[,2], 
                       med.lo = qts_lo[,1], med = avgs[,3], med.hi = qts_hi[,1], 
                       avg.lo = qts_lo[,2], avg = avgs[,4], avg.hi = qts_hi[,2],
                       exp.lo = qts_lo[,3], exp = avgs[,5], exp.hi = qts_hi[,3])
    }
    else {
      avgs = apply(flat[true, c('Dist', 'N', 'Metric.50', 'Metric.avg')],
                   2, function(x) tapply(x, round(flat$Dist[true], 4), mean))
      qts_lo = apply(flat[true, c('Metric.50', 'Metric.avg')], 2, function(x)
                     tapply(x, round(flat$Dist[true], 4), quantile, 0.025))
      qts_hi = apply(flat[true, c('Metric.50', 'Metric.avg')], 2, function(x)
                     tapply(x, round(flat$Dist[true], 4), quantile, 0.975))
      dat = data.frame(grain = as.numeric(sub('comm','',uniComms[i])),
                       Dist = avgs[,1], N = avgs[,2], 
                       med.lo = qts_lo[,1], med = avgs[,3], med.hi = qts_hi[,1], 
                       avg.lo = qts_lo[,2], avg = avgs[,4], avg.hi = qts_hi[,2])
    }
    if (i == 1)
      out = dat
    else
      out = rbind(out, dat)
  }
  return(out)
}

merge_drop = function(results) {
  result_names = names(results)
  to_drop = unlist(sapply(c('cocoli2','sherman2','sherman3'),
                          function(x) grep(x, result_names)))
  cocoli_results = (subset(results$cocoli1, select= -Comm) + 
                    subset(results$cocoli2, select= -Comm)) / 2
  cocoli_results = data.frame(cocoli_results, Comm=results$cocoli1[ , 'Comm'])
  sherman_results = (subset(results$sherman1, select= -Comm) + 
                     subset(results$sherman2, select= -Comm)) / 2
  sherman_results = data.frame(sherman_results, Comm=results$sherman1[ , 'Comm'])
  results = results[-to_drop]
  results[[grep('cocoli1', names(results))]] = cocoli_results
  results[[grep('sherman1', names(results))]] = sherman_results
  return(results)
}

get_file_names = function(path, strings, invert=FALSE) {
  indices = NULL
  for (i in seq_along(strings)) {
    indices = c(indices, unlist(sapply(strings[[i]], function(x) 
                grep(x, dir(path), invert=invert[i]))))
  }  
  indices = as.numeric(names(which(table(indices) == length(strings))))
  files = dir(path)[indices]
  return(files)
}

loadSimResults = function(S, N, dirPath, B=12, dataType='abu') {
  ## this function is for loading the simulation results associated with the 
  ## parameter space analysis
  ## Arguments:
  ## S: richness values
  ## N: abundance values
  ## dirPath: the directory where the results are stored
  ## B: the number of bisection associated with 
  ## dataType: the kind of data: 'abu' or 'binary'
  results = vector('list', length(S) * length(N))
  icount = 0
  files = dir(dirPath)
  for (s in S) {
    for (n in N) {
      icount = icount + 1
      fileSuffix = paste('_S', s, '_N', n, '_C200_B', B, '_grid_', dataType, sep='')
      fileName = files[grep(fileSuffix, files)]
      if (length(fileName) == 0) {
        next 
      }
      load(file.path(dirPath, fileName))
      results[[icount]] = metrics
      rm(metrics)
    }
  }
  return(results)
}

avgSimResults = function(results, metric, null=TRUE, perm_null=FALSE, bisect=FALSE) {
  ## This function averages the simulation results over all the different runs
  ## for each parameter combination
  ## Arguments
  ## results: the output of the function 'loadSimResults'
  ## metric: the metric analyzed: 'sorensen', 'varWithin', 'jaccard'
  simAvg = vector('list', length(results))
  names(simAvg) = names(results)
  for (i in seq_along(results)) {
    if (is.null(results[[i]]))
      next
    simRaw = reshapeResults(results[[i]], metric, null, perm_null, bisect, 
                            sim_result=TRUE)
    simAvg[[i]] = avgResults(simRaw, null)
  }
  return(simAvg)
}

getSimStats = function(results, S, N) {
  ## computes the slope, intercept, and R2 value for the exponential
  ## and power models across the range of S and N for the simulated results
  ## using Ordinary Least Squares and Weighted regression approaches
  ## the output is a multidimensional array:
  ## the dimensions are [model type, coef type, meth type, grain, S, N]
  ## Arguments:
  ## results: the results file which the models will be fit to
  ## S: the richness values used to parametrize the simulated results
  ## N: the abundance values used to parameterize the simulated results
  grains = unique(results[[1]]$grain)
  stats = array(NA, dim=c(2, 3, 2, length(grains), length(S), length(N)))
  dimnames(stats)[[1]] = c('exp', 'pwr')
  dimnames(stats)[[2]] = c('b0', 'b1', 'r2')
  dimnames(stats)[[3]] = c('ols', 'wtr')
  dimnames(stats)[[4]] = grains
  dimnames(stats)[[5]] = S
  dimnames(stats)[[6]] = N
  i = 0
  for (s in seq_along(S)) {
    for (n in seq_along(N)) {
      i = i + 1
      if (is.null(results[[i]]))
        next
      for (g in seq_along(grains)) {
        j = 1
        for (meth in c('ols', 'wtr')) {  ## wtr is for weighted regression
          if (meth == 'ols') {
            if (min(results[[i]]$avg) > 0) {
              expmod = lm(log10(avg) ~ Dist, subset = grain == grains[g], 
                          data=results[[i]])
              pwrmod = lm(log10(avg) ~ log10(Dist), subset = grain == grains[g],
                          data = results[[i]])
            }  
          }   
          if (meth == 'wtr') {
            if (min(results[[i]]$avg) > 0) {
              expmod = lm(log10(avg) ~ Dist, weights=N, subset = grain == grains[g],
                          data=results[[i]])
              pwrmod = lm(log10(avg) ~ log10(Dist), weights=N, 
                          subset = grain == grains[g], data = results[[i]])
            }  
          }
          if (min(results[[i]]$avg) > 0) {
            stats[1, 1:2, j, g, s, n] = coef(expmod)
            stats[1, 3, j, g, s, n] = summary(expmod)$r.squ
            stats[2, 1:2, j, g, s, n] = coef(pwrmod)
            stats[2, 3, j, g, s, n] = summary(pwrmod)$r.squ
          }  
          else 
            print(paste('skipped', S[s], N[n], sep=' ')) 
          j = j + 1
        }  
      }  
    }
  }
  return(stats)
}

get_ddr_resid = function(obs, exp) {
  ## computes distance decay relationship residuals
  ## this function is not very generic
  ## will only work for the METE spatial project data results
  for (i in seq_along(obs)) {
    index = grep(names(obs[i]), names(exp))
    if (length(index) == 0)
      next
    resid = matrix(NA, ncol=3, nrow=nrow(obs[[i]]))
    icol = 1
    for (metric in c('median', 'average', 'binomial')) {
      if (metric == 'median') {
        metric_column_obs = na.omit(match('Metric.50', names(obs[[1]])))[1]
        metric_column_exp = na.omit(match('Med', names(exp[[1]])))[1]
      }  
      if (metric == 'average') {
        metric_column_obs = na.omit(match('Metric.avg', names(obs[[1]])))[1]
        metric_column_exp = na.omit(match('Avg', names(exp[[1]])))[1]
      } 
      if (metric == 'binomial') {
        metric_column_obs = na.omit(match('Metric.avg', names(obs[[1]])))[1]
        metric_column_exp = na.omit(match('Exp.avg', names(obs[[1]])))[1]
      }
      obs_var = obs[[i]][ , metric_column_obs]
      if (metric != 'binomial') {
        exp_var = exp[[index]][ , metric_column_exp]
      }
      else {
        exp_var = obs[[i]][ , metric_column_exp]
      }  
      resid[ , icol] = obs_var - exp_var
      icol = icol + 1
    }
    resid = as.data.frame(resid)
    names(resid) = c('med.res', 'avg.res', 'exp.res')
    if (i == 1)
      out = data.frame(site=names(obs)[i], obs[[i]], resid)
    else
      out = rbind(out, data.frame(site=names(obs)[i], obs[[i]], resid))
  }  
  return(out)
}

getStats = function(results, metric='median') {
  ## computes the slope, intercept, and R2 value for the exponential
  ## and power models across the range of S and N for the simulated results
  ## using Ordinary Least Squares and Weighted regression approaches
  ## the output is a multidimensional array:
  ## the dimensions are [model type, coef type, meth type, grain, S, N]
  ## Arguments:
  ## results: the results file which the models will be fit to
  ## S: the richness values used to parametrize the simulated results
  ## N: the abundance values used to parameterize the simulated results
  if (metric == 'median')
    metric_column = na.omit(match(c('Med', 'Metric.50'), names(results[[1]])))[1]
  else if (metric == 'average')
    metric_column = na.omit(match(c('Avg', 'Metric.avg'), names(results[[1]])))[1]
  else
    stop('metric agruments should be either "median" or "average"')
  out = vector('list', length(results))
  names(out) = names(results)
  for (i in seq_along(results)) {
    grains = as.character(unique(results[[i]]$Comm))
    resp_var = results[[i]][ , metric_column]
    stats = array(NA, dim=c(2, 3, 2, length(grains)))
    dimnames(stats)[[1]] = c('exp', 'pwr')
    dimnames(stats)[[2]] = c('b0', 'b1', 'r2')
    dimnames(stats)[[3]] = c('ols', 'wtr')
    dimnames(stats)[[4]] = grains
    for (g in seq_along(grains)) {
      j = 1
      for (meth in c('ols', 'wtr')) {  ## wtr is for weighted regression
        if (meth == 'ols') {
          if (min(resp_var) > 0) {
            expmod = lm(log10(resp_var) ~ Dist, subset = as.character(Comm) == grains[g], 
                        data=results[[i]])
   
            pwrmod = lm(log10(resp_var) ~ log10(Dist), subset = as.character(Comm) == grains[g],
                          data = results[[i]])
          }  
        }   
        if (meth == 'wtr') {
          if (min(resp_var) > 0) { 
            expmod = lm(log10(resp_var) ~ Dist, weights=N, 
                        subset = as.character(Comm) == grains[g],
                        data=results[[i]])
            pwrmod = lm(log10(resp_var) ~ log10(Dist), weights=N, 
                        subset = as.character(Comm) == grains[g], data = results[[i]])
          }  
        }
        if (min(resp_var) > 0) {
          stats[1, 1:2, j, g] = coef(expmod)
          stats[1, 3, j, g] = summary(expmod)$r.squ
          stats[2, 1:2, j, g] = coef(pwrmod)
          stats[2, 3, j, g] = summary(pwrmod)$r.squ
        }  
        else {
          print(paste('skipped', names(results)[i], sep=' ')) 
        }  
        j = j + 1        
      }  
    }
    out[[i]] = stats
  }
  return(out)
}

plotEmpir = function(results, metric='median', expected=FALSE, log="", 
                     quants=FALSE, ylim=NULL, title=TRUE, sub="", 
                     alpha=1/3, add=FALSE, col=NULL, lwd=NULL, ...)
{
  ## Purpose: to plot the results, expects that the graphical window has been 
  ## setup appropriately 
  ## Arguments
  ## results: a list from which the results should be plotted
  ## metric: the metric of interest: median or average
  ## log: which axes are to be log transformed
  ## quants: if TRUE will plot quantiles as well
  ## ylim: a list of ylims that is the length of results
  ## title: if TRUE the the plot is given a main title from the names of results
  ## alpha: specifies the transpancy level on the quantiles
  ## col: what colors to use
  ## lwd: the line width of the lines
  ## ...: other arguments supplied to the function lines()
  if (metric == 'median')
    metric_column = na.omit(match(c('Med', 'Metric.50'), names(results[[1]])))[1]
  else if (metric == 'average')
    metric_column = na.omit(match(c('Avg', 'Metric.avg'), names(results[[1]])))[1]
  else
    stop('metric agruments should be either "median" or "average"')
  if (expected)
    expected_column = na.omit(match(c('Exp','Exp.avg'), names(results[[1]])))[1]
  y_is_log = regexpr('y', log)[1] != -1
  if (y_is_log) { 
    mins = unlist(lapply(results, function(x) min(x[ , metric_column])))
    if (any(mins == 0))
      print("Note: Points will be dropped where the metric equals zero becuase the y-axis is log transformed")
  }
  for(i in seq_along(results)){
    grains = results[[i]]$Comm
    unigrains = unique(grains)
    resp_var = results[[i]][ , metric_column]
    if(is.null(col))
      col = palette()[-1]
    if (quants[1]) {
      col_poly = apply(col2rgb(col), 2, 
                  function(x) rgb(x[1], x[2], x[3], 
                                  alpha=alpha * max(x), 
                                  maxColorValue = max(x)))
    }
    if (is.null(ylim)) {
     if (y_is_log)
         yRange = range(resp_var[resp_var > 0])
      else
        yRange = range(resp_var)
    }
    else {
      yRange = ylim[[i]]
    }  
    if (title)
      main = names(results)[i]
    else 
      main = ''
    if (is.null(lwd))
      lwd = 2
    if (!add)
      plot(resp_var ~ Dist, data = results[[i]], subset = resp_var > 0,
           xlim = range(Dist), ylim = yRange, type='n', log=log,
           main=paste(main, '\n', sub, sep=''), ylab=metric)
    for (j in seq_along(unigrains)) {
      dat = subset(results[[i]], Comm == unigrains[j])
      resp_var = dat[ , metric_column]
      if (quants[1]) {
        xvals = dat$Dist
        xvals = c(xvals, rev(xvals))
        yvals = dat[ , (metric_column - 1)]
        yvals = c(yvals, rev(dat[ , (metric_column + 1)]))
        if (y_is_log) { 
          xvals = xvals[yvals > 0]
          yvals = yvals[yvals > 0]
        }
        if (length(xvals) > 0)
          polygon(xvals, yvals, border=NA, col=col_poly[j])
      }
      if (y_is_log)
        lines(resp_var ~ Dist, data = dat, subset=resp_var > 0, col=col[j], lwd=lwd, ...)
      else
        lines(resp_var ~ Dist, data = dat, col=col[j], lwd=lwd, ...)
      if (expected) {
        exp_var = dat[ , expected_column]
        lines(exp_var ~ Dist, data = dat, subset=resp_var > 0, col=col[j], lty=2,
              lwd=lwd, ...)
      }  
    }  
  }
}

str_clean = function(str) {
  str = gsub('\n', '', str)
  str = gsub(' ', '', str)
  return(str)
}

mk_legend = function(...) {
  plot(1:10, 1:10, type='n', axes=F, frame.plot=F, xlab='', ylab='')
  legend(...)
}

comp_results = function(site, results, args=NULL, ...) {
  plot_names = names(results[[1]])
  if (!is.null(args))
    args = sapply(args, str_clean) 
  for (i in seq_along(plot_names)) {
    if (!grepl(site, plot_names[i])) 
     next
    for (j in seq_along(results)) {
      if (j == 1) {
        if (is.null(args[j])) 
          plotEmpir(results[[j]][i], ...)
        else 
          eval(parse(text=paste('plotEmpir(results[[j]][i],', args[j], ',...)',
                                sep='')))
      }
      else {
        index = grep(plot_names[i], names(results[[j]]))
        if (length(index) > 0) {
          if (is.null(args[j])) 
            plotEmpir(results[[j]][index], add=TRUE, ...)
          else 
            eval(parse(text=paste('plotEmpir(results[[j]][index], add=TRUE,', args[j],
                                  ',...)', sep='')))
        }
      }  
    }  
  }
}

addCI = function(x, y_lo, y_hi, col, data=NULL) {
  ## if data is not null then all of the arguments
  ## including data itself must be text strings
  if (!is.null(data)) {
    x = eval(parse(text=paste(data, '$', x, sep='')))
    y_lo = eval(parse(text=paste(data, '$', y_lo, sep='')))
    y_hi = eval(parse(text=paste(data, '$', y_hi, sep='')))
  }  
  xvals = c(x, rev(x))
  yvals = c(y_lo, rev(y_hi))
  polygon(xvals, yvals, border=NA, col=col)
}

addAxes = function(cex.axis=1.75, padj=.5, lwd=4, ...) {
  addAxis(side=1, cex.axis, padj, lwd, ...)
  addAxis(side=2, cex.axis, padj, lwd, ...)
}

addAxis = function(side=1, cex.axis=1.75, padj=.5, lwd=4, ...) {
  axis(side, cex.axis=cex.axis, padj=padj, lwd=lwd, ...)
}

addxlab = function(...) {
  mtext(side=1, cex=2, ...)
}

addylab = function(...) {
  mtext(side=2, cex=2, padj=-2, ...)
}


get_beta_sad_mle = function(S, N) {
  ## edited version of the function get_lambda_sad_mle by Xiao Xiao
  ## edits are only stylistic to improve readibility
  if (N <= 0) {
    print("Error: N must be greater than 0.")
    return(NA)
  }
  else if (S <= 0) {
    print("Error: S must be greater than 0.")
    return(NA)
  }
  else if (S / N >= 1) {
    print("Error: N must be greater than S.")
    return(NA)
  }
  else {
    ## Solve for lambda 1 in Harte et al. 2009 based on (3)
    m = seq(1 : N)  
    b = 10^-99
#    y = function(x) N / x - S / sum(x^m / m) * (1 - x^N) / (1 - x)
#    y = function(x) sum(x ^ m / N * S) - sum(x^m / m)
    y = function(x) sum(x ^ m) / sum(x^m / m) - (N / S)
    
    bound=10^-15  ## Define lower bound for x
    ## If same sign at both ends
    if (y(bound) * y(1 - bound) >= 0) {  
      ## Parameter for log-series
      p = uniroot(y, lower=bound, upper=1.1 - bound,
                  tol=.Machine$double.eps^0.5)$root
    }
    else{
      p = uniroot(y, lower=bound, upper=1 - bound,
                  tol=.Machine$double.eps^0.5)$root   
    }
    beta_sad = -log(p)
    return(list(p=p, beta=beta_sad))
  }
}

avg_site_results = function(dat, site_names) {
  indices = match(site_names, names(dat))
  name_sum = paste(paste('dat$', site_names, sep=''), collapse='+')
  dat[[indices[1]]] = eval(parse(text=paste('(',name_sum, ') /', length(site_names), sep='')))
  dat = dat[-indices[-1]]
  return(dat)
}

get_sar_resids = function(obs_dat, pred_dat, obs_field, pred_field) {
  ## arguments:
  ## obs_dat: the data frame that contains the observed sar
  ## perd_dat: the data frame that contains the predicted sar
  ## obs_field: the name of the field that has the observed richness
  ## pred_field: the name of the field that has the predicted richness
  site = NULL
  area = NULL
  res = NULL
  for (i in seq_along(obs_dat)) {
    site_nm = names(obs_dat)[i]
    site = c(site, rep(site_nm, nrow(obs_dat[[i]])))
    area = c(area, obs_dat[[i]]$area)
    obs_sr = eval(parse(text = paste('obs_dat[[i]]$', obs_field, sep='')))
    if (site_nm %in% names(pred_dat)) {
      index = match(site_nm, names(pred_dat))
      pred_sr = eval(parse(text = paste('pred_dat[[index]]$', pred_field, sep='')))
      len_diff = length(obs_sr) - length(pred_sr)
      pred_sr = c(rep(NA, len_diff), pred_sr)
      res = c(res, obs_sr - pred_sr)
    }
    else
      res = c(res, rep(NA, length(obs_sr)))
  }
  return(data.frame(site, area, res))
}  

get_R2 = function(obs, pred, na.rm=FALSE) {
  if (na.rm) {
    true = !(is.na(obs) | is.na(pred))
    obs = obs[true]
    pred = pred[true]
  }
  SSerr = sum((obs - pred)^2)
  SStot = sum((obs - mean(obs))^2)
  R2 = 1 - SSerr / SStot
  return(R2)
}

get_SSAD = function(comms) {
  grains = unique(comms[ , 1])
  for (g in seq_along(grains)) {
    tmp_comm = comms[comms[ , 1] == grains[g], -(1:3)]
    uni_abu = sort(unique(as.vector(tmp_comm)))
    ssad = apply(tmp_comm, 2, function(x) table(c(x, uni_abu)) - 1)
    ssad = cbind(grains[g], uni_abu, ssad)
    if (g == 1)
      out = ssad
    else
      out = rbind(out, ssad)
  }
  colnames(out)[1:2] = c('grain', 'abu')
  return(out)
}


match_index = function(x, table) {
  ## An extenstion of the function match() that returns the 
  ## index of the match rather than the matching element
  ## Note: not currently being used
  len_x = length(x)
  indices = 1:len_x
  return(indices[match(x, table, nomatch=0) == 1])
}

download_data = function(urls, delim, output_path) {
  ## Function that downloads flat data files from the web
  ## Arguments:
  ## urls: the web address of the flat data file
  ## delim: the column delimintors for each dataset, either
  ##   'comma' or 'tab'
  ## output_path: the directory that the data will be written to
  require(RCurl)
  for (i in seq_along(urls)) {
    temp = getURL(urls[i])
    if (delim[i] == 'comma')
      dat = read.csv(textConnection(temp))
    if (delim[i] == 'tab')
      dat = read.delim(textConnection(temp))
    namesplit = unlist(strsplit(urls[i], split='/'))
    filename = namesplit[length(namesplit)]
    if (delim[i] == 'comma')
      write.csv(dat, file = file.path(output_path, filename), row.names=FALSE)
    if (delim[i] == 'tab')
      write.table(dat, file = file.path(output_path, filename), sep='\t',
                  row.names=FALSE)
    print(paste(filename, ' was downloaded and written to ',
                file.path(output_path, filename), sep=''))
  }
}