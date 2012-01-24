## Purpose: to graphically summarize the empirical results

setwd('/home/danmcglinn/maxent/spat')

getResults = function(names,metric,dataType)
{
  results = vector('list',length=length(names))
  names(results) = names
  for(i in seq_along(results)){
    load(paste('./',metric,'/',metric,'_',names[i],'_',dataType,'.Rdata',
               sep=''))
    results[[i]] = eval(parse(text=metric))
  }
  return(results)
}

reshapeResults = function(results,metric){
  out= vector('list',length=length(results))
  names(out) = names(results)
  for(i in seq_along(results)){ ## the datset
    vExp = vector('list',length=length(results[[i]]))
    Dist = vExp
    n = vExp
    for(j in seq_along(results[[i]])){ ## the grain/community
      Dist[[j]] = results[[i]][[j]][[1]]$vario$Dist
      n[[j]] = results[[i]][[j]][[1]]$vario$n
      if(any(metric %in% c('sorensen','jaccard')))
        vExp[[j]] = 1 - results[[i]][[j]][[1]]$vario$exp.var
      else
        vExp[[j]] = results[[i]][[j]][[1]]$vario$exp.var
    }
    if(is.null(names(results[[i]])))
      commNames = 1:length(results[[i]])
    out[[i]] = data.frame(Dist = unlist(Dist),Metric = unlist(vExp),N = unlist(n))
    out[[i]] = data.frame(out[[i]],
               Comm = unlist(mapply(rep,commNames,each=sapply(vExp,length))))
  }
  return(out)
}  

avgResults = function(results,combine = NULL){
  ## combine is a matrix, each row specifies an index
  ## each column specifies a different one to average over
  out = vector('list',length(results))
  names(out) = names(results)
  for(i in seq_along(results)){
    if(is.na(combine[[i]][1]))
      combine[[i]] = results[[i]]$Comm
    unicombine = unique(combine[[i]])
    for(j in seq_along(unicombine)){
      true = combine[[i]] == unicombine[j]
      Dist = tapply(results[[i]]$Dist[true],round(results[[i]]$Dist[true],3),mean)
      Metric = tapply(results[[i]]$Metric[true],round(results[[i]]$Dist[true],3),mean)
      N = tapply(results[[i]]$N[true],round(results[[i]]$Dist[true],3),sum)
      Comm = rep(unicombine[j],length(Metric))
      if(j == 1)
        out[[i]] = data.frame(Dist,Metric,N,Comm)
      else
        out[[i]] = rbind(out[[i]],data.frame(Dist,Metric,N,Comm))
    }  
  }
  return(out)
}

shrtnames = c('bci','cocoli','sherman','serp')
empirBin = getResults(shrtnames,'sorensen','binary')
empirAbu = getResults(shrtnames,'sorensen','abu')
empirSorBin = reshapeResults(empirBin,'sorensen')
empirSorAbu = reshapeResults(empirAbu,'sorensen')

combine = vector('list',length=length(empirSorBin))
combine[[1]] = NA
combine[[2]] = ifelse(empirSorBin$cocoli$Comm > 6,empirSorBin$cocoli$Comm - 6,
                      empirSorBin$cocoli$Comm)
combine[[3]] = ifelse(empirSorBin$sherman$Comm > 6 & empirSorBin$sherman$Comm < 13,
                      empirSorBin$sherman$Comm - 6, empirSorBin$sherman$Comm)
combine[[4]] = NA

empirSorBinAvg = avgResults(empirSorBin,combine)
empirSorAbuAvg = avgResults(empirSorAbu,combine)

#pdf('spat_empir_bin_DD_curves.pdf')
par(mfrow=c(2,2))
for(i in seq_along(empirSorBinAvg)){
  unigrains = unique(empirSorBinAvg[[i]]$Comm)
  plot(Metric ~ Dist,data = empirSorBinAvg[[i]],subset= Comm == i,
     xlim = range(Dist), ylim = range(Metric),type='n',main=names(empirSorBinAvg)[i])
  for(j in seq_along(unigrains)){
    lines(Metric ~ Dist,data = empirSorBinAvg[[i]],subset= Comm == unigrains[j],
          col=j,lwd=2)
  }  
}


par(mfrow=c(2,2))
for(i in seq_along(empirSorBinAvg)){
  unigrains = unique(empirSorBinAvg[[i]]$Comm)
  plot(Metric ~ Dist,data = empirSorBinAvg[[i]],subset= Comm == i,log='xy',
     xlim = range(Dist), ylim = range(Metric),type='n',main=names(empirSorBinAvg)[i])
  for(j in seq_along(unigrains)){
    lines(Metric ~ Dist,data = empirSorBinAvg[[i]],subset= Comm == unigrains[j],
          col=j,lwd=2)
  }  
}

#dev.off()

#pdf('spat_empir_abu_DD_curves.pdf')
par(mfrow=c(2,2))
for(i in seq_along(empirSorAbuAvg)){
  unigrains = unique(empirSorAbuAvg[[i]]$Comm)
  plot(Metric ~ Dist,data = empirSorAbuAvg[[i]],subset= Comm == i,
     xlim = range(Dist), ylim = range(Metric),type='n',main=names(empirSorAbuAvg)[i])
  for(j in seq_along(unigrains)){
    lines(Metric ~ Dist,data = empirSorAbuAvg[[i]],subset= Comm == unigrains[j],
          col=j,lwd=2)
  }  
}


par(mfrow=c(2,2))
for(i in seq_along(empirSorAbuAvg)){
  unigrains = unique(empirSorAbuAvg[[i]]$Comm)
  plot(Metric ~ Dist,data = empirSorAbuAvg[[i]],subset= Comm == i,log='xy',
     xlim = range(Dist), ylim = range(Metric),type='n',main=names(empirSorAbuAvg)[i])
  for(j in seq_along(unigrains)){
    lines(Metric ~ Dist,data = empirSorAbuAvg[[i]],subset= Comm == unigrains[j],
          col=j,lwd=2)
  }  
}

#dev.off()


vGridExpAvg = apply(vGridExp,1,mean)
vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))


