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
    if(is.null(combine[[i]]))
      combine[[i]] = results[[i]]$Comm
    avg = tapply(results[[i]]$Metric,list(results[[i]]$Dist,combine[[i]]),mean)
    dist = as.numeric(row.names(avg))
    row.names(avg) = NULL
    out[[i]] = cbind(dist,avg)
  }
  return(out)
}

shrtnames = c('bci','cocoli','sherman','serp')
empirBin = getResults(shrtnames,'sorensen','binary')
empirAbu = getResults(shrtnames,'sorensen','abu')
empirSorBin = reshapeResults(empirBin,'sorensen')
empirSorAbu = reshapeResults(empirAbu,'sorensen')

combine = vector('list',length=length(empirSorBin))
combine$bci = NULL
combine$cocoli = ifelse(empirSorBin$cocoli$Comm > 6,empirSorBin$cocoli$Comm - 6,
                        empirSorBin$cocoli$Comm)
combine$sherman = ifelse(empirSorBin$sherman$Comm > 6 & empirSorBin$sherman$Comm < 13,
                         empirSorBin$sherman$Comm - 6, empirSorBin$sherman$Comm)
combine$serp = NULL

empirSorBinAvg = avgResults(empirSorBin,combine)

plot(empirSorBinAvg[[2]][,1],empirSorBinAvg[[2]][,2],ylim=range(empirSorBinAvg[[2]][,-1],na.rm=T))
plot(empirSorBinAvg[[2]][,1],empirSorBinAvg[[2]][,6])


vGridExpAvg = apply(vGridExp,1,mean)
vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))


