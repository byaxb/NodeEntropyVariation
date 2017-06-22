


###########################################
#
#  iGetSccSizes.r
#
#  byaxb (axb@bupt.edu.cn)
#
#  2017-6-14
#
###########################################




#iGetSccSizes
#get the strongly connected component size
#if nodes share the same importance scores
#(as in most of the cases)
#we will rank them randomly
#and get ntimes SCC size sequences

#S3
iGetSccSizes <- function(vImpSeq, ig, decreasing = TRUE, ntimes = 50,  ...) UseMethod("iGetSccSizes")
#for numeric vector
iGetSccSizes.numeric <- function(vImpSeq, ig, decreasing = TRUE, ntimes = 50, ...) {
  
  #to record the start time
  sp <- Sys.time()
  cat("\n[Start at:", as.character(sp))
  
  #function to get order list randomly
  getRandomOrderList <- function(x, decreasing = TRUE, ntimes = 50) {
    sorted_unique <- sort(unique(x), decreasing = decreasing)
    index_list <- list()
    for(i in seq_along(sorted_unique)) {
      index_list[[i]] <- which(x == sorted_unique[i])
    }
    maxLen <- max(unlist(lapply(index_list, length)))
    #if the number of permutaitons is bigger than ntimes
    #we will do the randomization ntimes
    #else we will do #permutaiton times
    ntimes <- min(ntimes, factorial(maxLen))
    orderList <- list()
    for(k in 1:ntimes) {
      for(i in seq_along(sorted_unique)) {
        index_list[[i]] <- index_list[[i]][sample(length(index_list[[i]]))]
      }
      orderList[[k]] <- unlist(index_list)
    }
    return(orderList)
  }
  
  #parallel
  library(foreach)
  library(doParallel)
  cl <- makeCluster(detectCores())
  registerDoParallel(cl, cores = detectCores())
  #to distribute ig to all the cores
  #since ig is one of the args inside this function
  #envir should be set as environment()
  clusterExport(cl, 
                varlist = c("ig"),
                envir=environment()) 
  #get the id
  vidxList <- getRandomOrderList(vImpSeq, decreasing = decreasing, ntimes = ntimes)
  scc_mxs_list <- list()
  for(k in 1:length(vidxList)) {
    vidx <- vidxList[[k]]
    #get the scc_mxs
    scc_mxs <- foreach(
      i = seq_along(vidx),
      .combine = "c",
      .packages = c("foreach", "doParallel", "igraph", "entropy")) %dopar% {
        tmpGraph <- delete.vertices(ig, V(ig)[vidx[1:i]])
        SCC <- clusters(tmpGraph, mode = "strong")
        max(SCC$csize)
      }
    scc_mxs_list[[k]] <- scc_mxs
  }
  
  #stop clusters
  stopCluster(cl)
  ep <- Sys.time()
  cat("\tFinised at:", as.character(ep), "]\n")
  cat("[Time Ellapsed:\t",
      difftime(ep, sp, units = "secs"),
      " seconds]\n")
  return(scc_mxs_list)
}
#for data.frame
iGetSccSizes.data.frame <- function(vImpDataFrame, ig, decreasing = TRUE, ntimes = 50,  ...) {
  igName <- deparse(substitute(ig))
  ire <- lapply(colnames(vImpDataFrame), function(x) {
    cat("SccSizes computation: ", x, " for ", igName)
    curColRe <- NA
    tryCatch( curColRe <- iGetSccSizes.numeric(vImpDataFrame[, x], ig,  decreasing, ntimes),
              error = function(err) {
                cat("\n\n Error while conduct computaion on ", x, "\n\n")
              })
    return(curColRe)
  })
  names(ire) <- colnames(vImpDataFrame)
  return(ire)
}
#default
iGetSccSizes.default <- function(vImpDataFrame, ig, decreasing = TRUE, ...) {
  cat("Currently only support numeric vectors or data.frame!\n")
}

