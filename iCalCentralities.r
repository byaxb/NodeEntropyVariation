


###########################################
#
#  iCalCentralities.r
#
#  byaxb (axb@bupt.edu.cn)
#
#  2017-6-14
#
###########################################



#To calculate the ExFs in a parallel mode
iCalExF <- function(ig) {
  #The next function, ExF, is from 
  #https://github.com/glennlawyer/ExpectedForce
  #Expected Forece is proposed by Lawyer in:
  #Lawyer, G. Understanding the influence of all nodes in a network.
  #Scientific Reports 2015, 5, 8665
  ExF <- function(qnode,graph){
    .FI <- function(graph,clust){
      return(length(unique(unlist(neighborhood(graph,clust,order=1)))) - 3) }
    ## Here, we know the cluster will have three elements, hence the "- 3"
    ## if the cluster size is unknown, better to use "- length(clust)"
    
    ## Get all neighbors of the querry node at distance one and two:
    neigh <- graph.bfs(graph,qnode,order=FALSE,dist=TRUE,
                       callback=function(graph,data,extra){
                         data["dist"]==3})$dist
    ## vector of nodes at distance one
    d.one.nodes <- which(neigh==1)
    n.d.one <- length(d.one.nodes)
    ## vector of nodes at distance two
    d.two.nodes <- which(neigh==2)
    
    ## pre-allocate the vector of FI values
    guestimated.numFI <- 2*sum(n.d.one*length(d.two.nodes))
    allFI <- numeric(guestimated.numFI+5)
    numFI <- 0;  totalFI <- 0
    
    ## The iteration is over all nodes at distance one from the source,
    ##     within this loop we consider both all remaining d.one.nodes
    ##     and all d.two.nodes reachable from the current d.one.node.
    for(i in 1:n.d.one){
      if(i<n.d.one){   ## all remaining d.one.nodes
        for(j in (i+1):n.d.one){
          ## Increase storage, if necessary: (code for optimization only)
          if(numFI>guestimated.numFI){
            guestimated.numFI <-  round(1.5 * guestimated.numFI)
            foo <- allFI
            allFI <- vector(mode="numeric",length=guestimated.numFI+5)
            allFI[1:numFI] <- foo[1:numFI]
          } ## END increase storage
          ## compute cluster FI
          clustFI <- .FI(graph, c(qnode, d.one.nodes[c(i,j)]))
          ## is there an edge between nodes i and j?
          mult <- 
            if(length(E(graph)[ d.one.nodes[i] %--% d.one.nodes[j] ])) 4 else 2
          ## store cluster FI the appropriate number of times
          allFI[seq(numFI+1,length.out=mult)] <- clustFI
          totalFI <- totalFI + mult*clustFI
          numFI <- numFI+mult
        }} ## end all remaining d.one.nodes
      for(dtn in d.two.nodes){   ## all d.two.nodes 
        if(length(E(graph)[ d.one.nodes[i] %--% dtn ])){
          ## If an edge to the current d.one.node
          ## increase storage, if necessary: (code for optimization only)
          if(numFI>guestimated.numFI){
            guestimated.numFI <-  round(1.5 * guestimated.numFI)
            foo <- allFI
            allFI <- vector(mode="numeric",length=guestimated.numFI+5)
            allFI[1:numFI] <- foo[1:numFI]
          } ## END increase storage
          ## compute cluster FI
          clustFI <- .FI(graph, c(qnode, d.one.nodes[i], dtn))
          numFI <- numFI+1
          allFI[numFI] <- clustFI
          totalFI <- totalFI + clustFI
        }}
    } ## end looping over all nodes at distance one
    ## calculate the entropy, note that this clips allFI at numFI
    norm <- allFI[1:numFI]/totalFI
    -sum(norm*log(norm))
    return(-sum(norm*log(norm)))
  }
  
  library(foreach)
  library(doParallel)
  cl <- makeCluster(detectCores())
  registerDoParallel(cl, cores = detectCores())
  clusterExport(cl,
                varlist = c("ig", "ExF"),
                envir=environment())
  ExFs <- foreach(
    x = 1:igraph::vcount(ig),
    .combine = "c",
    .packages = c("foreach", "doParallel", "igraph")) %dopar% {
      ExF(x, ig)
    }
  stopCluster(cl)
  #return the ExFs
  return(ExFs)
}

iCalEnV <- function(ig, mode = "in") {
  library(foreach)
  library(doParallel)
  cl <- makeCluster(detectCores())
  registerDoParallel(cl, cores = detectCores())
  clusterExport(cl,
                varlist = c("ig", "mode"),
                envir=environment())
  nrEntropies <- foreach(
    x = igraph::V(ig),
    .combine = "c",
    .packages = c("foreach", "doParallel", "igraph", "entropy")) %dopar% {
      tmpG <- igraph::delete_vertices(ig, x)
      if(mode %in% c("all", "in", "out")) {
        probs <- igraph::degree(tmpG, mode = mode) / sum(igraph::degree(tmpG, mode = mode))
      } else if(mode == "btw") {
        probs <- igraph::betweenness(tmpG) / sum(igraph::betweenness(tmpG))
      } else {
        #never be here
        #because only btw/all/in/out are allowed in methodsDF
      }
      return(entropy::entropy.plugin(probs))
    }
  stopCluster(cl)
  
  library(entropy)
  tmpG <- ig
  if(mode != "btw") {
    probs <- igraph::degree(tmpG, mode = mode) / sum(igraph::degree(tmpG, mode = mode))
  } else {
    probs <- igraph::betweenness(tmpG) / sum(igraph::betweenness(tmpG))
  }
  initEntropy <- entropy::entropy.plugin(probs)
  EnV <- initEntropy - nrEntropies
  return(EnV)
}


iCalCen <- function(ig, centralities = c("DEall","DEin","DEout",
                                         "COall","COin","COout",
                                         "HIall","HIin","HIout",
                                         "ECall","ECin","ECout",
                                         "CLall","CLin","CLout",
                                         "IN","BE","LO","ST",
                                         "SU","EI","AL","PR",
                                         "AU","HU",
                                         "ExFs", 
                                         "EnVbtw", "EnVall", "EnVout", "EnVin")) {

  igName <- deparse(substitute(ig))
  cat("\n\n###########################################")
  cat("\n\nNow calculating the centralities of ", igName)
  sp <- Sys.time()
  cat("\n[Start at:", as.character(sp), "]\n")
  
  #to redefine functions to get the final results of 
  #eigen_centrality, page_rank, authority_score,
  #hub_score and hIndex
  #which will be more convinent
  #to be called by do.call
  iCalEigen <- function(ig) {
    igraph::eigen_centrality(ig, directed = TRUE)$vector
  }
  iCalPR <- function(ig) {
    igraph::page_rank(ig)$vector
  }
  iCalAuthority <- function(ig) {
    igraph::authority_score(ig)$vector
  }
  iCalHu <- function(ig) {
    igraph::hub_score(ig)$vector
  }
  hIndex <- function(ig, mode = c("all")) {
    library(foreach)
    library(doParallel)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl, cores = detectCores())
    clusterExport(cl,
                  varlist = c("ig", "mode"),
                  envir=environment())
    hi <- foreach(
      x = igraph::V(ig),
      .combine = "c",
      .packages = c("foreach", "doParallel", "igraph", "entropy")
    ) %dopar% {
      agop::index_h(igraph::degree(ig, igraph::neighbors(ig, x), mode = mode))
    }
    stopCluster(cl)
    return(hi)
  }
  
  #some of the centralities will be computed 
  #with the add-on package sna
  #it is necessary to convert an igraph object
  #to a network object
  #the convertion will be done with intergraph
  library(sna)
  library(intergraph)
  inet <- asNetwork(ig)
  ep <- Sys.time()
  cat("\t\nPreparation done@", as.character(ep), "[", difftime(ep, sp, units = "secs"), "secs]")

  
  #a data.frame methodsDF will be defined
  #to mapp the centralities to their methods
  #along with some of its args
  #if more centralities are to be calculated,
  #just exetend the methodsDF
  #adding a new line will be enough
  #methodsDF comprises of three columns:
  #method, igraph or network, args
  #the centrality is defined as the row.names
  #thus the mapping is done
  centralityNames <- c("DEall","DEin","DEout",
                       "COall","COin","COout",
                       "HIall","HIin","HIout",
                       "ECall","ECin","ECout",
                       "CLall","CLin","CLout",
                       "IN","BE","LO","ST",
                       "SU","EI","AL","PR",
                       "AU","HU",
                       "ExFs", 
                       "EnVbtw", "EnVall", "EnVout", "EnVin")
  methodsDF <- data.frame(method = c("igraph::degree","igraph::degree","igraph::degree",
                                     "igraph::coreness","igraph::coreness","igraph::coreness",
                                     "hIndex","hIndex","hIndex",
                                     "igraph::eccentricity","igraph::eccentricity","igraph::eccentricity",
                                     "igraph::closeness","igraph::closeness","igraph::closeness",
                                     "sna::infocent","igraph::betweenness","sna::loadcent","sna::stresscent",
                                     "igraph::subgraph_centrality","iCalEigen","igraph::alpha_centrality","iCalPR",
                                     "iCalAuthority","iCalHu",
                                     "iCalExF", 
                                     "iCalEnV", "iCalEnV", "iCalEnV", "iCalEnV"),
                          graph = c(rep("ig", 15),
                                    "inet", "ig", "inet", "inet",
                                    rep("ig", 11)),
                          arg = c("all","in","out",
                                   "all","in","out",
                                   "all","in","out",
                                   "all","in","out",
                                   "all","in","out",
                                   NA,NA,NA,NA,
                                   NA,NA,NA,NA,
                                   NA,NA,
                                   NA, 
                                   "btw", "all", "out", "in"))
  row.names(methodsDF) <- centralityNames
  
  #define a function getfun 
  #to tell do.call to work in the right way
  #while passing "pkg::function" to it 
  getfun<-function(x) {
    if(length(grep("::", x))>0) {
      parts<-strsplit(x, "::")[[1]]
      getExportedValue(parts[1], parts[2])
    } else {
      x
    }
  }
  
  #all the centrality scores will be stored in allCen
  allCen <- NULL
  #which centralities will be finally selected
  #centralities stands for the ones passed to iCalCen
  #centralityNames stans for the supported ones
  selectedCen <- intersect(centralities, centralityNames)
  for(curCen in selectedCen) {
    tmpCen <- rep(NA, igraph::vcount(ig))
    tryCatch({
      if(is.na(methodsDF[curCen, "arg"])) {
        tmpCen <- do.call(getfun(methodsDF[curCen, "method"]), list(get(methodsDF[curCen, "graph"])))
      } else {
        tmpCen <- do.call(getfun(methodsDF[curCen, "method"]), list(get(methodsDF[curCen, "graph"]), 
                                                                    mode = methodsDF[curCen, "arg"]))
      }
    },            error = function(err){
      cat("\n\t!!!!!!Error while calculating ", curCen, ":\n", err$message)
    })
    allCen <- cbind(allCen, tmpCen)
    ep <- Sys.time()
    cat("\t\n", curCen, "  finished @", as.character(ep), "[", difftime(ep, sp, units = "secs"), "secs]")
  }
  allCen <- as.data.frame(allCen)
  colnames(allCen) <- selectedCen
  cat("\n\n###########################################\n")
  return(allCen)
}
