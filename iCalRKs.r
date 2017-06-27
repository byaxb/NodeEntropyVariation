


###########################################
#
#  iCalRKs.r
#
#  byaxb (axb@bupt.edu.cn)
#
#  2017-6-14
#
###########################################




iCalRKs <- function(sccSizes, ig) {
  library(igraph)
  igName <- deparse(substitute(ig))
  cat("r(k) scores computation:", igName)
  sp <- Sys.time()
  cat("\n[Start at:", as.character(sp))
  
  library(data.table)
  SCC <- clusters(ig, mode = "strong")
  initialSccSize <- max(SCC$csize)
  rk <- NULL
  for (curCen in names(sccSizes)) {
    tmp <- t(as.matrix(rbindlist(sccSizes[curCen])))
    tmp <- 1 - tmp / initialSccSize
    rk <- cbind(rk, apply(tmp, 2, mean))
  }
  rk <- as.data.frame(rk)
  colnames(rk) <- names(sccSizes)
  ep <- Sys.time()
  cat("\tFinised at:", as.character(ep), "]\n")
  cat("[Time Ellapsed:\t",
      difftime(ep, sp, units = "secs"),
      " seconds]\n")
  return(rk)
}


#Examples
#load igraph objects
load("inet.rda")
#For the Hens
#to calculate different centralities
#make sure that abbreviated names are correct
# DEall: Degree (all)
# DEin: Degree (in)
# DEout: Degree (out)
# COall: Coreness (all)
# COin: Coreness (in)
# COout: Coreness (out)
# HIall: H-index(all)
# HIin: H-index (in)
# HIout: H-index (out)
# ECall: Eccentricity (all)
# ECin: Eccentricity (in)
# ECout: Eccentricity (out)
# CLall: Closeness (all)
# CLin: Closeness (in)
# CLout: Closeness (out)
# IN: Information index
# BE: Betweenness centrality
# LO: Load centrality
# ST: Stress centrality
# SU: Subgraph centrality
# EI: Eigenvector centrality
# AL: Alpha centrality
# PR: Page Rank
# AU: HITs (Authority)
# HU: HITs (Hub)
# ExFs: Expected Force (ExF)
# EnVbtw: Entropy Variation (btw)
# EnVall: Entropy Variation (all)
# EnVout: Entropy Variation(out)
# EnVin: Entropy Variation(in)
HenCentralities <- iCalCen(Hens,c("DEall","DEin","DEout",
                                  "COall","COin","COout",
                                  "HIall","HIin","HIout",
                                  "ECall","ECin","ECout",
                                  "CLall","CLin","CLout",
                                  "IN","BE","LO","ST",
                                  "SU","EI","AL","PR",
                                  "ExFs", 
                                  "EnVbtw", "EnVall", "EnVout", "EnVin"))

HensSccSizes <- iGetSccSizes(HenCentralities, Hens)
HensRKs <- iCalRKs(HensSccSizes, Hens)
#r(k=10)
HensRKs[10, ]

