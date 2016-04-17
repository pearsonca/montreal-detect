#!/usr/bin/env Rscript

rm(list=ls())

require(data.table)
require(igraph)
require(parallel)

source("../montreal-digest/buildStore.R")

graphPartition <- function(res, mp, referenceCommunities, ulim=60, verbose=F) {
  setkey(referenceCommunities, user_id)
  allnewusers <- mp[user_id %in% setdiff(mp[res[,unique(c(user.a, user.b))], user_id], referenceCommunities$user_id), new_user_id]
  
  gg <- graph(t(res[,list(user.a, user.b)]), directed=F)
  E(gg)$weight <- res$score
  
  comps <- components(gg)
  
  leftovers <- which(comps$csize > ulim)
  completeCommunities <- (1:comps$no)[-leftovers] # components to treat as their own communities
  
  newComms <- comps$membership[allnewusers]
  
  inSmalls <- newComms %in% completeCommunities
  inBigs <- !inSmalls
  
  base <- if (sum(inSmalls)) {
    consensusComms <- sapply(allnewusers[inSmalls], function(newMember){
      comm <- comps$membership[newMember]
      others <- setkey(mp[setdiff(which(comps$membership == comm), newMember), list(user_id)], user_id)
      othercomms <- referenceCommunities[others]
      othercomms[is.na(community), community := -1]
      tmp <- othercomms[,.N,by=community]
      tmp[N==max(N),max(community)]
    })
    data.table(new_user_id=allnewusers[inSmalls], community=consensusComms)
  } else data.table(new_user_id=integer(0), community=integer(0))
  
  if (sum(inBigs)) {
    bigConsensusComms <- sapply(allnewusers[inBigs], function(newMember) {
      others <- setkey(mp[
        setdiff(targettedGraphPartitionOne(newMember, gg, comps), newMember),
        list(user_id)
        ], user_id)
      othercomms <- referenceCommunities[others]
      othercomms[is.na(community), community := -1]
      othercomms[,.N,by=community][N==max(N), max(community)]
    })
    # for each new users
    base <- rbind(base, data.table(new_user_id=allnewusers[inBigs], community=bigConsensusComms))
  }
  originalUserIDs(base, mp)
}

emptygraph <- data.table(user_id=integer(), community=integer())

perturbedPersistenceComms <- function(accPert, bgacc, bgpc, verbose) {
  if (dim(accPert[score > 1])[1]) {
    base.dt <- rbind(accPert[score > 1], bgacc[score > 1])
    ret <- with(relabeller(base.dt), graphPartition(res, mp, bgpc))
    ret
  } else emptygraph
}

for (i in 54:101) {
  args <- c(
    sprintf("input/background-clusters/spin-glass/agg-15-30/%03d.rds",i),
    sprintf("input/background-clusters/spin-glass/pc-15-30/%03d.rds",i),
    sprintf("output/matched/mid/lo/late/10/001-covert-0/%03d-acc.rds",i)
  )
  backgroundAgg <- readRDS(args[1])
  backgroundPC <- readRDS(args[2])
  foregroundAgg <- readRDS(args[3])
  newfile <- sub("-acc.rds",".rds", args[3], fixed = T)
  saveRDS(
    perturbedPersistenceComms(foregroundAgg, backgroundAgg, backgroundPC, verbose = T),
    newfile
  )
  cat("finished",i,"\n")
}

# parse foreground according to interval