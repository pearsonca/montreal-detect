#!/usr/bin/env Rscript

rm(list=ls())

require(data.table)
require(igraph)
require(parallel)

source("../montreal-digest/buildStore.R")

emptycomp <- data.table(new_user_id=integer(0), community=integer(0))
emptygraph <- data.table(user_id=integer(), community=integer())

decompose <- function(allnewusers, comps, gg, completeCommunities, referenceCommunities, mp) {
  newComms <- comps$membership[allnewusers]
  inSmalls <- newComms %in% completeCommunities
  largeCommunities(
    !inSmalls, allnewusers, comps, gg, referenceCommunities, mp,
    smallComponents(inSmalls, allnewusers, comps, referenceCommunities, mp)
  )
}

smallComponents <- function(inSmalls, allnewusers, comps, referenceCommunities, mp) if (!sum(inSmalls)) {
  emptycomp
} else {
  consensusComms <- sapply(allnewusers[inSmalls], function(newMember){
    comm <- comps$membership[newMember]
    others <- setkey(mp[setdiff(which(comps$membership == comm), newMember), list(user_id)], user_id)
    othercomms <- referenceCommunities[others]
    othercomms[is.na(community), community := -1]
    tmp <- othercomms[,.N,by=community]
    tmp[N==max(N),max(community)]
  })
  data.table(new_user_id=allnewusers[inSmalls], community=consensusComms)
}

largeCommunities <- function(inBigs, allnewusers, comps, gg, referenceCommunities, mp, base) if (!sum(inBigs)) base else {
  candidates <- allnewusers[inBigs]
  res <- emptycomp
  while (length(candidates)) {
    tar <- candidates[1]
    found <- targettedGraphPartitionOne(tar, gg, comps)
    others <- setkey(mp[
      setdiff(found, tar),
      list(user_id)
      ], user_id)
    othercomms <- referenceCommunities[others]
    othercomms[is.na(community), community := -1]
    fndcomm <- othercomms[,.N,by=community][N==max(N), max(community)]
    known <- intersect(candidates, found)
    res <- rbind(res, data.table(new_user_id = known, community = fndcomm))
    candidates <- setdiff(candidates, found)
  }
  rbind(base, res)
}

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
  
  base <- decompose(allnewusers, comps, gg, completeCommunities, referenceCommunities, mp)
  originalUserIDs(base, mp)
}

perturbedPersistenceComms <- function(accPert, bgacc, bgpc, verbose) {
  if (dim(accPert[score > 1])[1]) {
    base.dt <- rbind(accPert[score > 1], bgacc[score > 1])
    ret <- with(relabeller(base.dt), graphPartition(res, mp, bgpc))
    ret
  } else emptygraph
}

parse_args <- function(argv = commandArgs(trailingOnly = T)) {
  parser <- optparse::OptionParser(
    usage = "usage: %prog path/to/rawevents.rds path/locationslifetimes.rds path/to/precomputed path/to/simout path/to/target",
    description = "compute snapshot communities + persistence scores for a simulation",
    option_list = list(
      optparse::make_option(
        c("--verbose","-v"),  action="store_true", default = FALSE,
        help="verbose?"
      )
    )
  )
  req_pos <- list(
    agg.dt=readRDS,
    pc.dt=readRDS,
    pertpath=identity
  )
  parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
  parsed$options$help <- NULL
  result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
  if(result$verbose) print(result)
  result
}

resolve <- function(agg.dt, pc.dt, pertpath, verbose) {
  foregroundAgg <- readRDS(pertpath)
  newfile <- sub("-acc.rds",".rds", pertpath, fixed = T)
  saveRDS(
    perturbedPersistenceComms(foregroundAgg, agg.dt, pc.dt, verbose = T),
    newfile
  )
  if (verbose) cat(sprintf("finishing %s\n", newfile))
}

#for (i in 86:113) {
do.call(resolve, parse_args(
#  c(sprintf("input/background-clusters/spin-glass/agg-15-30/%03d.rds",i),sprintf("input/background-clusters/spin-glass/pc-15-30/%03d.rds",i),sprintf("output/matched/mid/lo/late/10/001-covert-0/%03d-acc.rds",i), "-v")
))
#}
#   args <- c(
#     sprintf("input/background-clusters/spin-glass/agg-15-30/%03d.rds",i),
#     sprintf("input/background-clusters/spin-glass/pc-15-30/%03d.rds",i),
#     sprintf("output/matched/mid/lo/late/10/001-covert-0/%03d-acc.rds",i)
#   )
#   backgroundAgg <- readRDS(args[1])
#   backgroundPC <- readRDS(args[2])
#   foregroundAgg <- readRDS(args[3])
#   newfile <- sub("-acc.rds",".rds", args[3], fixed = T)
#   saveRDS(
#     perturbedPersistenceComms(foregroundAgg, backgroundAgg, backgroundPC, verbose = T),
#     newfile
#   )
#   cat("finished",i,"\n")
# }
# 
# # parse foreground according to interval