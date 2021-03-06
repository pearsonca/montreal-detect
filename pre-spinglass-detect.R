#!/usr/bin/env Rscript

rm(list=ls())

# (1) raw background
# (2) location lifetimes to censor covert
# (3) reference persistence communities
# (4) reference simulation source

require(data.table)
require(igraph)
require(parallel)

source("../montreal-digest/buildStore.R")

loadBase <- function(srcpath) list.files(srcpath, pattern='\\d{3}\\.rds$', full.names = T)
  
parse_args <- function(argv = commandArgs(trailingOnly = T)) {
  parser <- optparse::OptionParser(
    usage = "usage: %prog path/to/rawevents.rds path/locationslifetimes.rds path/to/precomputedbase path/to/cc path/to/cu",
    description = "compute snapshot communities + persistence scores for a simulation",
    option_list = list(
      optparse::make_option(
        c("--verbose","-v"),  action="store_true", default = FALSE,
        help="verbose?"
      )
    )
  )
  req_pos <- list(
    raw.dt=readRDS,
    bgcomms=loadBase,
    foreground.dt=readRDS,
    intDays=as.integer,
    winDays=as.integer
  )
  parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
  parsed$options$help <- NULL
  result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
  if(result$verbose) print(result)
  result
}

slice <- function(dt, low, high) relabeller(
  dt[start < high*24*3600 & low*24*3600 < end,
     list(score=.N, covert.a=user.a < 0, covert.b=user.b < 0),
     by=list(user.a, user.b)]
)

emptycomp <- data.table(new_user_id=integer(), community=integer())
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
    found <- unique(c(targettedGraphPartitionOne(tar, gg, comps), tar))
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

resolve <- function(
  base.dt, intDays, winDays, bgcommunities, #outputdir, crs,
  mxint=NA, verbose=F, ...
) {
  st = base.dt[1, floor(start/60/60/24)]
  targets <- 1:length(bgcommunities)
  thing <- mapply(function(inc, bgcomm) {
    with(slice(base.dt, st + inc*intDays-winDays, st + inc*intDays), {
      store <- if (dim(res)[1] == 0) emptygraph else {
        coverts <- res[,unique(c(user.a[covert.a],user.b[covert.b]))]
        if (length(coverts)) {
          if (verbose) cat("something to do ", inc, length(coverts), mp[coverts, user_id], "\n") 
          graphPartition(res, mp, readRDS(bgcomm))
        } else emptygraph
      }
    })[, increment:=inc ]
  }, targets, bgcommunities, SIMPLIFY = F)
  rbindlist(thing)
}

saveRDS(with(parse_args(
# c("input/digest/raw/pairs.rds", "input/digest/background/15/15/base", "input/detection/high/hi/early/10/15/15/010/trim.rds", "15", "15")
),{

  temp <- rbind(foreground.dt[,-"increment",with=F], raw.dt)
  temp[
    user.b < user.a,
    `:=`(user.b = user.a, user.a = user.b)
  ]
  
  ref.dt <- setkey(temp[user.a != user.b], start, end)
  
  resolve(ref.dt, intDays, winDays, bgcomms)
}), pipe("cat","wb"))