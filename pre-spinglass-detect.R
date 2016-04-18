#!/usr/bin/env Rscript

rm(list=ls())

# (1) raw background
# (2) location lifetimes to censor covert
# (3) reference persistence communities
# (4) reference simulation source

require(data.table)
require(igraph)
require(parallel)

readForeground <- function(simpathroot) {
  cc.dt <- setkey(
    fread(paste0(sub("/$","",simpathroot),"-cc.csv"),
          col.names = c("user.a","user.b","location_id","start","end","type")
    ), location_id)
  cu.dt <- fread(paste0(sub("/$","",simpathroot),"-cu.csv"), col.names = c("user.a","user.b","location_id","start","end","type"))
  list(cc.dt=cc.dt, cu.dt=cu.dt)
}

filelister <- function(srcpath) list.files(srcpath, full.names = T)

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
    raw.dt=readRDS,
    location_lifetimes.dt=readRDS,
    srcpath=identity,
    foreground=readForeground,
    outpath=identity
  )
  parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
  parsed$options$help <- NULL
  result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
  if(result$verbose) print(result)
  result
}

source("../montreal-digest/buildStore.R")

slice <- function(dt, low, high) relabeller(
  dt[start < high*24*3600 & low*24*3600 < end,
     list(score=.N, covert.a=user.a < 0, covert.b=user.b < 0),
     by=list(user.a, user.b)]
)

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

with(parse_args(
#  c("input/raw-pairs.rds", "input/raw-location-lifetimes.rds", "input/background-clusters/spin-glass/base-15-30", "output/matched/mid/lo/late/10/001-covert-0", "output/matched/mid/lo/late/10/001-covert-0-base.rds")
),{
  maxuid <- raw.dt[,max(user.a, user.b)]
  bgcommunities <- filelister(srcpath)
  interval <- as.integer(sub(".+-(\\d+)-\\d+$","\\1", srcpath))
  window <- as.integer(sub(".+-\\d+-(\\d+)$","\\1", srcpath))
  
  trimforegroundcc.dt <- location_lifetimes.dt[foreground$cc.dt][
    end > arrive & start < depart,
    list(user.a, user.b, location_id, start, end, type)
  ]
  
  foreground.dt <- setcolorder(rbind(foreground$cu.dt, trimforegroundcc.dt)[,
      list(user.a, user.b, start, end)
    ], c("user.a", "user.b", "start", "end")
  )
  
  temp <- rbind(foreground.dt, raw.dt)
  temp[
    user.b < user.a,
    `:=`(user.b = user.a, user.a = user.b)
  ]
  
  ref.dt <- setkey(temp[user.a != user.b], start, end)
  
  perturbedCommunities <- resolve(ref.dt, interval, window, bgcommunities)
  
  saveRDS(perturbedCommunities, outpath)
})