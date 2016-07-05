#!/usr/bin/env Rscript

rm(list=ls())

require(data.table)
require(igraph)
require(parallel)

source("../montreal-digest/buildStore.R")

emptycomp <- data.table(new_user_id=integer(0), community=integer(0))
emptygraph <- data.table(user_id=integer(), community=integer())

decompose <- function(allnewusers, comps, gg, completeCommunities, referenceCommunities, mp, verbose) {
  newComms <- comps$membership[allnewusers]
  inSmalls <- newComms %in% completeCommunities
  largeCommunities(
    !inSmalls, allnewusers, comps, gg, referenceCommunities, mp,
    smallComponents(inSmalls, allnewusers, comps, referenceCommunities, mp, verbose),
    verbose
  )
}

smallComponents <- function(inSmalls, allnewusers, comps, referenceCommunities, mp, verbose = F) if (!sum(inSmalls)) {
  if (verbose) cat("no small communities\n", file=stderr())
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

largeCommunities <- function(inBigs, allnewusers, comps, gg, referenceCommunities, mp, base, verbose = F) if (!sum(inBigs)) {
    if (verbose) cat("no big communities\n", file=stderr())
    base
  } else {
  candidates <- allnewusers[inBigs]
  res <- emptycomp
  while (length(candidates)) {
    tar <- candidates[1]
    if (verbose) cat("partioning for", tar, "\n", file=stderr())
    found <- unique(c(targettedGraphPartitionOne(tar, gg, comps, verbose), tar))
    if (verbose) {
      cat("found:\n", found, file=stderr())
    }
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
    if (verbose) {
      cat("remaining...\n", candidates, file=stderr())
    }
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

  base <- decompose(allnewusers, comps, gg, completeCommunities, referenceCommunities, mp, verbose)
  originalUserIDs(base, mp)
}

perturbedPersistenceComms <- function(accPert, bgacc, bgpc, score_mode, verbose) {
  if (score_mode != "drop-only" & dim(accPert[score > 1])[1]) {
    base.dt <- rbind(accPert[score > 1, score, by=list(user.a, user.b)], bgacc[score > 1, score, by=list(user.a, user.b)])
    with(relabeller(base.dt), graphPartition(res, mp, bgpc, verbose))
  } else if (dim(accPert)[1]) {
    base.dt <- rbind(accPert[, score, by=list(user.a, user.b)], bgacc[, score, by=list(user.a, user.b)])
    with(relabeller(base.dt), graphPartition(res, mp, bgpc, verbose))
  } else emptygraph
}

listPC <- function(srcpath) list.files(srcpath, pattern = "\\d{3}\\.rds$", full.names = T)
listACC <- function(srcpath) list.files(srcpath, pattern = "\\d{3}\\.rds$", full.names = T)
listPacc <- function(srcpath) list.files(srcpath, pattern = "\\d{3}\\.rds$", full.names = T)

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
    bgaccs=listACC,
    bgpccommunities=listPC,
    pertaccs=listPacc,
    score_mode=identity
  )
  parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
  parsed$options$help <- NULL
  result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
  # if(result$verbose) cat(format(result), file=stderr())
  result
}

resolve <- function(bgaccs, bgpccommunities, pertaccs, score_mode, verbose) {
    n <- min(length(bgaccs), length(bgpccommunities))
    ret <- rbindlist(mapply(function(bgaccincfn, bgpcincfn, paccfn, accinc){
      agg.dt <- readRDS(bgaccincfn)
      pc.dt <- readRDS(bgpcincfn)
      sl <- readRDS(paccfn)
      res <- perturbedPersistenceComms(sl, agg.dt, pc.dt, score_mode, verbose)
      res[, increment := accinc ]
      if(verbose) cat("finishing increment ", accinc, "; size ",dim(res),"\n",file=stderr())
      res
    }, bgaccincfn=bgaccs[1:n], bgpcincfn=bgpccommunities[1:n], paccfn=pertaccs[1:n], accinc=1:n, SIMPLIFY = F))
    # browser()
    ret
}

#for (i in 1:135) {
# args <- c("input/digest/background/15/30/censor/acc", "input/digest/background/15/30/censor/pc", "input/detection/high/hi/early/10/15/30/censor/010/acc.rds","-v")
saveRDS(
  with(
    parse_args(
#      args #c(sprintf("input/background-clusters/spin-glass/agg-15-30/%03d.rds",i),sprintf("input/background-clusters/spin-glass/pc-15-30/%03d.rds",i),sprintf("output/matched/mid/lo/late/10/001-covert-0/%03d-acc.rds",i), "-v")
    ),
    resolve(bgaccs, bgpccommunities, pertaccs, score_mode, verbose)
  ),
  pipe("cat","wb")
)
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
