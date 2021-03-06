#!/usr/bin/env Rscript

rm(list=ls())

require(data.table)
require(igraph)
require(parallel)

filelister <- function(srcpath) list.files(srcpath, full.names = T)

loadBase <- function(srcpath) list.files(srcpath, pattern = "\\d{3}\\.rds$", full.names = T)

readFore <- function(wh) {
  res <- fread(wh, col.names = c("user.a","user.b","location_id","start","end","type"))
  res[
    user.b < user.a,
    `:=`(user.b = user.a, user.a = user.b)
    ]
  res
}


parse_args <- function(argv = commandArgs(trailingOnly = T)) {
  parser <- optparse::OptionParser(
    usage = "usage: %prog path/to/precomputed path/to/simoutcomms.rds path/to/target",
    description = "accumulate scores prior to persistence calculations",
    option_list = list(
      optparse::make_option(
        c("--verbose","-v"),  action="store_true", default = FALSE,
        help="verbose?"
      )
    )
  )
  req_pos <- list(
    bgcommunities=loadBase, perturbedCommunities=readRDS, trim.dt = readRDS,
    intDays = as.integer, winDays = as.integer, scoremode = identity, tardir = identity
  )
  parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
  parsed$options$help <- NULL
  result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
  if(result$verbose) print(result)
  result
}

## need bgcommunities, foreground pattern, perturbation base

emptyscore <- data.table(community=integer(0), user.a = integer(0), user.b=integer(0), score=integer(0), increment=integer(0))

scorePerturbations <- function(pertComms, bgcommunities, scoring, trim.dt) {
  subset(rbindlist(mapply(function(inc, bginc) {
    perturb.dt <- pertComms[increment == inc, list(user_id, community)]
    res <- if (dim(perturb.dt)[1]) {
      perturbed <- perturb.dt[,unique(community)]
      back.dt <- readRDS(bginc)[community %in% perturbed]
      input.dt <- rbind(back.dt, perturb.dt)
      if (dim(input.dt)[1]) input.dt[, {
        tmp <- combn(user_id, 2)
        list(user.a = tmp[1,], user.b = tmp[2,], score = 1, increment = inc)
      },
      by=list(community)
      ] else emptyscore
    } else emptyscore
#     don't want raw.pairs, want cc / cu pairs
    if ((scoring == "bonus") & dim(res)[1]) {
      pairs.dt <- setkey(trim.dt[increment == inc], user.a, user.b)
      setkey(res, user.a, user.b)
      tars <- res[pairs.dt][!is.na(community), list(increment=T), keyby=list(user.a, user.b)]
      res[tars, score := score + 1]
    }
    res
  }, 1:length(bgcommunities), bgcommunities, SIMPLIFY = F))
  , select=-community)
}

storePerturbedScores <- function(obj, fn) {
  saveRDS(obj, fn)
  obj
}

accumPerturbedScores <- function(perturbedScores, tarfiles, discount, censor, n) {
  censor_score <- discount^censor
  incs <- 2:n
  Reduce(
    function(prev, currentN) {
      newres <- rbind(perturbedScores[increment == currentN, list(user.a, user.b, score, increment)], data.table::copy(prev)[, score := score*discount ])
      storePerturbedScores(
        newres[,list(score = sum(score), increment = currentN), keyby=list(user.a, user.b)][score > censor_score],
        tarfiles[currentN]
      )
    },
    incs,
    storePerturbedScores(perturbedScores[increment == 1, list(user.a, user.b, score, increment)], tarfiles[1])
  )
  TRUE
}

with(parse_args(),{
  perturbedScores <- scorePerturbations(perturbedCommunities, bgcommunities, scoremode, trim.dt)
  tars <- sprintf("%s/%03d.rds", tardir, 1:length(bgcommunities))
  if (verbose) cat(paste0(tars, collapse = "\n"), file=stderr())
  accumPerturbedScores(
    perturbedScores, tars,
    discount=0.9, censor=6, n=length(bgcommunities)
  )
  # accumPerturbedScores(perturbedScores, 0.9, 6, length(bgcommunities))
  # rbindlist(accumulatedPerturbs)
})
