#!/usr/bin/env Rscript

rm(list=ls())

# (1) raw background
# (2) location lifetimes to censor covert
# (3) reference persistence communities
# (4) reference simulation source

require(data.table)
require(igraph)
require(parallel)

filelister <- function(srcpath) list.files(srcpath, full.names = T)

parse_args <- function(argv = commandArgs(trailingOnly = T)) {
  parser <- optparse::OptionParser(
    usage = "usage: %prog path/to/precomputed path/to/simoutcomms.rds path/to/target",
    description = "compute snapshot communities + persistence scores for a simulation",
    option_list = list(
      optparse::make_option(
        c("--verbose","-v"),  action="store_true", default = FALSE,
        help="verbose?"
      )
    )
  )
  req_pos <- list(
    bgcommunities=filelister,
    perturbedCommunities=readRDS,
    outpath=identity
  )
  parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
  parsed$options$help <- NULL
  result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
  if(result$verbose) print(result)
  result
}

## need bgcommunities, foreground pattern, perturbation base

emptyscore <- data.table(community=integer(0), user.a = integer(0), user.b=integer(0), score=integer(0), increment=integer(0))

scorePerturbations <- function(pertComms, bgcommunities) {
  subset(rbindlist(mapply(function(inc, bginc) {
    perturb.dt <- pertComms[increment == inc, list(user_id, community)]
    if (dim(perturb.dt)[1]) {
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
  }, 1:length(bgcommunities), bgcommunities, SIMPLIFY = F))
  , select=-community)
}

accumPerturbedScores <- function(perturbedScores, discount, censor, n) {
  censor_score <- discount^censor
  Reduce(
    function(prev, currentN) {
      newres <- rbind(perturbedScores[increment == currentN, list(user.a, user.b, score)], data.table::copy(prev)[, score := score*discount ])
      newres[,list(score = sum(score)), keyby=list(user.a, user.b)][score > censor_score]
    },
    2:n,
    perturbedScores[increment == 1, list(user.a, user.b, score)],
    accumulate = T
  )
}

with(parse_args(
#  c("input/background-clusters/spin-glass/base-15-30", "output/matched/mid/lo/late/10/001-covert-0-base.rds", "output/matched/mid/lo/late/10/001-covert-0")
),{
  perturbedScores <- scorePerturbations(perturbedCommunities, bgcommunities)
  accumulatedPerturbs <- accumPerturbedScores(perturbedScores, 0.9, 6, length(bgcommunities))
  mapply(function(accPerturbs, inc) {
    # browser()
    saveRDS(accPerturbs, sprintf("%s/%03d-acc.rds", outpath, inc))
  }, accumulatedPerturbs, 1:length(accumulatedPerturbs))
})
