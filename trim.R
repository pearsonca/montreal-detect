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
    location_lifetimes.dt=readRDS,
    cc.dt=readFore,
    cu.dt=readFore,
    intDays=as.integer,
    winDays=as.integer
  )
  parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
  parsed$options$help <- NULL
  result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
  if(result$verbose) print(result)
  result
}

slice <- function(dt, low, high) data.table::copy(
  dt[start < high*24*3600 & low*24*3600 < end]
)


saveRDS(with(parse_args(
#  c("input/raw-location-lifetimes.rds"
),{

  trimforegroundcc.dt <- location_lifetimes.dt[cc.dt][
    end > arrive & start < depart,
    list(user.a, user.b, location_id, start, end, type)
  ]
  
  foreground.dt <- setcolorder(rbind(cu.dt, trimforegroundcc.dt)[,
      list(user.a, user.b, start, end)
    ], c("user.a", "user.b", "start", "end")
  )
  
  res <- foreground.dt[user.a != user.b]
  
  ## need to assign increment based on raw.dt, intDays, winDays
  st <- raw.dt[1, floor(start/60/60/24)]
  mx <- raw.dt[1, ceiling((ceiling(end/60/60/24) - floor(start/60/60/24))/intDays)]
  setkey(rbindlist(lapply(1:mx, function(inc) {
    slice(res, st + inc*intDays-winDays, st + inc*intDays)[, increment := inc ]
  })), start, end)
}), pipe("cat","wb"))