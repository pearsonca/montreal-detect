#!/usr/bin/env Rscript

rm(list=ls())

require(data.table)
require(igraph)
require(parallel)

args <- c("input/raw-pairs.rds", "input/location-lifetimes.rds", "input/background-clusters/spin-glass/pc-30-30", "output/matched/mid/lo/late/10/001-covert-0")

raw.dt <- readRDS(args[1])
maxuid <- raw.dt[,max(user.a, user.b)]

loclifes.dt <- readRDS(args[2])

backgroundsrc <- args[3]
bgcommunities <- list.files(sub("pc","base",args[3]), full.names = T)
interval <- as.integer(sub(".+-(\\d+)-\\d+$","\\1", backgroundsrc))
window <- as.integer(sub(".+-\\d+-(\\d+)$","\\1", backgroundsrc))

# this will be something like XYZ(15|30)-(15|30), where first number is window, second is interval
foregroundpat <- args[4]
foregroundcc.dt <- setkey(fread(paste0(foregroundpat,"-cc.csv"), col.names = c("user.a","user.b","location_id","start","end","type")), location_id)

# trim this by location appearances
trimforegroundcc.dt <- loclifes.dt[foregroundcc.dt][end > arrive & start < depart, list(user.a, user.b, location_id, start, end, type)]

# no trim required - intersection with real users
foregroundcu.dt <- fread(paste0(foregroundpat,"-cu.csv"), col.names = c("user.a","user.b","location_id","start","end","type"))


foreground.dt <- rbind(foregroundcu.dt, trimforegroundcc.dt)[,list(user.a,user.b,start,end)]
setcolorder(foreground.dt,c("user.a", "user.b", "start", "end"))

temp <- rbind(foreground.dt, raw.dt)
temp[
  user.b < user.a,
  `:=`(user.b = user.a, user.a = user.b)
]

ref.dt <- setkey(temp[user.a != user.b], start, end)

source("../montreal-digest/buildStore.R")

slice <- function(dt, low, high) relabeller(
  dt[start < high*24*3600 & low*24*3600 < end,
     list(score=.N, covert.a=user.a < 0, covert.b=user.b < 0),
     by=list(user.a, user.b)]
)

# parse foreground according to interval

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

resolve <- function(
  base.dt, intDays, winDays, #outputdir, crs,
  mxint=NA, verbose=F, ...
) {
  st = base.dt[1, floor(start/60/60/24)]
  n <- min(ceiling(base.dt[,max(end)/60/60/24 - st]/intDays), mxint, na.rm = TRUE)
  targets <- 1:n
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

perturbedCommunities <- resolve(ref.dt, interval, window)

# refUsers <- data.table(user_id=rep(-1:-10, times=68), increment=rep(1:68, each=10), key="user_id")
# plot_ref <- merge(refUsers, perturbedCommunities, by=c("increment","user_id"), all=T)
# ggplot(plot_ref[user_id < 0]) + facet_grid(user_id ~ .) + aes(x=increment, y=(community==-1), color=(community==-1)) + geom_point()

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
      ] else data.table(community=integer(0), user.a = integer(0), user.b=integer(0), score=integer(0), increment=integer(0))
    } else data.table(community=integer(0), user.a = integer(0), user.b=integer(0), score=integer(0), increment=integer(0))
  }, 1:length(bgcommunities), bgcommunities, SIMPLIFY = F))
  , select=-community)
}

perturbedScores <- scorePerturbations(perturbedCommunities, bgcommunities)

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

accumulatedPerturbs <- accumPerturbedScores(perturbedScores, 0.9, 6, length(bgcommunities))

# have a directory from base file
# save 001, 002, etc for each list element

mapply(function(accPerturbs, inc) {
  # browser()
  saveRDS(accPerturbs, sprintf("%s/%03d-acc.rds", foregroundpat, inc))
}, accumulatedPerturbs, 1:length(accumulatedPerturbs))