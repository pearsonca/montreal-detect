#!/usr/bin/env Rscript

require(data.table)

args <- c("input/raw-pairs.rds","input/location-lifetimes.rds","input/background-clusters/spin-glass/pc-30-15", "output/matched/mid/lo/late/10/001-covert-0")

raw.dt <- readRDS(args[1])
maxuid <- raw.dt[,max(user.a, user.b)]

loclifes.dt <- readRDS(args[2])

backgroundsrc <- args[3]
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

ref.dt <- setkey(temp, start, end)

source("../montreal-digest/buildStore.R")

slice <- function(dt, low, high) relabeller(
  dt[start < high*24*3600 & low*24*3600 < end,
     list(score=.N, covert.a=user.a < 0, covert.b=user.b < 0),
     by=list(user.a, user.b)]
)

# parse foreground according to interval

graphPartition <- function(res, referenceCommunities, coverts, ulim=60, verbose=F) {
  gg <- graph(t(res[,list(user.a, user.b)]), directed=F)
  E(gg)$weight <- res$score

  comps <- components(gg)
  
  leftovers <- which(comps$csize > ulim)
  completeCommunities <- (1:comps$no)[-leftovers] # components to treat as their own communities
  
  if (verbose) {
    cat(sprintf("found %d components of size <= %d\n", length(completeCommunities), ulim))
    cat(sprintf("found %d components of size > %d\n", length(leftovers), ulim))
  }
  
  base <- if (length(completeCommunities)) {
    newuids <- which(comps$membership %in% completeCommunities)
    commap  <- rep.int(NA, max(completeCommunities))
    commap[completeCommunities] <- 1:length(completeCommunities)
    data.table(
      new_user_id=newuids,
      community=commap[comps$membership[newuids]],
      key="new_user_id"
    ) #  mp[newuids, list(user_id, community=newcoms)]
  } else data.table(new_user_id=integer(0), community=integer(0))
  
  list(gg=gg, comps=comps, leftovers=leftovers, base=base)
}

resolve <- function(
  base.dt, intDays, winDays, #outputdir, crs,
  mxint=NA, verbose, ...
) {
  st = base.dt[1, floor(start/60/60/24)]
  n <- min(ceiling(base.dt[,max(end)/60/60/24 - st]/intDays), mxint, na.rm = TRUE)
  targets <- 1:n
  lapply(targets, function(inc) {
    with(slice(base.dt, st + inc*intDays-winDays, st + inc*intDays), {
      coverts <- res[,unique(c(user.a[covert.a],user.b[covert.b]))]
      if (length(coverts)) {
        cat("something to do ", inc, length(coverts), mp[coverts, user_id], "\n") 
        #graphPartition(res, coverts)
        
      } else {
 
      }
    })
  })
  # completed <- as.integer(gsub(".rds","", list.files(outputdir, "rds") ))
  # want <- targets[c(-completed,-(n+1))]
  # #  system.time(
  # mclapply(want, function(inc) with(slice(base.dt, st + inc*intDays-winDays, st + inc*intDays), {
  #   store <- if (dim(res)[1] == 0) emptygraph else buildStore(res)
  #   resfile <- sprintf("%s/%03d.rds", outputdir, inc)
  #   if (verbose) cat("finishing", resfile,"\n")
  #   saveRDS(
  #     originalUserIDs(store, mp),
  #     resfile
  #   )
  # }), mc.cores = crs, mc.allow.recursive = F)
}

resolve(ref.dt, interval, window)
