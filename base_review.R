# review the base detection for some combination of params, across samples

require(ggplot2)
require(data.table)

args <- commandArgs(trailingOnly = T)
sampDir <- args[1]
backgroundDir <- args[2]
tar <- args[3]

foregroundN <- as.integer(sub(".+/(\\d+)/\\d+/\\d+$","\\1",sampDir))

bgs <- list.files(backgroundDir, pattern="\\d{3}.rds", full.names = T)
background <- rbindlist(lapply(bgs, function(fn) {
  inc <- as.integer(sub(".+/(\\d+)\\.rds", "\\1", fn))
  res <- readRDS(fn)
  res[, increment := inc ]
  res
}))

fgs <- list.files(sampDir, pattern = "base.rds", full.names = T, recursive = T)

saveRDS(rbindlist(lapply(fgs, function(fn) {
  # parse out covert size from fn
  samp <- as.integer(sub(".+/(\\d+)/base.rds", "\\1", fn))
  foreground <- try(readRDS(fn))
  if (class(foreground)=="try-error") stop(fn)

  counts <- rbind(background, foreground)[,
    list(total=.N, bg=sum(user_id >= 0), fg=sum(user_id < 0)),
    keyby=list(community, increment)
  ]

#   snap <- counts[foreground][user_id < 0]
#   
#   jn <- snap[,
#     list(found = total <= 30),
#     keyby=list(increment, user_id)
#   ]
  
  res <- counts[, list(
    snapFPR = sum(bg[total <= 30])/max(sum(bg)),
    snapTPR = sum(fg[total <= 30])/foregroundN
  ), keyby=increment]
  
  res[, sample_id := samp ]
  
  res
})), pipe("cat","wb"))

do.call(c(warnings(), file=stderr()))