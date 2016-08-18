# review the base detection for some combination of params, across samples

require(ggplot2)
require(data.table)

args <- commandArgs(trailingOnly = T)
sampDir <- args[1]
backgroundDir <- args[2]

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
  if (class(foreground)[1]=="try-error") stop(fn)

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

# wrn <- warnings()
# cat(paste(names(wrn), unlist(wrn), collapse = "\n"), file=stderr())
#
# plotres <- thing[,{ tar <- which.max(snapTPR); list(`peak detection TPR`=snapTPR[tar], increment=increment[tar])},by=sample_id]
# ggplot(plotres) + aes(x=increment, y=`peak detection TPR`, label=sample_id) + geom_text() + coord_cartesian(ylim=c(0,1)) + theme_bw()
#
# plotres2 <- thing[snapTPR != 0, c(list(`summary TPR`=sum(snapTPR)),as.list(quantile(snapTPR))),by=sample_id]
# ggplot(plotres2) + aes(y=`summary TPR`, x=sample_id) + geom_point() + theme_bw()
# setkey(plotres2,`100%`,`75%`,`50%`,`25%`,`0%`)
# plotres2[,list(new_sample_id=.GRP),by=list(`100%`,`75%`,`50%`,`25%`,`0%`,sample_id)]
# plotres3 <- melt(, measure.vars = c('0%','25%','50%','75%','100%'), variable.name = "quantile")
# ggplot(plotres3) + aes(y=value, x=sample_id, color=quantile) + geom_line() + theme_bw()