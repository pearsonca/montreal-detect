rm(list=ls())

# args are (1) base background (for increment),
# (2) base foreground (all)
# (3) pc background (for increment)
# (4) pc foreground (for increment)
# (5) target

# parse_args <- function(argv = commandArgs(trailingOnly = T)) {
#   parser <- optparse::OptionParser(
#     usage = "usage: %prog path/to/input-useratob-scores/ path/to/accumulated-useratob-scores/",
#     description = "convert (community X, user.a in community X, user.b in community X, 1) at interval k, to (...) cumulated up to interval k.",
#     option_list = list(
#       optparse::make_option(
#         c("--verbose","-v"),  action="store_true", default = FALSE,
#         help="verbose?"
#       )
#     )
#   )
#   req_pos <- list(inputfiles=filelister, outputdir=identity)
#   parsed <- optparse::parse_args(parser, argv, positional_arguments = length(req_pos))
#   parsed$options$help <- NULL
#   result <- c(mapply(function(f,c) f(c), req_pos, parsed$args, SIMPLIFY = F), parsed$options)
#   result$storeres <- function(dt, was) {
#     saveRDS(dt, sub(parsed$args[1], parsed$args[2], was))
#     dt
#   }
#   
#   if(result$verbose) print(result)
#   result
# }

# given community affiliations for background + foreground
#  for both regular and persistence communities X intervals
#   compute # of covert members in small comms vs large comms
#   compute # of non-covert members in small vs large comms
#   compute TPR, FPR

lf <- "mid"
pwr <- "lo"
pk <- "late"

# lf <- "mid"
# pwr <- "med"
# pk <- "middle"

foregroundSnapshots <- readRDS(sprintf("output/matched/%s/%s/%s/10/001-covert-0-base.rds",lf,pwr,pk))
foregroundN <- 10
n <- 135
res <- data.table(increment=1:n,
  snapFPR=numeric(n), snapTPR=numeric(n), snapAltTPR=numeric(n),
  pcFPR=numeric(n), pcTPR=numeric(n), pcAltTPR=numeric(n),
  key="increment"
)

covertStates <- c(abs="absent",big="present in large community",lit="present in small community")

bgcomp_user_ids <- readRDS("input/user.RData")[
  (lifetime_main == "mid" & pwr_main == "lo" & peak_main == "late"), user_id
  ]

altfgN <- min(foregroundN, length(bgcomp_user_ids))
bgsmp <- sample(bgcomp_user_ids, altfgN)

covertTimeLine <- data.table(
  increment=rep(1:n, each=foregroundN+altfgN), user_id=rep(c(-(1:foregroundN), bgsmp), times=n),
  outcome=factor(covertStates["abs"], levels=covertStates), key=c("increment","user_id")
)

covertPersistence <- data.table(
  increment=rep(1:n, each=foregroundN+altfgN), user_id=rep(c(-(1:foregroundN), bgsmp), times=n),
  outcome=factor(covertStates["abs"], levels=covertStates), key=c("increment","user_id")
)

covertSnapAccumulation <- data.table(
  increment=rep(1:n, each=foregroundN+altfgN), user_id=rep(c(-(1:foregroundN), bgsmp), times=n),
  count = 0, key=c("increment","user_id")
)

covertPersistAccumulation <- data.table(
  increment=rep(1:n, each=foregroundN+altfgN), user_id=rep(c(-(1:foregroundN), bgsmp), times=n),
  count = 0, key=c("increment","user_id")
)

## compute some random background comparison on mid/lo/late users

## compute 
for (i in 1:n) {
  args <- c(
    sprintf("input/background-clusters/spin-glass/base-15-30/%03d.rds",i),
    sprintf("input/background-clusters/spin-glass/pc-15-30/%03d.rds",i),
    sprintf("output/matched/mid/lo/late/10/001-covert-0/%03d.rds",i)
  )

  backgroundSnapshot <- readRDS(args[1])
  backgroundSampshot <- backgroundSnapshot[user_id %in% bgsmp, user_id, keyby=community]
  foregroundSnapshot <- foregroundSnapshots[increment == i, user_id, keyby=community]
  
  counts <- rbind(backgroundSnapshot, foregroundSnapshot)[,
    list(total=.N, bg=sum(user_id >= 0), fg=sum(user_id<0), altfg=length(intersect(user_id, bgsmp))),
    keyby=community
  ]
  
  snap <- rbind(
    counts[foregroundSnapshot][user_id < 0], # sometimes regular users are added due to covert perturbation
    counts[backgroundSampshot]
  )
  
  if (dim(snap)[1]) {
    jn <- setkey(snap[,
      list(increment = i, found = total <= 30),
      by=user_id
    ], increment, user_id)
    covertSnapAccumulation[jn, count := found + 0]
    covertTimeLine[jn, outcome := ifelse(found, covertStates["lit"], covertStates["big"])]
  }
  
  res[increment == i,
    `:=`(
      snapFPR = counts[total <= 30, sum(bg)]/max(counts[,sum(bg)],1),
      snapTPR = counts[total <= 30, sum(fg)]/foregroundN,
      snapAltTPR = counts[total <= 30, sum(altfg)]/altfgN
    )
  ]

  if (!file.exists(args[2])) {
    cat("missing background",i,"\n")
    next()
  }
  backgroundPersistent <- readRDS(args[2])
  tres <- try(readRDS(args[3]))
  foregroundPersistent <- if (class(tres)[1] == "try-error") {
    cat("missing",i,"\n")
    data.table(user_id=integer(), community=integer(), key="community")
  } else setkey(tres, community)

  counts <- rbind(backgroundPersistent, foregroundPersistent)[,
    list(total=.N, bg=sum(user_id >= 0), fg=sum(user_id<0), altfg=length(intersect(user_id, bgsmp))),
    keyby=community
  ]

  backgroundSampsist <- backgroundPersistent[user_id %in% bgsmp, user_id, keyby=community]
  
  snap <- rbind(
    counts[foregroundPersistent][user_id < 0], # sometimes regular users are added due to covert perturbation
    counts[backgroundSampsist]
  )
    
  #snap <- counts[foregroundPersistent][user_id < 0]

  if (dim(snap)[1]) {
    jn <- setkey(snap[,
      list(increment = i, found = total <= 30),
      by=user_id
    ], increment, user_id)
    covertPersistAccumulation[jn, count := found + 0]
    covertPersistence[jn, outcome := ifelse(found, covertStates["lit"], covertStates["big"])]
  }
  
  res[increment == i,
    `:=`(
      pcFPR = counts[total <= 30, sum(bg)]/counts[,sum(bg)],
      pcTPR = counts[total <= 30, sum(fg)]/foregroundN,
      pcAltTPR = counts[total <= 30, sum(altfg)]/altfgN
    )
  ]
  
}

snapAcc <- covertSnapAccumulation[,list(increment, cs=cumsum(count)), by=user_id]
persAcc <- covertPersistAccumulation[,list(increment, cs=cumsum(count)), by=user_id]

ggplot(snapAcc) + aes(y=cs,x=increment, color=factor(user_id)) + geom_step() + theme_bw()
ggplot(persAcc) + aes(y=cs,x=increment, color=factor(user_id)) + geom_step() + theme_bw()

require(ggplot2)
require(reshape2)

pltres <- melt.data.table(res, id.vars="increment", variable.name = "measure", value.name = "rate")

pltres[,
  `community analysis` := factor(c(pc="persistence", snap="snapshot", snapAlt="altSnap", pcAlt="altpc")[gsub("[FT]PR", "", measure)])
][,
  outcome := factor(gsub(".+([FT]PR)","\\1", measure))
]

ggplot(pltres) + theme_bw() + theme(panel.border=element_blank()) +
  aes(x=increment, y=rate, color=outcome, linetype=`community analysis`) + geom_line() +
  scale_color_manual(values=c(FPR='red',TPR='blue')) + labs(y="TPR & FPR") +
  scale_linetype_manual(values=c(persistence="solid",snapshot="dashed", altpc="dotted", altSnap="dotdash")) +
  ggtitle(sprintf("For %s %s %s", lf, pwr, pk))

ggplot(covertTimeLine) + theme_bw() + aes(x=increment, y=factor(user_id), fill=outcome) +
  geom_raster() + scale_fill_manual(values=c(absent="red",`present in large community`="yellow",`present in small community`="green")) +
  ggtitle(sprintf("For %s %s %s", lf, pwr, pk))
  
ggplot(covertPersistence) + theme_bw() + aes(x=increment, y=factor(user_id), fill=outcome) +
  geom_raster() + scale_fill_manual(values=c(absent="red",`present in large community`="yellow",`present in small community`="green")) +
  ggtitle(sprintf("For %s %s %s", lf, pwr, pk))
