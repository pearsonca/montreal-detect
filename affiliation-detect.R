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

foregroundSnapshots <- readRDS("output/matched/mid/lo/late/10/001-covert-0-base.rds")
foregroundN <- 10
n <- 61
res <- data.table(increment=1:n,
  snapFPR=numeric(n), snapTPR=numeric(n),
  pcFPR=numeric(n), pcTPR=numeric(n),
  key="increment"
)

for (i in 1:n) {
  args <- c(
    sprintf("input/background-clusters/spin-glass/base-15-30/%03d.rds",i),
    sprintf("input/background-clusters/spin-glass/pc-15-30/%03d.rds",i),
    sprintf("output/matched/mid/lo/late/10/001-covert-0/%03d.rds",i)
  )

  backgroundSnapshot <- readRDS(args[1])
  foregroundSnapshot <- foregroundSnapshots[increment == i, list(user_id, community)]
  counts <- rbind(backgroundSnapshot, foregroundSnapshot)[,
    list(total=.N, bg=sum(user_id >= 0), fg=sum(user_id<0)),
    keyby=community
  ]
  
  res[increment == i,
    snapFPR:=counts[total <= 30, sum(bg)]/max(counts[,sum(bg)],1),
  ]
  res[increment == i,
    snapTPR:=counts[total <= 30, sum(fg)]/foregroundN
  ]
  
  backgroundPersistent <- readRDS(args[2])
  foregroundPersistent <- readRDS(args[3])

  counts <- rbind(backgroundPersistent, foregroundPersistent)[,
    list(total=.N, bg=sum(user_id >= 0), fg=sum(user_id<0)),
    keyby=community
  ]
  
  res[increment == i,
    pcFPR := counts[total <= 30, sum(bg)]/counts[,sum(bg)]
  ]
  res[increment == i,
    pcTPR := counts[total <= 30, sum(fg)]/foregroundN
  ]
  
}

require(ggplot2)
require(reshape2)

pltres <- melt.data.table(res, id.vars="increment", variable.name = "measure", value.name = "rate")

pltres[,
  stage := factor(gsub("[FT]PR", "", measure))
][,
  outcome := factor(gsub(".+([FT]PR)","\\1", measure))
]

ggplot(pltres) + theme_bw() +
  aes(x=increment/2, y=rate, color=outcome, linetype=stage) + geom_line() +
  labs(x="increment") + scale_color_manual(values=c(FPR='red',TPR='blue')) +
  scale_linetype_manual(values=c(pc="solid",snap="dashed"))