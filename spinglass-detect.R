#!/usr/bin/env Rscript

require(data.table)

args <- c("input/raw-input.rds","input/background-clusters/spin-glass/pc-30-15", "output/matched/mid/lo/late/10/001-covert-0")

raw.dt <- readRDS(args[1])
backgroundsrc <- args[2]
interval <- as.integer(sub(".+-(\\d+)-\\d+$","\\1", backgroundsrc))
window <- as.integer(sub(".+-\\d+-(\\d+)$","\\1", backgroundsrc))
# this will be something like XYZ(15|30)-(15|30), where first number is window, second is interval
foregroundpat <- args[3]
foregroundcc.dt <- fread(paste0(foregroundpat,"-cc.csv"))
foregroundcu.dt <- fread(paste0(foregroundpat,"-cu.csv"))
# parse foreground according to interval

