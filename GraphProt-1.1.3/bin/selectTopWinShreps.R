library(plyr)
# parse command line arguments
args <- commandArgs(trailingOnly=T)
nt_margins <- args[1]
# get data
d <- read.table(nt_margins, col.names=c('shrep','position','score','min','max','mean','median','sum'))
# get best scored windows per shrep
d_top_wins <- ddply(d, .(shrep), function(df) {subset(df, mean==max(mean))})
# get best scored womdow farthest to the right per shrep
d_top_win <- ddply(d_top_wins, .(shrep), function(df) {subset(df, position==min(position))})
# write selection to stdout
write.table(d_top_win, '', quote=F, col.names=F, row.names=F, sep="\t")
