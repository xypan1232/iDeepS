library(stats)

args <- commandArgs(trailingOnly=T)

has_percentile <- as.numeric(args[1])
infile <- args[2]
outfile <- args[3]

# get data
d <- read.table(infile, col.names=c('seqid', 'seqpos', 'mean_margin'))

# calculate threshold
has_percentile_prob <- has_percentile / 100
has_threshold <- quantile(d[,'mean_margin'], has_percentile_prob)[1]

# subset threshold
dthresh <- subset(d, mean_margin > has_threshold)

# write data
write.table(dthresh, file=outfile, quote=F, sep = "\t", row.names=F, col.names=F)
