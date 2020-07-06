library(changepoint)

args <- commandArgs(trailingOnly=T)

rdata <- args[[1]]
load(rdata)

cpt <- cpt.meanvar(
  data = as.matrix(t(obs[, -1, drop=F])),
  method = "BinSeg",
  test.stat = "Exponential",
  minseglen = 500
)

cat("sample\tchangepoints\tmeans\n")

for (i in 1:length(cpt)) {
  sid <- names(obs)[i + 1] # sample ID
  pts <- paste(paste(obs[cpt[[i]]@cpts, 1], collapse = ","), ",", sep = "") # detected changepoints
  avg <- paste(paste(round(1 / cpt[[i]]@param.est$rate, 2), collapse = ","), ",", sep = "") # detected means
  line <- paste(sid, pts, avg, sep = "\t")
  cat(paste(line, "\n"))
}
