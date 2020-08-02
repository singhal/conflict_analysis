# non-parametric bootstrapping to see if loci in support are weird

# locus data
files = list.files("~/Desktop/activeWork/conflict_analysis/", pattern = "loci_profiling", full.names = T)
d1 = read.csv(files[1], stringsAsFactors = F)
d2 = read.csv(files[2], stringsAsFactors = F)
d3 = read.csv(files[3], stringsAsFactors = F)
d4 = read.csv(files[4], stringsAsFactors = F)
d5 = full_join(d1, d2)
d6 = full_join(d5, d3)
d7 = full_join(d6, d4)
d7$comphet_RCFV = as.numeric(d7$comphet_RCFV)
keeper = c("branch_outliers", "comphet_RCFV", "GC", "heterozygosity",
         "max_PI",  "mean_resid","missing", "root_tip_var", 
         "saturation_cval", "tree_len", "length",
         "when_max_PI", "occupancy", "SH")
pos = c("when_max_PI", "occupancy", "SH", "length", "tree_len")
d8 = d7[, c("locus", keeper)]

dd4 = read.csv("dLNL.csv", stringsAsFactors = F)
clades = unique(dd4$clade)
dd5 = inner_join(dd4, d8, by="locus")
nsim = 1000
dLNL_min = 0
allres = vector('list', length(clades))

pdf(paste("~/Desktop/test.pdf", sep=""))
for (i in 1:length(clades)) {
  tmp = dd5 %>% filter(clade == clades[i], num_support < 2)
  # identify the #1 and #2 topology
  w = tmp %>% group_by(topology) %>% filter(dLNL > dLNL_min) %>% 
    summarise(count = n()) %>% arrange(desc(count))
  w1 = w$topology[1]
  w2 = w$topology[2]
  
  l1 = tmp %>% filter(topology == w1, dLNL > 2)
  l2 = tmp %>% filter(topology == w2, dLNL > 2)
  pop = bind_rows(l1, l2)
  
  create_random <- function(i) {
    # create random vector of topology
    d = bind_rows(l1, l2)
    d$topology = sample(c(rep(w1, nrow(l1)), rep(w2, nrow(l2))))
    means = d %>% group_by(topology) %>% select(keeper, "topology") %>% 
      summarize_all(mean, na.rm = T)
    means = means[1, keeper] - means[2, keeper]
    return(means)
  }
  
  res = lapply(1:1000, create_random)
  res2 = as.data.frame(do.call(rbind, res))

  actdiff = l1 %>% select(keeper) %>% summarize_all(mean, na.rm = T) - 
    l2 %>% select(keeper) %>% summarize_all(mean, na.rm = T)
  
  res3 = as.data.frame(matrix(NA, nrow = length(keeper), ncol = 5), stringsAsFactors = F)
  names(res3) = c("topology", "mean", 
                  as.character(w1), as.character(w2), 
                  "significance")  
  res3$topology = clades[i]
  rownames(res3) = keeper
  
  
  par(mfrow=c(3, 5), mar=c(3, 3, 2, 0))
  for (j in 1:length(keeper)) {
    res3[keeper[j], "mean"] = mean(tmp[, keeper[j]], na.rm = T)
    res3[keeper[j], w1] = mean(l1[, keeper[j]], na.rm = T)
    res3[keeper[j], w2] = mean(l2[, keeper[j]], na.rm = T)
    
    xvals = res2[, keeper[j]]
    sig = sum(xvals > abs(pull(actdiff[keeper[j] ]))) / 1000
    res3[keeper[j], "significance"] = sig
    
    xrange = range(xvals, actdiff[keeper[j] ])
    hist(xvals, main = keeper[j], las = 1, xlab="", 
         breaks=30, col="gray", border="gray", xlim = xrange)
    if (res3[keeper[j], "significance"]  < 0.05) {
      abline(v = actdiff[keeper[j] ], col = "red") 
    } else {
      abline(v = actdiff[keeper[j] ], col = "black") 
    }
  }
  allres[[i]] = res3
  cat(i)
}
dev.off()

simplify_df <- function(x) {
  names(x) <- c("topology", "mean", "t1", "t2", "significance")
  x$metric = rownames(x)
  return(x)
}
allres2 = allres
allres2 = lapply(allres2, simplify_df)
allres2 = do.call("rbind", allres2)
allres2$pos = 1
allres2[allres2$metric %in% pos, "pos"] = -1
allres3 = allres2 %>% mutate(diff2 = (t1 - t2) * pos, diff = (t1 - t2))

rounds = c("mean", "t1", "t2", "diff", "diff2")
for (i in 1:length(rounds)) {
  allres3[, rounds[i]] = signif( allres3[, rounds[i]], 3 )
#   if (mean(allres3[, rounds[i]]) < 0) {
#     allres3[, rounds[i]] = formatC(allres3[, rounds[i]], format = "e", digits = 2)
#   }
}

allres3 = allres3 %>% filter(significance < 0.01)


clades= c("anniellidae", "anomalepidadae", "bolyeridae", 
          "cylindro_uropelt", "dibamidae", "diplodactylidae",
          "eublepharid", "gecko_skink", "gymnophthalmidae",
          "homalopsidae", "iguanids", "lanthanotidae",
          "rhineuridae", "skinks", "snakes", "typhlopidae",
          "xenosauridae")
cnames = c("Anniellidae", "Anomalepidadae", "Bolyeridae", 
           "Cylindrophiidae & Uropeltidae", "Dibamidae", "Diplodactylidae",
           "Eublepharidae", "Gekkota & Scincoidea", "Gymnophthalmidae",
           "Homalopsidae", "Iguania", "Lanthanotidae",
           "Rhineuridae", "Scincoidea", "Serpentes", "Typhlopidae",
           "Xenosauridae")
names(cnames) = clades
metrics = c("branch_outliers", "comphet_RCFV", "GC", "heterozygosity",
         "length", "max_PI", "when_max_PI", "mean_resid",
         "missing", "occupancy", "root_tip_var", "saturation_cval",
         "SH", "tree_len")
mnames = c("branch outliers", "comp. het. RCFV", "GC", "heterozygosity",
          "length", "max PI", "when max PI", "mean residuals",
          "missing", "occupancy", "root-tip variance", "saturation C-value",
          "SH", "tree length")
names(mnames) = metrics

allres3$topology = cnames[allres3$topology]
allres3$metric = mnames[allres3$metric]
allres3 = allres3[ , c("topology", "metric", "significance", "mean", "t1", "t2", "diff", "diff2")]
write.csv(allres3, paste(figdir, "loci_dLNL_", dLNL_min, ".csv", sep=""), 
          row.names = F)
