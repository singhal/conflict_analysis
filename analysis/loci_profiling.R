library(dplyr)
library(corrplot)
library(cowplot)
library(tidyr)

files = list.files("~/Dropbox (Personal)/conflict_analysis/", pattern = "loci_profiling", full.names = T)
d1 = read.csv(files[1], stringsAsFactors = F)
d2 = read.csv(files[2], stringsAsFactors = F)
d3 = read.csv(files[3], stringsAsFactors = F)
d4 = read.csv(files[4], stringsAsFactors = F)

d5 = full_join(d1, d2)
d6 = full_join(d5, d3)
d7 = full_join(d6, d4)
d7$comphet_RCFV = as.numeric(d7$comphet_RCFV)

d8 = d7 %>% gather("metric", "value", -locus)
d8 = d8 %>% filter(metric %in% keep)
d8$type = rep("gene", nrow(d8))
d8[grep("uce", d8$locus), "type"] = "uce"
d8[grep("AHE", d8$locus), "type"] = "AHE"

d9 = split(d8, d8$metric) 
aov10 = lapply(d9, function(x) {return(aov(x$value ~ x$type)) })

# pvals 
pvals1 = unlist(lapply(aov10, function(x) {return(summary(x)[[1]][["Pr(>F)"]][1])}))
names(pvals1) = keep1
# tukeyHSD
thsd = lapply(aov10, function(x) { return(TukeyHSD(x))})
names(thsd) = keep1
thsd[names(pvals1[pvals1 < (0.05 / 14)])]

minmax = d8 %>% group_by(metric) %>% summarise(min01 = quantile(value, 0.01, na.rm = T), 
                                      max99 = quantile(value, 0.99, na.rm = T))
d8out = inner_join(d8, minmax) %>% filter(value > min01 & value < max99)
names(keep1) = keep
d8out$metric = keep1[d8out$metric]
xx = ggplot(d8out, aes(type, value))  +
  geom_point(color=alpha("gray", 0.05), position = position_jitter()) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) +
  facet_wrap(~metric, scales = "free", ncol = 5, nrow = 3)
save_plot(paste(figdir, "locus_metrics_by_type.png", sep=""),
          xx, ncol = 5, nrow = 3, base_height = 1.8, base_width = 3)

M = cor(d7[ , 2:18], use="pairwise.complete.obs")
corrplot(M)

cols = colnames(M)
for (i in 1:length(cols)) {
  for (j in 1:length(cols)) {
    val = M[i, j]
    if (abs(val) >= 0.7 && i != j) {
      cat(cols[i], " ", cols[j], " ", round(val, 3), "\n")
    }
  }
}


#### this suggests to drop
# PICS, tree len, mean resid

keep = c("branch_outliers", "comphet_RCFV", "GC", "heterozygosity",
         "length", "max_PI", "when_max_PI", "mean_resid",
         "missing", "occupancy", "root_tip_var", "saturation_cval",
         "SH", "tree_len")
keep1 = c("branch outliers", "comp. het. RCFV", "GC", "heterozygosity",
          "length", "max PI", "when max PI", "mean residuals",
          "missing", "occupancy", "root-tip variance", "saturation C-value",
          "SH", "tree length")
d8 = d7[, keep]
rownames(d8) = d1$locus

d8$type = rep("gene", nrow(d8))
d8[grep("uce", rownames(d8)), "type"] = "uce"
d8[grep("AHE", rownames(d8)), "type"] = "AHE"

type = c("max", "max", "max", "max",
         "min", "max", "min", "max",
         "max", "min", "max", "max",
         "min", "max")
ditch = vector('list', length(keep))
chires = vector('list', length(keep))
# pdf("~/Desktop/activeWork/conflict_analysis/prelim_figures/molecular_evolution.pdf", width = 4, height = 3)
for (i in 1:length(keep)) {
  filter = type[i]
  vec = d8[, keep[i]]
  if (filter == 'max') {
    val = quantile(vec, 0.98, na.rm = T)
    ditch[[i]] = rownames(d8[which(d8[, keep[i]] >= val), ])
  } else {
    val = quantile(vec, 0.02, na.rm = T)
    ditch[[i]] = rownames(d8[which(d8[, keep[i]] <= val), ])
  }
  # do chisquare test
  tmp = d8[complete.cases(keep[i]), ]
  tmp$ditch = FALSE
  tmp[ditch[[i]], "ditch"] = TRUE
  chires[[i]] = chisq.test(table(tmp$ditch, tmp$type))
#  a = ggplot(d8, aes_string(keep[i])) + geom_histogram(bins = 100) + geom_vline(xintercept = val, color = "red")
#  print(a)
}
# dev.off()
# chisquare pvals
pvals = unlist(lapply(chires, function(x){return(x$p.value)}))
names(pvals) = keep1
# diff in expectations
diffexp = lapply(chires, function(x) { return(x$observed - x$expected) })
names(diffexp) = keep1
# significant differences
diffexp[names(pvals[pvals < (0.05 / 14)])]

to_drop = data.frame(loci = unique(unlist(ditch)), stringsAsFactors = F)
write_csv(to_drop, "~/Desktop/activeWork/conflict_analysis/molevol_outliers.csv")
