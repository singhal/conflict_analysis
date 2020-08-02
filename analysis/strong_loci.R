files = list.files("~/Dropbox (Personal)/conflict_analysis/", pattern = "loci_profiling", full.names = T)
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

d = d7[, c("locus", keeper)]
# pca = prcomp(d %>% select(-locus), scale. = T, center = T)
d = d[complete.cases(d), ]
for (i in 1:length(keeper)) {
  d[ ,paste(keeper[i], "_res", sep="")] = lm(d[, keeper[i]] ~ d$tree_len)$residuals
}
keeper2 = paste0(keeper, "_res")

t = read_csv("dLNL.csv") %>% group_by(locus) %>%
  summarize(mean_dLNL = mean(dLNL)) %>% ungroup()

# t = read_csv("dLNL.csv") %>% filter(clade == "gecko_skink")
# t = t %>% filter(topology == "gecko_skink_sister") %>% 
#   mutate(dLNL_strong = ifelse(dLNL > 2, TRUE, FALSE))
t2 = left_join(t, d)
# t3 = t2 %>% select(locus, dLNL_strong, keeper2) %>% tidyr::gather(metric, value, -locus, -dLNL_strong)
# View(t3 %>% group_by(metric, dLNL_strong) %>% 
#       summarize(mean(value, na.rm = T)) %>% filter(metric %in% c("SH_res", "length_res", # "max_PI_res", "tree_len_res", "when_max_PI_res")))
# ggplot(t3, aes(dLNL_strong, value)) + geom_boxplot() + facet_wrap(.~metric, scales = # "free")

t2 = as.data.frame(t2)
t2 = t2[complete.cases(t2$mean_dLNL), ]
t2 = t2[complete.cases(t2$tree_len), ]
t2$lnl_res = lm(t2$mean_dLNL ~ t2$tree_len)$residuals
for (i in 1:length(keeper)) {
  xx = cor.test(t2$lnl_res, t2[, keeper2[i]], method = "spearman")
  cat(keeper[i], round(xx$estimate, 4), round(xx$p.value, 4), "\n", sep = "\t")
}