library(readr)
library(tidyr)
library(dplyr)
library(patchwork)

setwd("~/Dropbox (Personal)/conflict_analysis/")
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
direction = c("upper", "upper", "upper", "upper",
              "both", "upper", "upper", "upper",
              "upper", "upper", "lower", "both", "lower", 
              "lower")
keep = c("branch_outliers", "comphet_RCFV", "GC", "heterozygosity",
         "length", "max_PI", "when_max_PI", "mean_resid",
         "missing", "occupancy", "root_tip_var", "saturation_cval",
         "SH", "tree_len")
keep1 = c("branch outliers", "comp. het. RCFV", "GC", 
          "heterozygosity",
          "length", "max PI", "when max PI", "mean residuals",
          "missing", "occupancy", "root-tip variance", 
          "saturation C-value",
          "SH", "tree length")
names(keep1) = keep
names(direction) = keeper
d = d7[, c("locus", keeper)]
d1 = d %>% tidyr::gather("metric", "value", -locus)
d1 = split(d1, d1$metric)

get_outliers <- function(df) {
  iqr = quantile(df$value, c(0.25, 0.75), na.rm = T)
  iqrrange = diff(iqr)
  df$upper= iqr[2] + 1.5 * iqrrange
  df$lower = iqr[1] - 1.5 * iqrrange
  return(df)
}
do1 = lapply(d1, get_outliers)
for (i in 1:length(do1)) {
  metric = names(do1)[i]
  dir = direction[metric]
  if (dir == "upper") {
    do1[[i]] = do1[[i]] %>% mutate(outlier = ifelse(value >= upper, TRUE, FALSE))
  } else if  (dir == "lower") {
    do1[[i]] = do1[[i]] %>% mutate(outlier = ifelse(value <= lower, TRUE, FALSE))
  } else {
    do1[[i]] = do1[[i]] %>% mutate(outlier = ifelse((value > lower & value < upper), FALSE, TRUE))
  }
}
do2 = do.call("rbind", do1)
do2 %>% filter(outlier == TRUE) %>% group_by(metric) %>% summarize(n())
do3 = do2 %>% filter(outlier == TRUE)
length(unique(do3$locus))

t2 = read_csv(file = "family_lnl.csv")
t3 = full_join(t2, d)

t3 = t3[complete.cases(t3), ]
for (i in 1:length(keeper)) {
  t3[ ,paste(keeper[i], "_res", sep="")] = lm(pull(t3[, keeper[i]]) ~ t3$tree_len)$residuals
}
t3$diff = t3$gene_ll - t3$concat_ll
t3$diff_res = lm(t3$diff ~ t3$tree_len)$residuals

t3$type = "gene"
t3[grep("uce", t3$locus), "type"] = "uce"
t3[grep("AHE", t3$locus), "type"] = "AHE"

############
# how does fit of loci to species tree change across loci?
############
t3$log_diff = log(t3$diff)
t3$log_diff_res = log(t3$diff_res - min(t3$diff_res) + 1)
aa = aov(t3$log_diff_res ~ t3$type)
summary(aa)
alphas = c(0.5, 0.8, 0.1)
names(alphas) = c("AHE", "gene", "uce")
t3$alpha = alphas[t3$type]
fig1 = ggplot(t3, aes(type, log_diff)) + 
  geom_point(color = "gray", aes(alpha = alpha), 
             position = position_jitter(), 
             size = 0.7, shape = 16) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) + 
  ylim(quantile(t3$log_diff, c(0.01, 0.99))) + 
  ylab("log(gene tree LnL - \nconcat LnL)") + 
  theme(legend.position = "none")
save_plot(paste(figdir, "AHE_UCE_gene1.pdf", sep=""), fig1,
          base_height = 3.5, base_width = 4.5)

res = data.frame(matrix(NA, nrow = length(keeper), ncol = 3))
for (i in 1:length(keeper)) {
  cortest = cor.test(pull(t3[ ,paste(keeper[i], "_res", sep="")]), t3$diff_res, method = "spearman")
  res[i, 1] = keeper[i]
  res[i, 2] = cortest$estimate
  res[i, 3] = cortest$p.value
}

t4 = t3 %>% dplyr::select(locus, branch_outliers_res, 
                          comphet_RCFV_res, GC_res, 
            heterozygosity_res, max_PI_res, mean_resid_res, 
            missing_res, root_tip_var_res, saturation_cval_res, 
            tree_len_res, length_res,      
            when_max_PI_res, occupancy_res, SH_res) %>% 
  tidyr::gather("metric", "value", -locus)
t5 = full_join(t4, t3 %>% dplyr::select(locus, diff_res))

t6 = t5 %>% group_by(metric) %>% 
  mutate(min_q = quantile(value, 0.005), 
         max_q = quantile(value, 0.995)) %>% ungroup()
t7 = t6 %>% filter(value > min_q, value < max_q)

tt4 = t3 %>% dplyr::select(locus, comphet_RCFV,
                           SH, mean_resid, missing) %>% 
  tidyr::gather("metric", "value", -locus)
tt4$type = "gene"
tt4[grep("uce", tt4$locus), "type"] = "uce"
tt4[grep("AHE", tt4$locus), "type"] = "AHE"
tt6 = tt4 %>% group_by(metric) %>% 
  mutate(min_q = quantile(value, 0.005), 
         max_q = quantile(value, 0.995)) %>% ungroup()
tt7 = tt6 %>% filter(value > min_q, value < max_q)

# tt7$metric = keep1[tt7$metric]
tt8 = tt7 %>% select(-min_q, -max_q) %>% 
  tidyr::spread(metric, value)
alphas = c(0.3, 0.8, 0.03)
names(alphas) = c("AHE", "gene", "uce")
tt8$alpha = alphas[tt8$type]
vals = unique(tt7$metric)
plots = vector('list', length(vals))
for (i in 1:length(plots)) {
  plots[[i]] = ggplot(tt8, aes_string("type", vals[i])) + 
    geom_point(position = position_jitter(), 
               aes(alpha = alpha), 
               color = "gray", size = 0.7, shape = 16) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    ylab(keep1[vals[i]]) + xlab("") +
    theme(legend.position = "none",
          axis.title.x=element_blank())
    
}
fig2 = (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])
layout <- "
#B
AB
AB
#B
"
both = fig1 + fig2 + plot_layout(widths = c(1, 2), design = layout)
save_plot(paste(figdir, "AHE_UCE_gene.pdf", sep=""), both,
          base_height = 4, base_width = 10)
save_plot(paste(figdir, "AHE_UCE_gene.png", sep=""), both,
          base_height = 4, base_width = 10)

# fig2 = ggplot(tt7, aes(type, value)) + 
#   geom_point(position = position_jitter(), alpha = 0.05, color = "gray70") +
#   geom_boxplot(outlier.shape = NA, fill = NA) + 
#   facet_wrap(~metric, scales = "free") +
#   theme(strip.text = element_text(size=12)) +   
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(colour="white", fill="white"),
#         panel.border = element_rect(colour = "black"))
# save_plot(paste(figdir, "AHE_UCE_gene1.pdf", sep=""), fig1,
#            base_height = 2, base_width = 2.5)



t7 = t7 %>% filter(metric != 'tree_len_res')
t7$metric = keep1[gsub("_res$", "", t7$metric)]

t8 = split(t7, t7$metric)
lm_func <- function(x) {
  tlm = lm(x$diff_res ~ x$value)
  pval = summary(tlm)$coefficients[2, 4]
  rsq = summary(tlm)$adj.r.squared
  coef = as.numeric(coefficients(tlm))
  corr = cor.test(x$diff_res, x$value, method = "spearman")
  
  return(c(coef, pval, rsq, corr$p.value, corr$estimate))
}
res = lapply(t8, lm_func)
res = as.data.frame(do.call(rbind, res))
names(res) = c("intercept", "slope", "pval", "rsq", "corr_pval", "corr_rho")
res$metric = rownames(res)
res1 = res[res$corr_pval < 0.05, ]
  
xx = ggplot(t7, aes(value, diff_res)) + 
  geom_point(color = alpha("gray30", 0.07)) + 
  ylim(quantile(t7$diff_res, c(0.005, 0.995))) + 
  ylab("gene tree LnL - concat LnL") + xlab("") +
  geom_abline(data = res1, aes(slope = slope, intercept = intercept), color="red") + 
  facet_wrap(~metric, scales = "free", ncol = 5, nrow = 3, strip.position="bottom") +
  theme(strip.text = element_text(size=12)) +   
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"),
        strip.placement = "outside",
        panel.spacing = unit(1.5, "lines"))
save_plot(paste(figdir, "lnl_molevol.png", sep = ""), 
          xx, ncol = 5, nrow = 3,  base_height = 2.5, base_width = 3)

# plots = vector("list", length(t8))
# for (i in 1:length(t8)) {
#   x = t8[[i]]
#   x$title <- unique(x$metric)
#   plots[[i]] = ggplot(x, aes(value, diff_res)) + geom_point() + 
#     geom_smooth(span = 0.5) + xlab("") + ylab("") + facet_grid(. ~ title)
# }