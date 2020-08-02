library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
theme_set(theme_cowplot())

setwd("~/Dropbox (Personal)/conflict_analysis/")
figdir = "~/Dropbox/conflict_analysis/manuscript/figures_v2/"

clades= c("anniellidae", "anomalepidadae", "bolyeridae", 
          "cylindro_uropelt", "dibamidae", "diplodactylidae",
          "eublepharid", "gecko_skink", "gymnophthalmidae",
          "homalopsidae", "iguanids", "lanthanotidae",
          "rhineuridae", "skinks", "snakes", "typhlopidae",
          "xenosauridae")
cnames = c("Anniellinae", "Anomalepididae", "Bolyeriidae", 
           "Cylindrophiidae & Uropeltidae", "Dibamidae", "Diplodactylidae",
           "Eublepharidae", "Gekkota & Scincoidea", "Gymnophthalmidae",
           "Homalopsidae", "Iguania", "Lanthanotidae",
           "Rhineuridae", "Scincoidea", "Serpentes", "Typhlopidae",
           "Xenosauridae")
names(cnames) = clades
dd5 = read_csv("dLNL.csv")
dd5$clade2 = cnames[dd5$clade]

win = read_csv("conflict_analysis_19April19.csv")
win2 = win %>% group_by(clade) %>% top_n(1, dLNL) %>% 
  ungroup() %>% select(clade, win = topology)
dd6 = left_join(dd5, win2) %>% mutate(support = ifelse(topology == win, TRUE, FALSE)) %>% select(clade, locus, dLNL, support)
dd6$dLNL2 = ifelse(dd6$support == TRUE, 1, -1) * dd6$dLNL

corrs = data.frame(clade1 = combn(unique(dd6$clade), 2)[1, ],
                   clade2 = combn(unique(dd6$clade), 2)[2, ],
                   corr1  = NA, corr2 = NA, stringsAsFactors = F)
overlap = data.frame(clade1 = combn(unique(dd6$clade), 2)[1, ],
                   clade2 = combn(unique(dd6$clade), 2)[2, ],
                   overlap  = NA, average = NA, pval = NA,
                   stringsAsFactors = F)
for (i in 1:nrow(corrs)) {
  sub1 = dd6 %>% filter(clade %in% corrs[i, 1:2]) %>% 
    select(-support, -dLNL2) %>% 
    filter(dLNL != 0) %>% tidyr::spread(clade, dLNL)
  sub2 = dd6 %>% filter(clade %in% corrs[i, 1:2]) %>% 
    select(-support, -dLNL) %>% 
    filter(dLNL2 != 0) %>% tidyr::spread(clade, dLNL2)
  
  tt = dd6 %>% filter(clade %in% corrs[i, 1:2]) %>% 
    select(-dLNL, -dLNL2)
  dups = tt %>% group_by(clade, locus) %>% 
    summarize(count = n()) %>% ungroup() %>% 
    filter(count > 1) %>% select(locus)
  tt1 = tt %>% filter(!(locus %in% pull(dups))) %>% 
    tidyr::spread(clade, support) 
  tt1 = tt1 %>% filter(complete.cases(tt1))
  over = nrow(tt1[which(tt1[ ,2] == FALSE & tt1[ ,3] == FALSE), ])
  
  sim = rep(NA, 100)
  for (runme in 1:100) {
    v1 = sample(c(rep(TRUE, length(which(tt1[ ,2] == TRUE))),
                rep(FALSE, length(which(tt1[ ,2] == FALSE)))))
    v2 = sample(c(rep(TRUE, length(which(tt1[ ,3] == TRUE))),
                  rep(FALSE, length(which(tt1[ ,3] == FALSE)))))
    sim[runme] = length(which(v1 == FALSE & v2 == FALSE))
  }
  
  overlap[i, "overlap"] = over / nrow(tt1)
  overlap[i, "average"] = mean(sim) / nrow(tt1)
  overlap[i, "pval"] = sum(sim > over) / 100
  
  corrs[i, "corr1"] = cor.test(pull(sub1[ , 2]), pull(sub1[, 3]))$estimate
  corrs[i, "corr2"] = cor.test(pull(sub2[ , 2]), pull(sub2[, 3]))$estimate
}

a = ggplot(corrs, aes(corr1)) + geom_histogram(bins = 50) +
  xlab("correlation in DLNL")
#   xlab(expression("correlation in " * italic("D"[LNL])))
overlap$diff = overlap$overlap - overlap$average
b = ggplot(overlap, aes(diff)) + geom_histogram(bins = 50) +
  xlab("overlap in conflicting genes")
ab = plot_grid(a, b, labels = c("A", "B"), align = "v")
save_plot(paste0(figdir, "conflict_overlap.pdf"), ab, 
          base_height = 2.8, base_width = 3.75, ncol = 2)
