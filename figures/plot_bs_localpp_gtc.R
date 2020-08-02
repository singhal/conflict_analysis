library(ape)
library(phytools)
library(phangorn)
library(ggplot2)
library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(cowplot)
theme_set(theme_cowplot())
library(patchwork)

setwd("~/Dropbox/conflict_analysis/")
# examl tree
et = read.tree("bootstrap/RAxML_bipartitions.boot.no_out.tre")
etdf = data.frame(nodenum = seq(Ntip(et) + 1, Ntip(et) + Nnode(et)),
                  bootstrap = as.numeric(et$node.label))

# gt conflict
gt = read.csv("concat/concat_conflict_counts.csv", stringsAsFactors = F)
gt$nodenum = NA
for (i in 1:nrow(gt)) {
  tips = strsplit(gt[i, "tips"], ":")[[1]]
  gt$nodenum[i] = findMRCA(et, tips)
}
gt$per_gene = gt$support / gt$total_trees
etdf = left_join(etdf, gt %>% select(nodenum, per_gene))

# astral tree
at = read.tree("astral/trees_0.01_0.05_10_100_5e-05_SH.all.rooted.tre")
atdf = data.frame(nodenum = NA,
                  local_pp = as.numeric(at$node.label))
for (i in 1:Nnode(at)) {
  tips = at$tip.label[ Descendants(at, i + Ntip(at), type = "tips")[[1]] ]
  anc = findMRCA(et, tips)
  tips2 = et$tip.label[ Descendants(et, anc, type = "tips")[[1]] ]
  if (length(tips2) == length(tips) & sum(tips2 %in% tips) == length(tips)) {
    atdf$nodenum[i] = findMRCA(et, tips)
  }
}
etdf = left_join(etdf, atdf)


# dnl
nn = read_xlsx("questionable_nodes.xlsx")
dlnl = read_csv("conflict_analysis_19April19.csv")
dlnl2 = dlnl %>% group_by(clade) %>% 
  mutate(tot_dLNL = sum(dLNL)) %>% 
  ungroup() %>% mutate(frac_dLNL = dLNL / tot_dLNL) %>%
  group_by(clade) %>% top_n(1, frac_dLNL) %>% 
  select(clade, frac_dLNL) %>% ungroup()
dlnl3 = left_join(nn, dlnl2, by = c("topology" = "clade"))

# don't know if i really trust this
d = read.csv("~/Dropbox (Personal)/squamate_tree/samples/squamate_phylogenomics_v11.csv", 
             stringsAsFactors = F)
fams = data.frame(tips = et$tip.label,  
                  family = d[match(et$tip.label, d$sample), "family"],
                  stringsAsFactors = F)

nodenum = rep(NA, nrow(dlnl3))
for (i in 1:nrow(dlnl3)) {
  nn2 = pull(dlnl3[i, "involved_families"])
  climb =  pull(dlnl3[i, "up_down"]) + 1
  tips = fams[fams$family == nn2, "tips"]
  tips = tips[tips %in% et$tip.label]
  if (length(tips) == 1) {
    anc = Ancestors(et, which(et$tip.label == tips))[climb]
  } else {
    anc1 = findMRCA(et, tips)
    anc = c(anc1, Ancestors(et, anc1))[climb + 1]
  }
  nodenum[i] = anc
}
dlnl3$nodenum = nodenum
dlnl4 = dlnl3 %>% select(nodenum, frac_dLNL)
etdf = left_join(etdf, dlnl4)

a = ggplot(etdf, aes(bootstrap, local_pp)) + 
  geom_point() + 
  geom_jitter() +
  ylab("local pp")
b = ggplot(etdf, aes(bootstrap, per_gene)) + 
  geom_point() + 
  geom_jitter() +
  ylab("gene tree support")
c = ggplot(etdf, aes(bootstrap, frac_dLNL)) + 
  geom_point() + 
  geom_jitter() +
  ylab("dLNL support")
d = ggplot(etdf, aes(per_gene, frac_dLNL)) + 
  geom_point() + 
  geom_jitter() +
  xlab("gene tree support") +
  ylab("dLNL support")
abc = (a | b | c | d)
save_plot("~/Desktop/support_values.pdf", abc, ncol = 3, base_height = 3, base_width = 4)
