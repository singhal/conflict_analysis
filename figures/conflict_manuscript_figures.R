library(ape)
library(corrplot)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(phangorn)
library(phytools)
library(dplyr)
library(readxl)

############
d = read.csv("~/Dropbox (Personal)/squamate_tree/samples/squamate_phylogenomics_v11.csv", 
             stringsAsFactors = F)
setwd("~/Dropbox (Personal)/conflict_analysis/")
fam = read.csv("final_species_list_12July19.csv", stringsAsFactors = F)
outs = c("allMis1", "chrPic1", "galGal5", "hg38", "taeGut2")
figdir = "manuscript/figures_v2//"
ours = c("rablab_brazil", "pyron_sqcl2", "pyron_sqcl2_2", 
         "rablab_sqcl2", "rablab_rapid")

root_tree <- function(tree, drop = TRUE) {
  outs = c("hg38", "galGal5", "taeGut2", 
           "allMis1", "chrPic1")
  tree = root(tree, outs[1], resolve.root = TRUE)
  if (drop) {
    tree = drop.tip(tree, outs)
  } 
  return(tree)
}

count_trees <- function(d1) {
  d1 = as.data.frame(d1)
  d1$total_trees = rep(NA, nrow(d1))
  
  for (i in 1:nrow(d1)) {
    tips = strsplit(as.character(d1[i, 'tips']), ":")[[1]]
    cts = unlist(lapply(trees, tips_in_tree, tips))
    d1[i, "total_trees"] = length(cts[cts > 1])
  }
  return(d1)
}

tips_in_tree <- function(tree, tips) {
  return(sum(tips %in% tree$tip.label))
}
###########

###########
# tip data
###########

get_vals <- function(dd, x, tips, name) {
  dd[ ,x] = as.numeric(dd[, x])
  par(mar=c(6, 2, 1, 2), tck=-0.02)
  xvals = pretty(range(dd[, x], na.rm = T), 2)
  plot.new()
  plot.window(xlim=range(xvals, na.rm = T), ylim=c(1, nrow(dd)))

  for (i in 1:nrow(dd)) {
    if (rownames(dd)[i] %in% no) {
      points(dd[i, x], i, pch=16, col="red", cex=1)
    } else {
      points(dd[i, x], i, pch=16, col="black", cex=1)
    }
  }
  
  axis(1, at=xvals, labels=NA, line=1)
  axis(1, at=xvals, labels=xvals, line=0.8, lwd = 0, cex.axis=1.5)
  mtext(name, side=1, line=3.7, cex=1.3)
}

t = read.tree("concat/ExaML_result.concat_ind0.01_loci0.05_all_n5343")
tr = root(t, "hg38", resolve.root = T)
tro = drop.tip(tr, outs)
trod = chronopl(tro, 0.01)
trod = read.tree(text = write.tree(ladderize(trod)))

q = read.csv("tip_data.csv", stringsAsFactors = F) 
rownames(q) = q$ind 
q = q[trod$tip.label, ]

# identify added data
no = trod$tip.label[!(d[match(trod$tip.label, d$sample), "type"] %in% ours)]
cols = rep("black", length(trod$edge))
cols[trod$edge[, 2] %in% match(no, trod$tip.label)] = "red"

fams = data.frame(tips = trod$tip.label,  
                  family = d[match(trod$tip.label, d$sample), "family"],
                  stringsAsFactors = F)
# nodes
nn = read_xlsx("questionable_nodes.xlsx")
nodenum = rep(NA, nrow(nn))
for (i in 1:nrow(nn)) {
  nn2 = pull(nn[i, "involved_families"])
  climb =  pull(nn[i, "up_down"]) + 1
  tips = fams[fams$family == nn2, "tips"]
  tips = tips[tips %in% trod$tip.label]
  if (length(tips) == 1) {
    anc = Ancestors(trod, which(trod$tip.label == tips))[climb]
  } else {
    mrca = findMRCA(trod, tips)
    anc = Ancestors(trod, mrca)[climb]
  }
  nodenum[i] = anc
}
clades = nn$short_name


# locus count, missing, coverage, heterozygosity
pdf(paste(figdir, "data_quality.pdf", sep=""), width=7, height=4, pointsize=8)
par(mfrow=c(1,5))
par(mar=c(6,0,1,0))
plot(trod, show.tip.label = F, edge.width=0.6, cex=1.2, edge.color = cols)
for (i in 1:length(nodenum)) {
  nodelabels("", nodenum[i], frame ="none", pch = 21, bg = "black")
}
get_vals(q, "num_total", rownames(q), "number of loci")
get_vals(q, "avg_length", rownames(q), "avg. locus length")
get_vals(q, "avg_cov", rownames(q), "coverage")
get_vals(q, "heterozygosity", rownames(q), "heterozygosity")
dev.off()

###########
# tip table
###########

# name, museum voucher, species, family, source, SRA, num loci, tot len
q = read.csv("tip_data.csv", stringsAsFactors = F) 
q$avg_length = round(q$avg_length, 0)
q$species = d[match(q$ind, d$sample), "species"]
q$family = d[match(q$ind, d$sample), "family"]
q$type = d[match(q$ind, d$sample), "type"]
q$voucher = rep(NA, nrow(q))
q$SRA = rep(NA, nrow(q))
q = q[,c("ind", "voucher", "species", "family", "type", "num_total", "avg_length", "SRA")]
q = q[order(q$family), ]
write.csv(q, paste(figdir, "S1_sample_data.csv", sep=""), row.names = F)

#########
# figure -- phylogeny
#########

t2 = read.tree("bootstrap/RAxML_bipartitions.boot.no_out.dated.tre")
t2 = read.tree(text = write.tree(ladderize(t2)))

c2 = read.tree("compare_astral_concat/astral_to_concat.concon.tre")[[1]]
c2 = minRotate(c2, setNames(1:Ntip(t2), t2$tip.label))

# (A) Gekkota, (B) Scincoidea, 
# (C) Lacertoidea,
# (E) Serpentes, 
# (F) Anguimorpha, and (G) Iguania.
branchcols = rep("black", length(t2$edge.length))
ctips = list(c("QMJ57120", "CHUNB56581"),
             c("SEW6684", "UMMZ_237560"),
             c("CHUNB56842", "USNM576222"),
             c("UMMZ_244201_ra_wait", "UMFS_10680"),
             c("UMFS_10293", "UMFS_11782"),
             c("MVZ_230099", "CHUNB62191"))
cols = c("#e41a1c", "#377eb8", "#4daf4a",
         "#984ea3", "#ff7f00", "#a65628")
for (i in 1:length(ctips)) {
  mnode = findMRCA(t2, tips = ctips[[i]], type = "node")
  desc = getDescendants(t2, mnode)
  branchcols[ which.edge(t2, desc) ] = cols[i]
}

t2$tip.label = paste(fam[match(t2$tip.label, fam$original_sample), "final_family"], 
                     fam[match(t2$tip.label, fam$original_sample), "final_species"])

pdf(paste(figdir, "concat.pdf", sep = ""), width = 5, height = 10)
par(xpd = T, mar = c(0, 1, 0, 12))

plot(t2, cex = 0.7, show.tip.label = F, edge.color = branchcols)
tiplabels(gsub("_", " ", t2$tip.label), 
          cex = 0.7, frame = "none", adj = 0, font = 3)
for (i in 1:length(c2$node.label)) {
  conflict = as.numeric(c2$node.label[i])
  if (is.na(conflict)) {
    conflict = 0
  }
  bs = as.numeric(t2$node.label[i])
  if (is.na(bs)) {
    bs = 100
  }
  
  if (conflict == 0 & bs > 95) {
    nodelabels("", i + Ntip(t2), frame = "none", pch = 21, 
             bg = "gray20", cex = 1)
  }  else if (conflict == 0  &  bs < 95) {
    nodelabels("", i + Ntip(t2), frame = "none", pch = 21, 
               bg = "gray", cex = 1)
  } else if (conflict == 1  &  bs < 95) {
    nodelabels("", i + Ntip(t2), frame = "none", pch = 21, 
               bg = "white", cex = 1)
  } 
}
add.scale.bar()
dev.off()

#########
# figure  -- astral vs. concat
#########


t1 = read.tree("astral/trees_0.01_0.05_10_100_5e-05_SH.all.tre")
t1 = root_tree(t1)
t1 = read.tree(text = write.tree(ladderize(t1)))
t2 = read.tree("bootstrap/RAxML_bipartitions.boot")
t2 = root_tree(t2)
t2 = chronopl(t2, 0.1)

c1 = read.tree("compare_astral_concat/concat_to_astral.concon.tre")[[1]]
c2 = read.tree("compare_astral_concat/astral_to_concat.concon.tre")[[1]]

t2 = minRotate(t2, setNames(1:Ntip(t1), t1$tip.label))
c1 = minRotate(c1, setNames(1:Ntip(t1), t1$tip.label))
c2 = minRotate(c2, setNames(1:Ntip(t2), t2$tip.label))


# (A) Gekkota, (B) Scincoidea, 
# (C) Lacertoidea, (D) Episquamata,
# (E) Serpentes, 
# (F) Anguimorpha, and (G) Iguania.
clades = c("A", "B", "C", "D", "E", "F", "G")
ctips = list(c("QMJ57120", "CHUNB56581"),
             c("SEW6684", "UMMZ_237560"),
             c("CHUNB56842", "USNM576222"),
             c("SEW6684", "UMFS_9997"),
             c("UMMZ_244201_ra_wait", "UMFS_10680"),
             c("UMFS_10293", "UMFS_11782"),
             c("MVZ_230099", "CHUNB62191"))
cnodes1 = rep(NA, length(clades))
for (i in 1:length(ctips)) {
  cnodes1[i] = findMRCA(t1, tips = ctips[[i]], type = "node")
}
cnodes2 = rep(NA, length(clades))
for (i in 1:length(ctips)) {
  cnodes2[i] = findMRCA(t2, tips = ctips[[i]], type = "node")
}

t1$tip.label = paste(fam[match(t1$tip.label, fam$original_sample), "final_family"], 
                     fam[match(t1$tip.label, fam$original_sample), "final_species"])
t2$tip.label = paste(fam[match(t2$tip.label, fam$original_sample), "final_family"], 
                     fam[match(t2$tip.label, fam$original_sample), "final_species"])

pdf(paste(figdir, "astral_concat_comparison.pdf", sep = ""), width = 10, height = 11)
par(xpd = T, mar = c(0, 1, 0, 12), mfrow=c(1, 2))

plot(t2, cex = 0.7, show.tip.label = F)
tiplabels(gsub("_", " ", t2$tip.label), 
          cex = 0.7, frame = "none", adj = 0, font = 3)
nodelabels("", frame = "none", pch = 21, 
           bg = ifelse(c2$node.label %in% c("", "1"), "black", "red"),
           cex = ifelse(c2$node.label %in% c("", "1"), 0.2, 1.2))
nodelabels(ifelse(as.numeric(t2$node.label) < 95, "*", ""), 
           frame = "none", col = "dodgerblue1", adj = c(1, 1), cex = 2)
for (i in 1:length(cnodes)) {
  nodelabels(clades[i], cnodes2[i], adj = c(1.4, -0.4), frame ="none")
}

par(mar = c(0, 12, 0, 0))
plot(t1, cex = 0.7, show.tip.label = F, direction = "leftwards")
tiplabels(gsub("_", " ", t1$tip.label), 
          cex = 0.7, frame = "none", adj = 1, font = 3)
nodelabels("", frame = "none", pch = 21, 
           bg = ifelse(c1$node.label %in% c("", "1"), "black", "red"),
           cex = ifelse(c1$node.label %in% c("", "1"), 0.2, 1.2))
nodelabels(ifelse(as.numeric(t1$node.label) < 0.95, "*", ""), 
           frame = "none", col = "dodgerblue1", adj = c(0, 1), cex = 2)
for (i in 1:length(cnodes)) {
  nodelabels(clades[i], cnodes1[i], adj = c(-0.7, -0.4), frame ="none")
}
dev.off()

#########
# astral phylogeny w conflict
#########

draw_pp <- function(t, pp) {
  vals = matrix(0, Nnode(t), 5)
  
  pp = pp %>% mutate(supportP = support / total_trees, 
                     max_conflictP = max_conflict / total_trees,
                     other_conflictP = (conflict - max_conflict) / total_trees,
                     non_infP = (total_trees - total_inf ) / total_trees,
                     total = total_trees)
  
  for (i in 1:nrow(pp)) {
    sps = strsplit(as.character(pp[i, 'tips']), ":")[[1]]
    # sps = sps[!( sps %in% c("MVZ_163062", "MVZ161183"))]
    sps = sps[sps %in% t$tip.label]
    if (length(sps) > 1) {
      node = findMRCA(t, sps) - Ntip(t)
      row = as.vector(t(pp[i, c("supportP", "max_conflictP", 
                                "other_conflictP", "non_infP", "total")]))
      vals[node, ] = row
    }
  }
  return(vals)
}

trees = read.tree("rooted_21June18.trees")
d1 = read_csv("astral_edges80.csv")
d1 = count_trees(d1)

t = read.tree("astral/trees_0.01_0.05_10_100_5e-05_SH.all.tre")
t = root_tree(t, drop = TRUE)
t = read.tree(text = write.tree(ladderize(t)))

vals1 = draw_pp(t, d1)
vals1[1, ] = c(0, 0, 0, 1, 0)

pdf(paste(figdir, "astral_conflict.pdf", sep=""), height=12, width=8)
par(lwd=0.5, mar=c(0, 2, 0, 12), xpd = T)
tips = paste(fam[match(t$tip.label, fam$original_sample), "final_family"], 
             fam[match(t$tip.label, fam$original_sample), "final_species"], sep="_")
plot(t, show.tip.label =  F)
cols = c("#4daf4a", "#e41a1c", "#377eb8", "gray")
tiplabels(gsub("_", " ", tips), cex = 0.7, frame = "none", adj = 0, font = 3)
nodelabels(pie = vals1[ , 1:4], cex = 0.7, piecol=cols)
legend(x=0.02, y = 12,
       legend=c("support", "highest conflict", "other conflict", "non-informative"), 
       fill = cols, bty="n", cex=1.2)
dev.off()

#########
# concat phylogeny w conflict
#########

trees = read.tree("rooted_21June18.trees")
d1 = read_csv("concat_edges80.csv")
d1 = count_trees(d1)

t = read.tree("concat/ExaML_result.concat_ind0.01_loci0.05_all_n5343")
t = root_tree(t, drop = TRUE)
t = read.tree(text = write.tree(ladderize(t)))
t = chronopl(t, 0.01)

fams = data.frame(tips = t$tip.label,  
                  family = d[match(t$tip.label, d$sample), "family"],
                  stringsAsFactors = F)
# pick one species per family
keep = pull(fams %>% group_by(family) %>% 
  filter(row_number(tips) < 2) %>% ungroup() %>% select(tips))
t2 = drop.tip(t, setdiff(t$tip.label, keep))

vals2 = draw_pp(t2, d1)
vals_full = draw_pp(t, d1)
vals2[1, ] = c(0, 0, 0, 1, 0)

# level of support for family
famnode = rep(FALSE, nrow(vals_full))
for (i in 1:nrow(vals_full)) {
  nn = i + Ntip(t)
  tt = t$tip.label[Descendants(t, nn, type = "tips")[[1]]]
  ff = unique(fams[fams$tips %in% tt, "family"])
  if (length(ff) == 1) {
    famnode[i] = TRUE
  }
}

# nodes
nn = read_xlsx("questionable_nodes.xlsx")
nodenum = rep(NA, nrow(nn))
for (i in 1:nrow(nn)) {
  nn2 = pull(nn[i, "involved_families"])
  climb =  pull(nn[i, "up_down"]) + 1
  tips = fams[fams$family == nn2, "tips"]
  tips = tips[tips %in% t2$tip.label]
  anc = Ancestors(t2, which(t2$tip.label == tips))[climb]
  nodenum[i] = anc
}
clades = nn$short_name

pdf(paste(figdir, "concat_conflict.pdf", sep=""), height=12, width=10)
par(lwd=0.5, mar=c(0, 2, 0, 6), xpd = T)
tips = d[match(t2$tip.label, d$sample), "family"]
plot(t2, show.tip.label =  F)
cols = c("#4daf4a", "#e41a1c", "#377eb8", "gray")
tiplabels(gsub("_", " ", tips), cex = 1, frame = "none", adj = 0, font = 3)
nodelabels(pie = vals2[ , 1:4], cex = 1, piecol=cols, lty = 0)
for (i in 1:length(nodenum)) {
  nodelabels(clades[i], nodenum[i], 
             adj = c(1.7, -0.4), frame ="none")
 }
legend(x=0.01, y = 8,
       legend=c("support", "highest conflict", "other conflict", "non-informative"), 
       fill = cols, bty="n", cex=1.2)
dev.off()

#######
# lnl comparison
#######

lnl = read_csv("family_lnl.csv")
lnl = lnl %>% mutate(lnl_diff = concat_ll - astral_ll)

t = read.tree("concat/concat1_to_concat2.concon.tre")[[1]]
t = read.tree(text = write.tree(ladderize(t)))
tips = paste(fam[match(t$tip.label, fam$original_sample), "final_family"], 
             fam[match(t$tip.label, fam$original_sample), "final_species"], sep="_")

pdf(paste(figdir, "concat_v1_v2.pdf", sep=""), width = 8, height = 8)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1, 2))

par(mar = c(14, 4, 12, 2))
hh = hist(lnl$lnl_diff, las = 1, main = "", 
     xlab = "", ylab = "", breaks = 30, col = "gray20", 
     border = "gray20", axes = F)
xvals =  pretty(lnl$lnl_diff, 3)[ pretty(lnl$lnl_diff, 3) < max(lnl$lnl_diff)]
axis(1, at =  pretty(lnl$lnl_diff, 3), labels = FALSE, 
     line = 0, tcl = -0.2)
mtext(xvals, side = 1, line = 0.2, 
      at = xvals)
mtext("concat ll - astral ll", side = 1, line = 1.2)
yvals = pretty(hh$counts, 3)[pretty(hh$counts, 3) < max(hh$counts)]
axis(2, at = yvals, 
     labels = FALSE, line = 0, tcl = -0.2)
mtext(yvals, side = 2, line = 0.2, las = 1, 
      at = yvals)
mtext("frequency", side = 2, line = 2.3)
# box(bty='L')

par(mar = c(0, 0, 0, 12), xpd = TRUE)
plot(t, cex = 0.7, show.tip.label = F)
tiplabels(gsub("_", " ", tips), 
          cex = 0.5, frame = "none", adj = 0, font = 3)
nodelabels("", frame = "none", pch = 21, 
           bg = ifelse(t$node.label %in% c("", "1"), "black", "red"),
           cex = ifelse(t$node.label %in% c("", "1"), 0.2, 0.7))
dev.off()

#################
# AHE vs. UCE
################

t1 = read.tree("UCE_AHE/AHE_to_UCE.astral.pp.tre.concon.tre")[[1]]
t1 = read.tree(text = write.tree(ladderize(t1)))
tips1 = paste(fam[match(t1$tip.label, fam$original_sample), "final_family"], 
              fam[match(t1$tip.label, fam$original_sample), "final_species"], sep="_")

t2 = read.tree("UCE_AHE/AHE_to_UCE.concat.pp.tre.concon.tre")[[1]]
t2 = minRotate(t2, setNames(1:Ntip(t1), t1$tip.label))
tips2 = paste(fam[match(t2$tip.label, fam$original_sample), "final_family"], 
              fam[match(t2$tip.label, fam$original_sample), "final_species"], sep="_")

pdf(paste(figdir, "AHE_vs_UCE.pdf", sep=""), width = 10, height = 8)
par(mar = c(0, 0, 0, 12), xpd = TRUE, mfrow = c(1, 2))

plot(t2, cex = 0.7, show.tip.label = F)
tiplabels(gsub("_", " ", tips2), 
          cex = 0.5, frame = "none", adj = 0, font = 3)
nodelabels("", frame = "none", pch = 21, 
           bg = ifelse(t2$node.label %in% c("", "1"), "black", "red"),
           cex = ifelse(t2$node.label %in% c("", "1"), 0.2, 0.7))
plot(t1, cex = 0.7, show.tip.label = F)
tiplabels(gsub("_", " ", tips1), 
          cex = 0.5, frame = "none", adj = 0, font = 3)
nodelabels("", frame = "none", pch = 21, 
           bg = ifelse(t1$node.label %in% c("", "1"), "black", "red"),
           cex = ifelse(t1$node.label %in% c("", "1"), 0.2, 0.7))

dev.off()

###############
# locus level plots
############### 

files = list.files(".", pattern = "loci_profiling", full.names = T)
d1 = read.csv(files[1], stringsAsFactors = F)
d2 = read.csv(files[2], stringsAsFactors = F)
d3 = read.csv(files[3], stringsAsFactors = F)
d4 = read.csv(files[4], stringsAsFactors = F)

d5 = full_join(d1, d2)
d6 = full_join(d5, d3)
d7 = full_join(d6, d4)
d7$comphet_RCFV = as.numeric(d7$comphet_RCFV)

keep = c("branch_outliers", "comphet_RCFV", "GC", "heterozygosity",
         "length", "max_PI", "when_max_PI", "mean_resid",
         "missing", "occupancy", "root_tip_var", "saturation_cval",
         "SH", "tree_len")
keep1 = c("branch outliers", "comp. het. RCFV", "GC", "heterozygosity",
          "length", "max PI", "when max PI", "mean residuals",
          "missing", "occupancy", "root-tip variance", "saturation C-value",
          "SH", "tree length")
d8 = d7 %>% select(keep)
names(d8) = keep1[ match(names(d8), keep) ]
d8$locus = seq(1, nrow(d8))
d9 = d8 %>% gather("metric", "value", -locus)

a = ggplot(d9, aes(value)) + geom_histogram(bins = 100, fill = "gray70") + 
  facet_wrap(~ metric, scales = 'free') + 
  theme(strip.text.x = element_text(size = 11)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))
save_plot(paste(figdir, "locus_metrics_histogram.pdf", sep=""),
          a, ncol = 3, nrow = 5, base_height = 1.5, base_width = 3)

###############
# locus level correlations
###############

M = cor(d8, use="pairwise.complete.obs")
colnames(M) = keep1
pdf(paste(figdir, "locus_metrics_correlation.pdf", sep = ""), width = 10, height = 10)
corrplot.mixed(M, tl.cex = 0.5, order = "hclust")
dev.off()

#######################
# dLNL histograms
#######################

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
dd6 = left_join(dd5, win2) %>% mutate(support = ifelse(topology == win, TRUE, FALSE))

# out = dd5 %>% group_by(clade) %>% summarise(outlier = quantile(dLNL, 0.99))
# out = inner_join(dd5, out)
# out = out %>% filter(dLNL < outlier)
short = nn$short_name
names(short) = nn$topology
dd6$clade3 = paste( short[dd6$clade], ": ", dd6$clade2, sep = "")

nn$clade3 = paste(nn$short_name, ": ", cnames[nn$topology], sep ="")
vars <- pull(nn %>% arrange(short_name) %>% select(clade3))

dd6$clade4 = factor(dd6$clade3, 
                    levels = vars, 
                    ordered = TRUE)
a = ggplot(dd6, aes(dLNL)) + 
  geom_histogram(bins = 100, aes(fill = support),
                 alpha=0.5, position="identity") +
  scale_fill_manual(values = c("#d7191c", "#2c7bb6")) +
  facet_wrap(~ clade4,scales = 'free_y') + 
  scale_x_log10() + geom_vline(xintercept=2, color = "black", linetype = "dotted") + 
  theme(strip.text.x = element_text(size = 12)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))
save_plot(paste(figdir, "dLNL.pdf", sep = ""), a, ncol = 5, nrow = 4,
          base_height = 2, base_width = 3)