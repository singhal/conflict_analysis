t = read.tree("~/Dropbox/conflict_analysis/anomaly_zone/anomaly_zone/output.tre")
branchcols = rep("black", length(t$edge.length))
an = which(t$node.label != "\"0\"") + Ntip(t)
branchcols[ t$edge[, 2] %in% an ] = "red"
pdf("~/Dropbox/conflict_analysis/manuscript/figures_v2/anomaly_zone.pdf", 
    height = 10, width = 9)
par(mar = c(0, 0, 0, 14), xpd = TRUE)
plot.phylo(t, edge.color = branchcols, cex = 0.7, show.tip.label = F)
tiplabels(gsub("_", " ", t$tip.label), 
          cex = 0.7, frame = "none", adj = 0, font = 3)
dev.off()