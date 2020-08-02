library(phangorn)
library(ape)

t = read.tree("~/Desktop/activeWork/conflict_analysis/concat/ExaML_result.concat_ind0.01_loci0.05_all_n5446.no_out.tre")
tips = t$tip.label

x = read.csv("~/Dropbox/squamate_tree/samples/squamate_phylogenomics_v11.csv", stringsAsFactors = F)
x = data.frame(tips = tips, family = x[match(tips, x$sample), "family"], stringsAsFactors = F)
x1 = x[match(tips, x$sample), ]
internal = c("pyron_sqcl2", "pyron_sqcl2_2", "rablab_brazil", "rablab_rapid", "rablab_sqcl2")

d = read.csv("~/Desktop/activeWork/conflict_analysis/edges_to_test_v2.csv", stringsAsFactors = F)

u = unique(d$name)

for (i in 1:length(u)) {
  d1 = d[d$name == u[i], ]
  v = unique(d1$topology)
  for (j in 1:length(v)) {
  
    clade1 = x[x$family %in% d1[d1$topology == v[j], "family"], "tips"]
    clade2 = tips[!(tips %in% clade1)]
    clade2 = c(clade2, "outgroup")
  
    if (length(clade1) > 0 & length(clade2) > 0) {
      clade1 = paste(clade1, collapse=",")
      clade2 = paste(clade2, collapse=",")
      sink(paste('/Users/sonal/Desktop/edges/', u[i], '.', v[j], '.tre', sep=""))
      tree = paste('((', clade1, '),', clade2, ');', sep="")
      cat(tree)
      sink()
    } else {
      print(i)
    }
  }
}  