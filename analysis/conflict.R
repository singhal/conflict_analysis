library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(ape)
library(phangorn)
library(Hmisc)

sort_vals <- function(x) {
  xsort = sort(x, decreasing = T)
  if (length(xsort) > 1) {
    return(xsort[1] - xsort[2])
  } else {
    return(0)
  }
}

read_res <- function(f) {
  f = read_csv(f, 
           col_names =  c("clade", "topology", "locus", "type", "value"),
           guess_max = 1000000)
  return(f)
}

t = read.tree("~/Dropbox (Personal)/conflict_analysis/rooted_21June18.out.trees")
## to identify tips with large amounts of missing data
# tt = lapply(t, function(tree) {return(tree$tip.label)})
# tt1 = table(unlist(tt))
# View(x[match(names(tt1[tt1 < 4000]), x$sample), c("sample", "species", "family")])
tnames = read.csv("rooted_21June18.csv")
names(t) = tnames$trees

# filter out loci that aren't informative b/c
# don't have at least two taxa in each group
# or if one of focal taxa isn't included

# key inds
key_inds = data.frame(inds = c("DWCT_R228", "LSUMZH9546", "CAS195955", "ISIS_393113"),
                      clade = c("bolyeridae", "dibamidae", "rhineuridae", "lanthanotidae"),
                      stringsAsFactors = F)

# do the filtering
edges = list.files("edges/", full.names = T)
evaluate_edge <- function(edgefile) {
  edge = read.tree(edgefile)
  name = gsub("^.*\\/", "", edgefile)
  name = gsub(".tre", "", name)
  name = strsplit(name, "\\.")
  
  clade1 = edge$tip.label[Descendants(edge, Ntip(edge) + 3, type = c("tips"))[[1]]]
  clade2 = setdiff(edge$tip.label, clade1)
  sum1 = sapply(t, function(tree) { return(sum(tree$tip.label %in% clade1))})
  sum2 = sapply(t, function(tree) { return(sum(tree$tip.label %in% clade2))})
  if (name[[1]][1] %in% key_inds$clade) {
    key_ind = key_inds[key_inds$clade == name[[1]][1], "inds"]
    sum3 = sapply(t, function(tree) { return(sum(tree$tip.label %in% key_ind))})
    res = sum1 > 1 & sum2 > 1 & sum3 > 0
  } else {
    res = sum1 > 1 & sum2 > 1
  }
  
  
  d = data.frame(clade = name[[1]][1],
                 topology = name[[1]][2],
                 locus = names(res),
                 complete = res,
                 stringsAsFactors = F)
  rownames(d) <- NULL
  return(d)
}
edgesres = lapply(edges, evaluate_edge)
edgesres = do.call(rbind, edgesres) 
ee = edgesres %>% group_by(clade, locus) %>% summarise(tot_complete = sum(complete), 
                                                       total = length(complete))
ee = ee %>% mutate(keep = tot_complete == total) 
ee = ee %>% select(clade, locus, keep)

x = read.csv("~/Dropbox/squamate_tree/samples/squamate_phylogenomics_v11.csv", 
             stringsAsFactors = F)

ff = list.files("results/", full.names = T)
d = lapply(ff, read_res)
d = do.call(rbind, d) 
dd = tidyr::spread(d, type, value)
dd$LNL = as.numeric(dd$LNL)

support = dd %>% group_by(clade, locus) %>% 
          summarise(num_support = (sum(CONFLICT == TRUE, na.rm = TRUE)))
dd2 = inner_join(dd, support)
dd2 = inner_join(dd2, ee)

# for each possible conflict
# for each topology
xx = dd2  %>% filter(num_support < 2, keep == TRUE) %>%
  group_by(clade, topology) %>%
  summarise(support = (sum(CONFLICT == TRUE, na.rm = TRUE)),
            conflict = (sum(CONFLICT == FALSE, na.rm = TRUE)))

# dlnl
dd3 = dd2 %>% group_by(clade, locus) %>% summarise(dLNL = sort_vals(LNL))
dd4 = dd2 %>% group_by(clade, locus) %>% filter(LNL == max(LNL))
dd5 = inner_join(dd3, dd4)
write_csv(dd5, "dLNL.csv")

dd6 = dd5 %>% filter(dLNL > 2) %>% group_by(clade, topology) %>% 
  summarise(dLNL_2 = sum(dLNL), num_genes_2 = n())
dd7 = dd5 %>% group_by(clade, topology) %>% 
  summarise(dLNL = sum(dLNL), num_genes = n())

out = dd5 %>% group_by(clade) %>% summarise(outlier = quantile(dLNL, 0.99))
out = inner_join(dd5, out)
dd10 = out %>% group_by(clade, topology) %>% filter(dLNL < outlier) %>%
  summarise(dLNL_out = sum(dLNL), num_genes_out = n())
dd8 = inner_join(dd7, dd6)
dd9 = inner_join(xx, dd8)
dd11 = inner_join(dd9, dd10)

ml = read_delim("~/Desktop/activeWork/conflict_analysis/conflicts_ML.csv", delim="\t")
dd12 = inner_join(ml, dd11)

write_csv(dd12 %>% arrange(clade, -dLNL), 
          "~/Desktop/activeWork/conflict_analysis/conflict_analysis_19April19.csv")


############  by type
dd2$locus_type = "gene"
dd2[grep("uce", dd2$locus), "locus_type"] = "uce"
dd2[grep("AHE", dd2$locus), "locus_type"] = "AHE"

# first look at how number of genes supporting changes
suptype = dd2 %>% filter(num_support < 2, keep == TRUE) %>% 
  group_by(clade, topology, locus_type) %>% 
  dplyr::summarise(support = sum(CONFLICT == TRUE, na.rm = T), conflict = sum(CONFLICT == FALSE, na.rm = T)) %>% 
  ungroup() %>% mutate(support_per = support / (support + conflict), 
                       conflict_per =conflict / (support + conflict) )

# do LNL analyses
dd5$locus_type = "gene"
dd5[grep("uce", dd5$locus), "locus_type"] = "uce"
dd5[grep("AHE", dd5$locus), "locus_type"] = "AHE"

dd6t = dd5 %>% filter(dLNL > 2) %>% group_by(clade, topology, locus_type) %>% 
  summarise(dLNL_2 = sum(dLNL), num_genes_2 = n())
dd7t = dd5 %>% group_by(clade, topology, locus_type) %>% 
  summarise(dLNL = sum(dLNL), num_genes = n())
dd8t = inner_join(dd7t, dd6t)
dd9t = dd8t %>% arrange(clade, locus_type, -dLNL) %>% filter(locus_type != "gene")
write.csv(dd9t, "conflict_analysis_by_type_12July19.csv")

