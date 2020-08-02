library(phytools)

# from phytools blog
cladebox<-function(tree,node,color=NULL,...){
  if(is.null(color)) color<-make.transparent("yellow",0.2)
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  h<-max(nodeHeights(tree))
  parent<-tree$edge[which(tree$edge[,2]==node),1]
  x0<-max(c(obj$xx[node]+obj$xx[parent])/2,obj$xx[node]-0.05*h)
  x1<-obj$x.lim[2]
  dd<-getDescendants(tree,node)
  y0<-min(range(obj$yy[dd]))-0.5
  y1<-max(range(obj$yy[dd]))+0.5
  polygon(c(x0,x1,x1,x0),c(y0,y0,y1,y1),col=color,
          border=0)
}

# pp = read_csv("~/Desktop/edges80.txt")
pp = read_csv("~/Desktop/iguanids/wo_oplurid/edges80.txt")
x = read_csv("~/Dropbox/squamate_tree/samples/squamate_phylogenomics_v10.csv")
# t = read.tree("~/Dropbox/squamate_tree/data/family_phylogeny/ExaML_result.concat_ind0.01_loci0.05_all_n5446.no_out.tre")
t = read.tree("~/Desktop/iguanids/wo_oplurid/concat_ind0.01_loci0.05_all_n5446.no_out.tre")
t = read.tree(text = write.tree(ladderize(t)))

t1 = t
t1$tip.label = paste(pull(x[match(t$tip.label, x$sample), "family"]), pull(x[match(t$tip.label, x$sample), "species"]), sep="_")

outdir = '/Users/sonal/Desktop/iguanids/wo_oplurid/edges//'
for (i in 1:nrow(pp)) {
  nodename = pull(pp[i, "nodename"])
  sps = strsplit(pull(pp[i, "tips"]), ":")[[1]]
  node = findMRCA(t, sps)
  
  conf = read_csv(paste(outdir, nodename, ".csv", sep=""))
  conf =  conf %>% arrange(desc(percent))
  
  out = pdf(paste(outdir, nodename, '.pdf', sep=""), height=5, width=5)
  par(mar=c(0,0,3,0))
  if (nrow(conf) > 0) {
    for (j in 1:nrow(conf)) {
      tipcols = rep("black", Ntip(t))
      csps1 = strsplit(pull(conf[j, "species1"]), ":")[[1]]
      csps2 = strsplit(pull(conf[j, "species2"]), ":")[[1]]
      tipcols[t$tip.label %in% csps1] ="red"
      tipcols[t$tip.label %in% csps2] ="blue"
      plot(t1, tip.color=tipcols, cex=0.7)
      cladebox(t1, node, make.transparent("steelblue",0.2))
      title(paste("support:", round(pull(conf[j, "percent"]), 3)))
    }
  }  else {
      tipcols = rep("black", Ntip(t))
      plot(t1, tip.color=tipcols, cex=0.7)
      cladebox(t1, node, make.transparent("steelblue",0.2))
      title("no significant conflict")
  }
  dev.off()
}