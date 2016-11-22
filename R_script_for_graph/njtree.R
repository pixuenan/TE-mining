data<- read.delim("M:/Thesis/data/tipTABLE/Rider_16gn.tsv", header=TRUE)
test <- as.matrix(data)
mobilitydata <- test[,6:22];
mobilitydata[mobilitydata[,1:16]=="confirmed"]<-0;
mobilitydata[mobilitydata[,1:16]=="unconfirmed"]<-1;
class(mobilitydata)<-'numeric'

test_matrix <- t(mobilitydata[,1:16])

f <- function(xx) nj(dist.gene(xx))
nj_tree <- f(test_matrix)
nj_tree_root <- root(nj_tree,1,r=T)
nj_tree_root_boot <- boot.phylo(nj_tree, FUN = f, test_matrix,rooted=T,B=100)
#20 simulations of the tree
for (i in c(1:20))
{print(boot.phylo(nj_tree, FUN = f, test_matrix,rooted=T,B=100))}
#make the graph
library(ape)
library(phangorn)
library(Rphylip)
pal=brewer.pal(10,'Set3')
png("njtree.png",res=400,height=1500,width=3600)
par(mar=c(4,4,3,1),omi=c(0.1,0.1,0.1,0.1),mgp=c(3,0.5,0),
    las=1,mex=1.5,cex.main=1.6,cex.lab=1.5,cex.axis=1.5)
plot(nj_tree_root,cex=1.2,label.offset=1,no.margin=TRUE,show.tip.label = T,use.edge.length= T,edge.width = 1)
add.scale.bar(x = 100.003, y=2.4)
nodelabels(nj_tree_root_boot,adj = c(-0.3,0),frame = "none")
dev.off()
