name <- c("SL_RT_F159_16gn","SL_RT_F17_16gn","Rider_16gn")

unique_TIPs = function(n){
data<- read.delim(sprintf("~/Downloads/tipTABLE/%s.tsv",n), header=TRUE);
data_matrix <- as.matrix(data);
mobilitydata <- data_matrix[,6:22];
mobilitydata[mobilitydata[,1:16]=="confirmed"]<-0;
mobilitydata[mobilitydata[,1:16]=="unconfirmed"]<-1;
mobilitydata[mobilitydata[,1:16]=="i"]<-1;
class(mobilitydata)<-'numeric';
unique_per_access <- NULL;
for(j in 1:16){
unique_per_access<-c(unique_per_access,nrow(subset(mobilitydata,mobilitydata[,j]==1 & mobilitydata[,17]==1)))
};
return(unique_per_access)
}
TIPs_list <- NULL
for (n in name){
unique_per_access <- unique_TIPs(n);
TIPs_list <- c(TIPs_list,unique_per_access)
}
TIPs_list_matrix <- matrix(TIPs_list,nrow=16,ncol=3);
#make the graph
library(RColorBrewer)
pal=brewer.pal(9,'Set3')
barplot(t(TIPs_list_matrix),ylim=c(0,250),beside=T,border="white",las=2,
        col=pal[4:6],names.arg=colnames(mobilitydata)[1:16],
        legend.text=c('SL_RT_F159','SL_RT_F17','Rider'),
        ylab="Number of TE locations \nunique for the accession",
        args.legend=list(x = "topleft",bty="n",cex=1.2,horiz=F))
png("3TEs_mobility.png",res=400,height=2000,width=3800)
par(mar=c(4,4,3,1),omi=c(0.1,0.1,0.1,0.1),mgp=c(2,0.5,0),
    las=1,mex=1.5,cex.main=1.6,cex.lab=1.5,cex.axis=1.5)
dev.off()
