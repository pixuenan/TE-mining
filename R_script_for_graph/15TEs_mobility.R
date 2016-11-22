name <- c("SL_RT_F159_4gn_1000_4000","SL_RT_F164_4gn_4000_10000","SL_RT_F170_4gn_4000_10000",
          "SL_RT_F108_4gn_4000_10000","SL_RT_F120_4gn_4000_10000","SL_RT_F132_4gn_8000_14000",
          "SL_RT_F17_4gn_12000_16000","SL_RT_F272_4gn_4000_10000","SL_RT_F274_4gn_4000_10000",
          "SL_RT_F50_4gn_4000_10000","SL_RT_F94_4gn_14000_20000","TGRE1_4gn_4000_10000",
          "TARE1_4gn_4000_10000","Rider_4gn_4000_10000","Jinling_4gn_8000_14000")

for (n in name){
data<- read.delim(sprintf("M:/Thesis/data/tipTABLE/mobility/%s.tsv.ITAG",n), header=TRUE)
test <- as.matrix(data)
mobilitydata <- test[,6:10];
mobilitydata[mobilitydata[,1:4]=="confirmed"]<-0;
mobilitydata[mobilitydata[,1:4]=="unconfirmed"]<-1;

a<-NULL;
for(j in 1:4){for(i in c('1','2','3','4'))
{a<-c(a,nrow(subset(mobilitydata,mobilitydata[,j]=="1" & mobilitydata[,5]==i)))}};
d <- a;
for (i in c(1:length(a))){if (i%%4){d[i] <- a[i]+a[4]}}
b <- matrix(a,nrow=4,ncol=4);
if (length(strsplit(n,'_')[[1]])>4){
titlename <- paste(strsplit(n,'_')[[1]][1:3],collapse='_')
}else {
  titlename <- paste(strsplit(n,'_')[[1]][1])
}

max_num <- 40
if (titlename=='Rider'){
  max_num <-140
}
min_num <- b[4,1]
barplot(b[4:1,],ylim=c(min_num,min_num+max_num),main=titlename,
        col=pal,names.arg=colnames(mobilitydata)[1:4],xpd=FALSE
        )
legend("topleft",c('3','2','1'),bty="n",cex=1.5,horiz=TRUE,fill=pal[2:4])
}
#make the graph
library(RColorBrewer)
pal = brewer.pal(8,'Set2')[c(4,1,2,3)]
png("15TEs_mobility.png",res=200,height=1500,width=3600)
par(mar=c(4,4,3,1),omi=c(0.1,0.1,0.1,0.1),mgp=c(3,0.5,0),
    las=1,mex=1.5,cex.main=1.6,cex.lab=1.5,cex.axis=1.5)
par(mfrow=c(3,5))
dev.off()
