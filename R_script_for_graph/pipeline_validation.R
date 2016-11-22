library("ggplot2")
coverage_data<- read.delim("M:/Thesis/data/coverage/coverage.tsv",header=TRUE)

#function to creat the formula from linear regression results
lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}
p <- ggplot(data = coverage_data, aes(x = estimate.number, y = total)) + 
  labs(x = "Number of TEs detected by coverage estimation", y = "Number of TEs detected by the pipeline")+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point(size=5)+theme(axis.text=element_text(size=20, color="black"),axis.title=element_text(size=20,face="bold"))
p1 <- p + geom_text(aes(x = 2000, y = 3000, label = lm_eqn(lm(coverage_data$total 
    ~ coverage_data$estimate.number, log10ed))), parse = TRUE,size=7)
#make the graph
png("pipeline_validation.png",res=200,height=1500,width=3600)
par(mar=c(4,4,3,1),omi=c(0.1,0.1,0.1,0.1),mgp=c(3,0.5,0),
    las=1,mex=1.5,cex.main=1.6,cex.lab=1.5,cex.axis=1.5)
p1
dev.off()
