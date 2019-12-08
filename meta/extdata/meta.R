setwd("/home/gsc/houston/ABCG2Meta")

library(rmeta)
data<-read.table("rawdata.txt",as.is=F,head=T,sep="\t")
data

allel<-data.frame(ncase=2*(data$TT+data$GT+data$GG),ncon=2*(data$GGCon+data$GTCon+data$TTCon),evencase=2*data$TT+data$GT, evencon=2*data$TTCon+data$TTCon)
TTvsGG<-data.frame(ncase=(data$TT+data$GG),ncon=(data$GGCon+data$TTCon),evencase=data$TT, evencon=data$TTCon)
TTvsGT<-data.frame(ncase=(data$TT+data$GT),ncon=(data$TTCon+data$GTCon),evencase=data$TT, evencon=data$TTCon)
GTvsGG<-data.frame(ncase=(data$GT+data$GG),ncon=(data$GGCon+data$GTCon),evencase=data$GT, evencon=data$GTCon)

allelcum<-cummeta(ncase,ncon,evencase,evencon,data=allel,statistic="OR",method="meta.MH")
TTvsGGcum<-cummeta(ncase,ncon,evencase,evencon,data=TTvsGG,statistic="OR",method="meta.MH")
TTvsGTcum<-cummeta(ncase,ncon,evencase,evencon,data=TTvsGT,statistic="OR",method="meta.MH")
GTvsGGcum<-cummeta(ncase,ncon,evencase,evencon,data=GTvsGG,statistic="OR",method="meta.MH")

summary(allelcum)
summary(TTvsGGcum)
summary(TTvsGTcum)
summary(GTvsGGcum)



allelcum<-cummeta(ncase,ncon,evencase,evencon,data=allel,statistic="OR",method="meta.DSL")
TTvsGGcum<-cummeta(ncase,ncon,evencase,evencon,data=TTvsGG,statistic="OR",method="meta.DSL")
TTvsGTcum<-cummeta(ncase,ncon,evencase,evencon,data=TTvsGT,statistic="OR",method="meta.DSL")
GTvsGGcum<-cummeta(ncase,ncon,evencase,evencon,data=GTvsGG,statistic="OR",method="meta.DSL")

summary(allelcum)
summary(TTvsGGcum)
summary(TTvsGTcum)
summary(GTvsGGcum)





d <- cummeta.summaries(b$logs, b$selogs, names=b$names,method="random", logscale=TRUE)
plot(d,summary.conf=TRUE)
summary(d)