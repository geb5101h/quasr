rm(list=ls())
if(Sys.getenv("USER")=="ebjan"){
  setwd("~/trw")
}else{
  setwd("~/Dropbox/trw")
}


library(fastICA)
library(mclust)
#m=3
m1=12
m2=12
source("startup.r")
source("trw_fit.r")
setwd("~")
source("~/Forest Density/treeLoader.R");
loadTree()
####################
set.seed(3498)
####################
len=500
ntest=300

limit=400
out_len=5

source("~/quasr/suff_stats_quasr.r")

dat<-matrix( scan(file="~/trw/MEG_art",what=double()) , ncol=122, byrow=F)[1:(len+ntest),]
#save(file="meg.rdata",list=c("dat"))
#source("meg.rdata")
dat<-scale(dat,center=TRUE,scale=TRUE)
dat[abs(dat)>6.5]=6.5
#dat[which(abs(dat)>6)]<- sign(dat[which(abs(dat)>6)])*6
#data2<-fastICA(data2[1:1000,],n.comp=122)$S
dat<- apply(dat,2,function(x)(x-min(x)+.01)/(max(x)-min(x)+.02))

sam<-sample(1:(len+ntest))
train<-dat[sam[1:len],]
test<-dat[sam[(len+1):(len+ntest)],]

d<-dim(train)[2]
####################
source("~/quasr/admm.r")

setwd("~/quasr")
d=122
fitsave=admmSel(train,seq(from=.1,to=1.3,length.out=10),2,2)

save.image()