set.seed(2141)
source("~/quasr/suff_stats_quasr.r")
source("~/trw/startup.r")

source("~/quasr/admm.r")

library(parallel)
library(magrittr)


library(huge)

d=20
n=100
#m1=2
#m2=1
m1=2
m2=1
#source("~/trw/startup.r")
setwd("~/quasr")
rep=10
#h=huge.generaor(n=n,d=d,graph="random",prob=.9)
h=huge.generator(n=2*rep*n,d=d,graph="random", prob=.2)
#h=huge.generator(n=2*rep*n,d=d,graph="hub",g=5)
h$data<-scale(h$data,center=T,scale=T)/8+0.5
h$data<- sign(h$data-.5)*(abs(h$data-.5))^.6
h$data<-scale(h$data,center=T,scale=T)/4+0.5

h$data[h$data>=1] = .95
h$data[h$data<=0] = .05

roc=rochuge=list()

hout=h$theta[row(h$theta)>=col(h$theta)+1]
spar=sum(hout)

lammax=1.5
lammax=max(abs(unlist(Kvec2(h$data,m1,m2))))*1.5
lamseq=exp(seq(from=log(.01*lammax),to=log(lammax),length.out=25))

out=mclapply(1:rep,function(r){
#for(r in 1:rep){

#print(paste("rep number ",r))
ncount=1
#dat=h$data[((r-1)*max(nlist)+1):(r*max(nlist)),]
dat=h$data[((r-1)*n+1):(r*n),]
test=h$data[rep*n + ((r-1)*n+1):(r*n) ,]
#test=h$data[(rep*max(nlist)+(r-1)*max(nlist)+1):(rep*max(nlist)+r*max(nlist)),]
ncount=ncount+n

#dat=matrix(runif(d*n),n,d)




g<-graph.full(d,directed=F)
edgelist<- get.edgelist(g)
#g=graph.adjacency(h$theta,mode="undirected")
#edgelist=get.edgelist(g)
# clust<- clusters(g)$membership
edge_mean <- edge.means(dat, edgelist, m2, cores)
node_mean <- node.means(dat , m1,cores)
edge_mean_test <- edge.means(test, edgelist, m2, cores)
node_mean_test <- node.means(test , m1,cores)
# #weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])
# edge_coef<- array(0,dim=c(m2,m2,dim(edgelist)[1]))
# node_coef<- matrix(0,d,m2)
# #lam=.1
# alpha_start=1
# sub=g
# tree_obj<-edge.Start(sub,burn.in=10)
# wt<-E(tree_obj)$weight
# lambda_seq<- lapply(1:dim(edgelist)[1],
#  	                    function(x){
# 	                    norm(edge_mean[,,x]-tcrossprod(node_mean[edgelist[x,1],],node_mean[edgelist[x,2],]),"F")
#    	                 })
# lambda_seq<-simplify2array(lambda_seq)
# seqlp<-round( seq(from=1,to=min(1e10,length(lambda_seq)),length.out=25) )
# lp<-(lambda_seq%>%sort(decreasing=T))[seqlp]
# #lp=seq(from=1e-1,to= 1.5,length.out=25)
# lp=lamseq
# lpout=lapply(lp,function(i){
# 	eind<- which(lambda_seq > i)
# 	elist_new<-edgelist[eind,,drop=F]
# 	sub<-subgraph.edges(g,eids=eind,delete.vertices=F)
# 	return(sub[])
# })


admm=admmSel(dat,lam=lamseq,m1=m1,m2=m2,rho=1,elist=edgelist)

Q = huge.npn(dat)
lshuge=seq(from=.001,to=1.5,length.out=100)
hugeout=huge(Q,lambda=lshuge,method="glasso")

# roc=lapply(lpout,function(x){
# 	a=1 - sum( ( x[row(x)>=col(x)+1]==0)  & (hout!=0) )/spar
# 	b=1- sum( (x[row(x)>=col(x)+1]!=0 ) & (hout==0 ) )/(d*(d-1)/2-spar)
# 	return(list(a,b))
# })

roc=lapply(admm,function(x){
	a=1 - sum( (apply(x$ze,3,function(i) 1*any(i!=0)) ==0 )  & (hout!=0) )/spar
	b=1- sum( (apply(x$ze,3,function(i) 1*any(i!=0)) !=0 ) & (hout==0 ) )/(d*(d-1)/2-spar)
	return(list(a,b))
})

rochuge=lapply(hugeout$path,function(x){
	a=1 - sum( (x[row(x)>=col(x)+1] ==0)  & (hout!=0) )/spar
	b=1- sum( (x[row(x)>=col(x)+1] !=0) & (hout==0 ) )/(d*(d-1)/2-spar)
	return(list(a,b))
})

roc=matrix(unlist(roc),nrow=2)
rochuge=matrix(unlist(rochuge),nrow=2)
#roc2=matrix(unlist(roc2),nrow=2)

return(list(roc,rochuge))
},mc.cores=rep)


rocplot1= Reduce("+",lapply(out,function(x)x[[1]]))/rep
rocplot2= Reduce("+",lapply(out,function(x)x[[2]]))/rep
#rocplot3= Reduce("+",lapply(out,function(x)x[[3]]))/rep



lapply(out,function(x){
	ind=which.max(apply(x[[1]],2,sum))
	ot= x[[1]][,ind]
	return(ot)
})




pdf("rocquasrERnongauss3.pdf")

plot(rocplot1[1,],rocplot1[2,],type="l",lty=6,xlim=c(0,1),ylim=c(0,1),col="red",
	xlab="TN%",ylab="TP%")
points(rocplot2[1,],rocplot2[2,],type="l",lty=1,col="blue")
#points(rocplot3[1,],rocplot3[2,],type="l",lty=2)

legend("bottomleft",c("quasr","skeptic"),lty=c(6,1),col=c("red","blue"))
dev.off()