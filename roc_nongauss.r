set.seed(2141)
source("~/quasr/suff_stats_quasr.r")
source("~/trw/startup.r")

source("~/quasr/admm.r")

library(parallel)



library(huge)

d=10
n=50
m1=3
m2=3
#source("~/trw/startup.r")
setwd("~/quasr")
rep=10
#h=huge.generaor(n=n,d=d,graph="random",prob=.9)
h=huge.generator(n=2*rep*n,d=d,graph="random",prob=2/d)
h$data<-scale(h$data,center=T,scale=T)/8+0.5
h$data<- sign(h$data-.5)*(abs(h$data-.5))^.6
h$data<-scale(h$data,center=T,scale=T)/4+0.5

h$data[h$data>=1] = .99
h$data[h$data<=0] = .001

roc=rochuge=list()

hout=h$theta[row(h$theta)>=col(h$theta)+1]
spar=sum(hout)

lammax=max(abs(unlist(Kvec2(h$data,m1,m2))))*1.5
lamseq=exp(seq(from=log(.01*lammax),to=log(lammax),length.out=25))

out=mclapply(1:rep,function(r){
#for(r in 1:rep){

#print(paste("rep number ",r))
ncount=1
#dat=h$data[((r-1)*max(nlist)+1):(r*max(nlist)),]
dat=h$data[((r-1)*n+1):(r*n),]
#test=h$data[(rep*max(nlist)+(r-1)*max(nlist)+1):(rep*max(nlist)+r*max(nlist)),]
ncount=ncount+n

#dat=matrix(runif(d*n),n,d)




g<-graph.full(d,directed=F)
edgelist<- get.edgelist(g)
#g=graph.adjacency(h$theta,mode="undirected")
#edgelist=get.edgelist(g)
# clust<- clusters(g)$membership
# #edge_mean <- edge.means(dat, edgelist, m2, cores)
# #node_mean <- node.means(dat , m1,cores)
# #edge_mean_test <- edge.means(test, edgelist, m2, cores)
# #node_mean_test <- node.means(test , m1,cores)
# #weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])
# edge_coef<- array(0,dim=c(m2,m2,dim(edgelist)[1]))
# node_coef<- matrix(0,d,m2)
# #lam=.1
# alpha_start=1
# sub=g
# tree_obj<-edge.Start(sub,burn.in=10)
# wt<-E(tree_obj)$weight



admm=admmSel(dat,lam=lamseq,m1=m1,m2=m2,rho=1,elist=edgelist)

Q = huge.npn(dat)
hugeout=huge(Q,lambda=lamseq,method="glasso")

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

return(list(roc,rochuge))
},mc.cores=20)

rocplot1= Reduce("+",lapply(out,function(x)x[[1]]))/rep
rocplot2= Reduce("+",lapply(out,function(x)x[[2]]))/rep








plot(rocplot1[1,],rocplot1[2,],type="l",lty=6,xlim=c(0,1),ylim=c(0,1))
points(rocplot2[1,],rocplot2[2,],type="l",lty=1)
dev.off()