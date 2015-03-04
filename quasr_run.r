set.seed(123441)
source("~/quasr/suff_stats_quasr.r")
source("~/quasr/admm.r")


library(parallel)



library(huge)

nlist= floor(exp(seq(from=log(40),to=log(5000),length.out=25)))
d=20
m1=4
m2=4
source("~/trw/startup.r")
setwd("~/quasr")
rep=5
#h=huge.generaor(n=n,d=d,graph="random",prob=.9)
h=huge.generator(n=2*rep*max(nlist),d=d,graph="scale-free")
h$data<-scale(h$data,center=T,scale=T)/8+0.5
#h$data<- pbeta(h$data,2,1)
#h$data<-h$data-colMeans(h$data)
#h$data<-qexp(pnorm(h$data),rate=12)
h$data<- sign(h$data-.5)*(abs(h$data-.5))^.6
h$data<-scale(h$data,center=T,scale=T)/4+0.5

h$data[h$data>=1] = .99
h$data[h$data<=0] = .001

rpt=rpq=matrix(0,rep,length(nlist))
out=list()

for(r in 1:rep){
print(paste("rep number ",r))
ncount=1
kk=0
for(n in nlist){
kk=kk+1
dat=h$data[((r-1)*max(nlist)+1):(r*max(nlist)),]
test=h$data[(rep*max(nlist)+(r-1)*max(nlist)+1):(rep*max(nlist)+r*max(nlist)),]
ncount=ncount+n

#dat=matrix(runif(d*n),n,d)

lammax=max(abs(unlist(Kvec2(dat,m1,m2))))
lamseq=exp(seq(from=log(.01),to=log(lammax),length.out=5))

#admm=admmSel(dat,lam=lamseq,m1=m1,m2=m2,rho=1,n=n)
# Rprof(line.profiling=T)
# admm=admmSel(dat,lam=lamseq,m1=m1,m2=m2,rho=1,n=n)
# Rprof(NULL)
# summaryRprof(lines="show")


# ind=apply(admm[[14]]$ze,3,function(x)any(x!=0))
# g=graph.full(d)
# quasr=subgraph.edges(g,E(g)[ind])[]
# spar= sum(h$theta!=0)/2
# 
# roc=array(0,dim=c(2,25,10))
# for(i in 1:10){
# admm=admmSel(dat[((i-1)*n+1):(i*n),],lam=lamseq,m1=m1,m2=m2,rho=1,n=n)
# 	roc[1,,i]<- unlist(lapply(1:25,function(i){
# 	indq=apply(admm[[i]]$ze,3,function(x)any(x!=0))
# 	quasr=subgraph.edges(g,E(g)[indq],delete.vertices=F)[]
# 	ind= row(quasr)<=(col(quasr)-1)
# 	1-sum( (h$theta[ind]!=0) & (quasr[ind]==0) )/spar
# 	}))
# 	roc[2,,i]<- unlist(lapply(1:25,function(i){
# 	indq=apply(admm[[i]]$ze,3,function(x)any(x!=0))
# 	quasr=subgraph.edges(g,E(g)[indq],delete.vertices=F)[]
# 	ind= row(quasr)<=(col(quasr)-1)
# 	1- sum( (quasr[ind]!=0) & (h$theta[ind]==0) )/(d*(d-1)/2-spar)
# 	}))
# }
# plot(apply(roc[1,,],1,mean),apply(roc[2,,],1,mean),xlim=c(0,1),ylim=c(0,1),
# 	type="l");dev.off()
# 
# unlist(lapply(admm,function(x) sum(apply(x$ze,3,function(y)any(y!=0)))))
# Rprof(NULL)
# summaryRprof(lines="show")
# 
# 
# indq=apply(admm[[10]]$ze,3,function(x)any(x!=0))
# 	quasr=subgraph.edges(g,E(g)[indq],delete.vertices=F)[]


g<-graph.full(d,directed=F)
edgelist<- get.edgelist(g)
g=graph.adjacency(h$theta,mode="undirected")
edgelist=get.edgelist(g)
lammax=max(abs(unlist(Kvec2(dat,m1,m2,edgelist))))
lam=lammax/2
clust<- clusters(g)$membership
edge_mean <- edge.means(dat[1:n,], edgelist, m2, cores)
node_mean <- node.means(dat[1:n,] , m1,cores)
edge_mean_test <- edge.means(test[1:n,], edgelist, m2, cores)
node_mean_test <- node.means(test[1:n,] , m1,cores)
#weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])
edge_coef<- array(0,dim=c(m2,m2,dim(edgelist)[1]))
node_coef<- matrix(0,d,m2)
#lam=.1
alpha_start=1
sub=g
tree_obj<-edge.Start(sub,burn.in=10)
wt<-E(tree_obj)$weight
#wt=rep(1,d*(d-1)/2)
# trwout=parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
#                    edgelist,wt,lambda=lam,alpha_start,
#                    sub,cores=cores)


#apply(trwout$edge_c,3,function(x) any(x!=0))
#apply(admm$ze,3,function(x)any(x!=0))



edgelist=get.edgelist(graph.adjacency(h$theta))
lamseq=exp(seq(from=log(.001),to=log(lammax),length.out=25))
admm=admmSel(dat[1:n,],lam=lamseq,m1=m1,m2=m2,rho=1,elist=edgelist)
#trwpath<-trw_path(edge_coef,node_coef,node_mean,edge_mean,
#		edgelist,g=g,out_len=5,limit=d-2)
		#weights=NULL,alpha_start=alpha_start,
		#	g,cores=cores,limit=limit,relax=F,edge.opt=T,out_len=out_len)




risk_path_quasr<-lapply(admm,function(x){
	den1=rwrapper2(x$ze, x$zn,
		edge_mean,node_mean,edgelist,
    	wt, legen.Vec, legen.array, lam, alpha_start, alg = 0L)
    
	-trw_lik(node_mean_test,edge_mean_test,x$zn,x$ze,den1$part_fn)
})
   

risk_path_trw<-lapply(lamseq,function(l){
	ptf=parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
                   edgelist,wt,lambda=l,alpha_start,
                   sub,cores=cores)
	return(list(lik=-trw_lik(node_mean_test,edge_mean_test,ptf$vert_c,ptf$edge_c,ptf$part_fn),trw=ptf))
})

rpt[r,kk]=min(unlist(lapply(risk_path_trw,function(x)x$lik)))
rpq[r,kk]=min(unlist(risk_path_quasr))

out[[kk]]=list(quasr=admm[[which.min(unlist(risk_path_quasr))]],
	trw=risk_path_trw[[which.min(unlist(lapply(risk_path_trw,function(x)x$lik)))]]
	)

}
}


save.image("scorematchpath.rdata")

# 
# 
# 
# ptf=parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
#                    edgelist,wt,lambda=lamseq[which.min(risk_path_trw)],alpha_start,
#                    sub,cores=cores)
#                    
# den1=rwrapper2(admm[[which.min(risk_path_quasr)]]$ze, admm[[which.min(risk_path_quasr)]]$zn,
# 		edge_mean_test,node_mean_test,edgelist,
#     	wt, legen.Vec, legen.array, lam, alpha_start, alg = 0L)
# 
pdf("riskpath.pdf")
plot(nlist,apply(rpt,2,mean,na.rm=T),col="red",type="l",
	ylim=range(apply(rpt,2,mean,na.rm=T),apply(rpq,2,mean,na.rm=T)[-3]),log="xy",
	xlab="Number of samples",
	ylab="Held-out neg log lik",lty=1)
points(nlist[-3],apply(rpq,2,mean,na.rm=T)[-3],type="l",col="blue",lty=6)
dev.off()
# pdf("bivden.pdf")
# for(i in 1:dim(edgelist)[1]){
# 	par(mfcol=c(1,2))
# 	contour(den1$edge_den[,,i],nlevels=15);
# 	points(dat,cex=.01)
# 	contour(ptf$edge_den[,,i],nlevels=15);
# 	points(dat,cex=.01)
# }
# dev.off()
pdf("bivdenrisk.pdf",width=10,height=7)
par(mfcol=c(1,2))
contour(out[[1]]$trw$trw$edge_den[,,10])
den1=rwrapper2(out[[1]]$quasr$ze, out[[1]]$quasr$zn,
 		edge_mean,node_mean,edgelist,
     	wt, legen.Vec, legen.array, lam, alpha_start, alg = 0L)
contour(den1$edge_den[,,10])

contour(out[[10]]$trw$trw$edge_den[,,10])
den1=rwrapper2(out[[10]]$quasr$ze, out[[10]]$quasr$zn,
 		edge_mean,node_mean,edgelist,
     	wt, legen.Vec, legen.array, lam, alpha_start, alg = 0L)
contour(den1$edge_den[,,10])


con