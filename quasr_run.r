set.seed(12341)
source("~/quasr/suff_stats_quasr.r")
source("~/quasr/admm.r")
library(parallel)



library(huge)
d=100;n=2500
h=huge.generator(n=n,d=d,graph="random",prob=.1)
dat=(h$data )/8 + .5
dat[dat<0]=.001
dat[dat>1]=.999
#dat=matrix(runif(d*n),n,d)
m1=4
m2=3
Rprof(line.profiling=T)
admm=admmSel(dat,lam=.1,m1=m1,m2=m2,rho=1)
Rprof(NULL)
summaryRprof(lines="show")

source("~/trw/startup.r")

g<-graph.full(d,directed=F)
edgelist<- get.edgelist(g)
clust<- clusters(g)$membership
edge_mean <- edge.means(dat, edgelist, m2, cores)
node_mean <- node.means(dat , m1,cores)

#weights<- rep(1, dim(edgelist)[1])  #rep(1,dim(edgelist)[1])
edge_coef<- array(0,dim=c(m2,m2,dim(edgelist)[1]))
node_coef<- matrix(0,d,m2)
lam=.1
alpha_start=1
sub=g
tree_obj<-edge.Start(sub,burn.in=10)
wt<-E(tree_obj)$weight
trwout=parallelTrwFit(edge_coef,node_coef,edge_mean,node_mean,
                   edgelist,wt,lambda=lam,alpha_start,
                   sub,cores=cores)
contour(trwout$edge_den[,,1]);dev.off()


apply(trwout$edge_c,3,function(x) any(x!=0))
apply(admm$ze,3,function(x)any(x!=0))


trwout$vert_c
admm$zn