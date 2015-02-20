source("~/quasr/suff_stats_quasr.r")
library(parallel)
thresh<-function(x,lam){
	# s=solve(A)%*%x
# 	if(lam==0) return(s)
# 	max(1-lam/crossprod(x,s)^(1/2),0)*s
	if(lam==0) return(x)
	max(1-lam/sum(x^2)^(1/2) , 0)*x
}

admmSel<-function(dat,lam,m1,m2,rho=1,cores=1){
	K = Kvec2(dat,m1,m2)
	G = Gamma.ret2(dat,m1,m2,cores=dim(dat)[2])
	d=dim(dat)[2]
	Ginv=mclapply(1:d,function(x) G[[x]]$v %*% ((1/(G[[x]]$d^2+rho)) * t(G[[x]]$v)),mc.cores=d)
	dd=dim(G[[1]]$v)[1]
	znode = matrix(0,d,m1)
	elist= get.edgelist(graph.full(d))
	zedge = array(0,dim=c(m2,m2,dim(elist)[1]))
	tol=Inf
	x=xold=y=zz=lapply(1:d,function(x)rep(0,dd))
	while(tol>1e-4){
	# x step
	xold=x
	x=mclapply(1:d,function(i){
		Ginv[[i]]%*%(-K[,i] - y[[i]] + rho*zz[[i]]) 
	},mc.cores=cores)
	#z-y step
	for(i in 1:d){
		for(j in i:d){
			if(i==j){	
				c= (i-1)*m2^2+1
				ind=c: (c+m1-1)
				znode[i,] = .5*(x[[i]][ind ] + y[[i]][ind]/rho )
				znode[i,] = thresh(znode[i,],lam/rho)
				y[[i]][ind] = y[[i]][ind] + rho*( x[[i]][ind] - znode[i,])
				zz[[i]][ind] = znode[i,]
			}else{
				e= which(elist[,1]%in%c(i,j) & elist[,2]%in%c(i,j))
				#c1 = ifelse(j>i , m1 + (j-2)*m2+1 ,(j-1)*m2+1 )
				#c2 = ifelse(j<i , m1 + (i-2)*m2+1 ,(i-1)*m2+1 )
				c1 = m1+(j-2)*m2^2+1
				c2 = (i-1)*m2^2+1
				ind1 =c1:(c1+m2^2-1)
				ind2 =c2:(c2+m2^2-1)
				zedge[,,e] = .5* (x[[i]][ind1] + x[[j]][ind2] +
								y[[i]][ind1]/rho + y[[j]][ind2]/rho) 
				zedge[,,e]=matrix( thresh(zedge[,,e],lam/rho) , m2 , m2)

				y[[i]][ind1] = y[[i]][ind1] + rho*(x[[i]][ind1] - c(zedge[,,e]))
				y[[j]][ind2] = y[[j]][ind2] + rho*(x[[j]][ind2] - c(zedge[,,e]))
				zz[[i]][ind1] = c(zedge[,,e])
				zz[[j]][ind2] = c(zedge[,,e])
			}
		}
	}
	#print(mean(abs(unlist(x)-unlist(zz))))
	tol = .5*mean(abs(unlist(xold)-unlist(x))) + .5*mean(abs(unlist(x)-unlist(zz)))
	print(paste("tol:",tol))
	}
	return(list(znode=znode,zedge=zedge))#list(z=zz,x=x,y=y,znode=znode,zedge=zedge))
}


library(huge)
d=5;n=250
dat=(huge.generator(n=n,d=d,graph="random",prob=.1)$data )/8 + .5
dat[dat<0]=.001
dat[dat>1]=.999
#dat=matrix(runif(d*n),n,d)
m1=3
m2=4

#Rprof(line.profiling=T)
admm=admmSel(dat,lam=.2,m1=m1,m2=m2,rho=1)
Rprof(NULL)
summaryRprof(lines="show")



################################


# coorDesc<-function(dat,lam,m1,m2){
# 	gr=Gamma.ret(dat,m1,m2)
# 	gn=Gamma.node(dat,m1,m2)
# 	kvn=Kvec.node(dat,m1)
# 	kve=Kvec.edge(dat,m2)
# 	tol=Inf
# 	d=dim(dat)[2]
# 	out.node.old=out.node=matrix(0,d,m1)
# 	out.edge.old=out.edge=array(0,dim=c(m2,m2,d*(d-1)/2))
# 	elist=get.edgelist(graph.full(d))
# 	while(tol>1e-3){
# 		for(i in 1:d){
# 			for(j in i:d){
# 				A = (gn[[i]][[j]]+gn[[j]][[i]])
# 				if(i!=j){
# 					ee= which(elist[,1]==i & elist[,2]==j)
# 					ind11 = which( (elist[,1]==j | elist[,2]==j) & elist[,1]!=i)
# 					ind12 = which( (elist[,1]==i | elist[,2]==i) & elist[,2]!=j)
# 					xx=crossprod(gr[[i]][[ee]] , out.node[i,])+
# 					    gr[[j]][[ee]] %*% out.node[j,] +
# 					    c(kve[,,ee])
# 				}else{
# 					ind11 = which( (elist[,1]==j | elist[,2]==j) )
# 					ind12 = which( (elist[,1]==i | elist[,2]==i) )
# 					xx = kvn[i,]					
# 				}
# 				#print(paste(i," " ,j))
# 				#print(elist[ind11,])
# 				#print(elist[ind12,])
# 				xx=xx+Reduce("+",lapply(1:length(ind11), function(x){
# 					 if(elist[ind11[x],1]==j){
# 					 	gr[[i]][[ind11[x]]] %*% c(out.edge[,,ind12[x]]) 
# 					 }else{
# 					 	t(gr[[i]][[ind11[x]]]) %*% c(out.edge[,,ind12[x]]) 
# 					 }
# 					 }))
# 				xx=xx+Reduce("+",lapply(1:length(ind11), function(x){
# 				 if(elist[ind12[x],1]==i){
# 					 gr[[j]][[ind12[x]]] %*% c(out.edge[,,ind11[x]]) 
# 				}else{
# 					t(gr[[j]][[ind12[x]]]) %*% c(out.edge[,,ind11[x]]) 
# 				}
# 				 
# 				 }))
# 				if(i==j){
# 					out.node[i,]=thresh(-xx,A,lam)
# 				}else{
# 					out.edge[,,ee] = matrix( thresh(-xx,A,lam) , m2,m2)
# 				}
# 			}
# 		}
# 		#print(out.edge)
# 		#print(out.node)
# 		#break
# 		tol = mean(abs(out.edge-out.edge.old))
# 		tol= tol + mean(abs(out.node.old-out.node))
# 		print(tol)
# 		out.node.old = out.node
# 		out.edge.old=out.edge
# 	}
# 	return(list(edge=out.edge,node=out.node))
# }