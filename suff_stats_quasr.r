library(orthopolynom)
library(igraph)
library(parallel)
legendre_Pl_array<-function(m,y){
  pf<-polynomial.functions(legendre.polynomials(m,normalized=T))
  lv<-sapply(X=1:(m+1),FUN=function(x){
    pf[[x]](y) #*sqrt(2*(x-1)+1)
  })
  return(t(lv)*sqrt(2))
}

len_deriv<-function(m,y){
	#y varies from [0,1]
	#returns x(1-x)deriv
	out=matrix(0,m+1,length(y))
	out[1,]=0
	lpa<-legendre_Pl_array(m,2*y-1)
	for(i in 2:(m+1)){
		k=i-1
		out[i,]<- -(k/2)*( (2*y-1)*lpa[i,]-sqrt((2*k+1)/(2*k-1))*lpa[i-1,] )
	}
	return(out)
} 	

len_deriv2<-function(m,y){
	#y varies from [0,1]
	#x(1-x)deriv2
	out=matrix(0,m+1,length(y))
	out[1,]=0
	out[2,]=0
	#lpa<-legendre_Pl_array(m,2*y-1)
	ld = len_deriv(m,y)
	for(i in 3:(m+1)){
		k=i-1
		#out[i,]<- sqrt(2*(i-1)+1)*(out[i-2,]/sqrt(2*(i-3)+1)+(1/2)*(2*(i-2)+1)*ld[i-1,]/sqrt(2*(i-2)+1) )
		out[i,] = out[i-2,]*sqrt((2*k+1)/(2*k-3)) + 2*(2*k-1)*sqrt(2*k+1)*ld[i-1,]/sqrt(2*k-1)
	}
	return(out)
}

Kvec.node<-function(dat,m){
	d=dim(dat)[2]
	out=matrix(0,d,m)
	out=lapply(1:d,function(i){
		ld1= t(-2*len_deriv(m,dat[,i])[-1,])*(2*dat[,i]-1)
		ld2= dat[,i]*(1-dat[,i])*t(len_deriv2(m,dat[,i])[-1,])
		rowMeans(t(ld1+ld2))
	})
	t(simplify2array(out))
}

Kvec.edge<-function(dat,m,edgelist=NULL){
	d=dim(dat)[2]
	g=graph.full(d)
	if(is.null(edgelist)) edgelist=get.edgelist(g)
	e=dim(edgelist)[1]
	if(d>=100){
		cores=detectCores()
	}else{
		cores=1
	}
	ret=mclapply(1:e,function(i){
		s=edgelist[i,1]
		t=edgelist[i,2]
		#out= crossprod(2*len_deriv(m,dat[,s])[-1,]*(2*dat[,s]-1) ,legendre_Pl_array(m,2*dat[,t]-1)[-1,] )
		#out = out + crossprod(legendre_Pl_array(m,2*dat[,s]-1)[-1,], 2*len_deriv(m,dat[,t])[-1,]*(2*dat[,t]-1) )
		#out= out+ crossprod(2*len_deriv2(m,dat[,s])[-1,]*(2*dat[,s]-1) ,legendre_Pl_array(m,2*dat[,t]-1)[-1,] )
		#out = out + crossprod(legendre_Pl_array(m,2*dat[,s]-1)[-1,], 2*len_deriv2(m,dat[,t])[-1,]*(2*dat[,t]-1) )
		out1= crossprod(-2*t(len_deriv(m,dat[,s])[-1,])*(2*dat[,s]-1) , t(legendre_Pl_array(m,2*dat[,t]-1)[-1,]) )
		out2 = legendre_Pl_array(m,2*dat[,s]-1)[-1,] %*% (-2*t(len_deriv(m,dat[,t])[-1,])*(2*dat[,t]-1) )
		out1= out1+ crossprod(dat[,s]*(1-dat[,s])*t(len_deriv2(m,dat[,s])[-1,]) ,t(legendre_Pl_array(m,2*dat[,t]-1)[-1,]) )
		out2 = out2 + legendre_Pl_array(m,2*dat[,s]-1)[-1,] %*% (dat[,t]*(1-dat[,t])*t(len_deriv2(m,dat[,t])[-1,]) )
		return(list(out1=out1/dim(dat)[1],out2=out2/dim(dat)[1]))
	},mc.cores=cores)
	return(list(array(unlist(lapply(ret,function(x)x$out1)),dim=c(m,m,e)),
		array(unlist(lapply(ret,function(x)x$out2)),dim=c(m,m,e))))
}

# Gamma.ret<-function(dat,m1,m2){
# 	d=dim(dat)[2]
# 	g=graph.full(d)
# 	ge=get.edgelist(g)	
# 	nedge=dim(ge)[1]
# 	out=mclapply(1:d,function(i){
# 		eout=mclapply(1:nedge,function(e){
# 			s=ge[e,1];t=ge[e,2]
# 			if(s==i){
# 				l1= len_deriv(m1,dat[,i])[-1,]
# 			}else{
# 				l1=numeric()
# 				lpa=legendre_Pl_array(m2,2*dat[,s]-1)[-1,]
# 				for(l in 1:m2){
# 					if(i<s){
# 						l1=rbind(
# 							l1,
# 							t(t(len_deriv(m2,dat[,i])[-1,]) * lpa[l,])
# 							#t(apply(len_deriv(m2,dat[,i])[-1,],1,function(x) x*lpa[l,,drop=F]))
# 						)
# 					}else{
# 						l1=rbind(
# 								l1,
# 								t(t(lpa) * len_deriv(m2,dat[,i])[-1,][l,])
# 								#t(apply(lpa,1,function(x) x*len_deriv(m2,dat[,i])[-1,][l,,drop=F]))
# 							)
# 					}
# 				}
# 			}
# 			if(t==i){
# 				l2= len_deriv(m1,dat[,i])[-1,]
# 			}else{
# 				l2=numeric()
# 				lpa=legendre_Pl_array(m2,2*dat[,t]-1)[-1,]
# 				for(l in 1:m2){
# 				if(i<t){
# 					l2=rbind(
# 						l2,
# 						t(t(len_deriv(m2,dat[,i])[-1,]) * lpa[l,])
# 
# 						#t(apply(len_deriv(m2,dat[,i])[-1,],1,function(x) x*lpa[l,,drop=F]))
# 						)
# 				}else{
# 					l2=rbind(
# 						l2,
# 						t(t(lpa) * len_deriv(m2,dat[,i])[-1,][l,])
# 						#t(apply(lpa,1,function(x) x*len_deriv(m2,dat[,i])[-1,][l,,drop=F]))
# 						)
# 				}
# 				}
# 			}
# 			tcrossprod(l1,l2)/dim(dat)[1]
# 		},mc.cores=detectCores())
# 		return(eout)
# 	},mc.cores=detectCores())
# 	return(out)
# 
# }

# Gamma.node<-function(dat,m1,m2){
# 	d=dim(dat)[2]
# 	out=mclapply(1:d,function(i){
# 		eout=mclapply(1:d,function(e){
# 			if(i==e){
# 				l1= len_deriv(m1,dat[,e])[-1,]
# 			}else{
# 				l1=numeric()
# 				lpa=legendre_Pl_array(m2,2*dat[,e]-1)[-1,]
# 				for(l in 1:m2){
# 				if(i<e){
# 					l1=rbind(
# 						l1,
# 						t(apply(len_deriv(m2,dat[,i])[-1,],1,function(x) x*lpa[l,,drop=F]))
# 					)
# 				}else{
# 					l1=rbind(
# 						l1,
# 						t(apply(lpa,1,function(x) x*len_deriv(m2,dat[,i])[-1,]))
# 					)
# 				}
# 				}
# 			}
# 			tcrossprod(l1,l1)/dim(dat)[1]
# 		},mc.cores=1)
# 		return(eout)
# 	},mc.cores=detectCores())
# 	return(out)
# }



Gamma.ret2<-function(dat,m1,m2,cores=1,ge=NULL){
	d=dim(dat)[2]
	if(is.null(ge))ge=get.edgelist(graph.full(d))	
	g=graph.edgelist(ge,directed=F)
	nedge=dim(ge)[1]
	out=mclapply(1:d,function(i){
		XX=numeric()
		for(j in sort(c(i,neighbors(g,i)))){	
				if(i==j){
					ld=len_deriv(m1,dat[,i])[-1,]
					XX=cbind(XX, t(ld) )
				}else if(i<j){
					e=which(ge[,1]%in%c(i,j) & ge[,2]%in%c(i,j))
					ld=len_deriv(m2,dat[,i])[-1,]
					lpa=legendre_Pl_array(m2,2*dat[,j]-1)[-1,]
					for(l in 1:m2){
						XX=cbind(XX, t(ld) * lpa[l,])
					}
				}else{
					e=which(ge[,1]%in%c(i,j) & ge[,2]%in%c(i,j))
					ld=len_deriv(m2,dat[,i])[-1,]
					lpa=legendre_Pl_array(m2,2*dat[,j]-1)[-1,]
					for(l in 1:m2){
						XX=cbind(XX, t(lpa) * ld[l,])
					}
				}
		}		
		#return(crossprod(XX,XX)/dim(dat)[1])
		#return(XX)
		return(svd(XX,nu=0))
 	},mc.cores=cores)
	return(out)

}

Kvec2<-function(dat,m1,m2,elist=NULL){
	if(is.null(elist)) elist=get.edgelist(graph.full(dim(dat)[2]))
	g=graph.edgelist(elist,directed=F)
	kvn=Kvec.node(dat,m1)
	kve=Kvec.edge(dat,m2,elist)
	#YY=numeric()
	YY=list()
	for(i in 1:d){
		XX=numeric()
		for(j in sort(c(i,neighbors(g,i)))){
			if(i==j){
				XX=c(XX,kvn[i,])
			}else if(i<j){
				e=which(elist[,1]%in% c(i,j) & elist[,2]%in%c(i,j))
				XX= c(XX, kve[[1]][,,e])
			}else{
				e=which(elist[,1]%in% c(i,j) & elist[,2]%in%c(i,j))
				XX= c(XX, kve[[2]][,,e])
			}
		}
		##YY=cbind(YY,XX)
		YY[[i]]=XX
	}
	return(YY)
}

# 
# Gamma.node(dat,3,3)
# d=5;n=40
# dat=matrix(runif(d*n),n,d)
# system.time()
# gr=Gamma.ret(dat,3,3)
# kvn=Kvec.node(dat,3)
# kve=Kvec.edge(dat,3)
# 
# 
# #Gamma: edges: d* d*(d-1)/2 m_2^2 x m_2^2 matrices
# #diag: d * d * matrices
# 
# Kvec.node(matrix(runif(100),10,10),5)
# 
# g=graph.full(10)
# Kvec.edge(matrix(runif(100),10,10),5,get.edgelist(g))
# 
# len_deriv(5,seq(from=0,to=1,length.out=50))
# 
# #for(i in 1:4) plot(legendre_Pl_array(10,2*seq(from=0,to=1,length.out=1000)-1)[i,])
# for(i in 1:20){
# plot(
# 	len_deriv2(25,seq(from=0,to=1,length.out=1000))[i,]
# )}
# #plot(
# #	len_deriv2(25,seq(from=0,to=1,length.out=1000))[5,]
# #)
# dev.off()
# 
