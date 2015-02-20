
library(huge)
library(clime)
library(glasso)
library(Rcpp)
library(clime)
library(parallel)
sourceCpp("./sparseprec.cpp")
set.seed(129488374)
len=20
err<-matrix(0,len,3)
syst=matrix(0,len,3)
sparsity<-matrix(0,len,3)


n<-100
#d=50
dlist<- c(30,50,100,120) #,120,150)
reps=3
sequ=log(rev(seq(from=exp(.06),to=exp(0.4),length.out=len)))
err=syst=sparsity=array(0,dim=c(len,3,length(dlist),reps))
roc=array(0,dim=c(len,3,2,length(dlist),reps))



dnum=0
for(d in dlist){
print(paste("dim : ",d))
dnum=dnum+1
h<-huge.generator(2*n*reps,d,graph="random",.1)
#h<-huge.generator(2*n*reps,d,graph="scale-free")

sequ=log(rev(seq(from=exp(1e-1),to=exp(0.5),length.out=len)))
# if(d>n){
# sequ=log(rev(seq(from=exp(.2),to=exp(0.7),length.out=len)))
# }
x<-h$data
	for(r in 1:reps){
		print(paste("rep number " ,r))
		#x<- matrix(rexp(n*d),n,d)
		x1<-x[((r-1)*n+1) :(r*n),]
		x2<-x[(reps*n+((r-1)*n+1)) : (reps*n+(r*n)),]
		sigma=cor(x1)
		sigma2<-cor(x2)
		spar=h$sparsity*d*(d-1)/2
		k=0
		#for(i in sequ){
		out=mclapply(sequ,function(i){
			k=k+1
			syst=err=sparsity=numeric(3)
			roc=matrix(0,3,2)
			# quasr<-sparseprec(sigma,i)
			# gl<-glasso(sigma,i)$wi
			# clime<-clime(sigma,i)$Omega[[1]]
			syst[1]<-system.time(quasr<-sparseprec(sigma,i))[1]
			syst[2]<-system.time(gl<-glasso(sigma,i)$wi)[1]
			syst[3]<-system.time(clime<-clime(sigma,lambda=i*.8,sigma=T)$Omega[[1]])[1]
			err[1]<- sum(diag(quasr%*%sigma2))-log(det(quasr))
			err[2]<-sum(diag(gl%*%sigma2))-log(det(gl))
			err[3]<-sum(diag(clime%*%sigma2))-log(det(clime))
			#TN
			ind= row(quasr)<=(col(quasr)-1)
			roc[1,1]<-1-sum( (h$theta[ind]!=0) & (quasr[ind]==0) )/spar
			#TP
			roc[1,2]<- 1- sum( (quasr[ind]!=0) & (h$theta[ind]==0) )/(d*(d-1)/2-spar)
			roc[2,1]<-1-sum( (h$theta[ind]!=0) & (gl[ind]==0) )/spar
			roc[2,2]<- 1- sum( (gl[ind]!=0) & (h$theta[ind]==0) )/(d*(d-1)/2-spar)
			roc[3,1]<-1-sum( (h$theta[ind]!=0) & (abs(clime[ind])<1e-2) )/spar
			roc[3,2]<- 1- sum( (abs(clime[ind])>1e-2) & (h$theta[ind]==0) )/(d*(d-1)/2-spar)
			sparsity[1]<-(sum(quasr[ind]!=0))
			sparsity[2]<-(sum(gl[ind]!=0))
			sparsity[3]<-(sum(abs(clime[ind])>1e-2))
		return(list(syst=syst,err=err,roc=roc,sparsity=sparsity))
		},mc.cores=detectCores())
		k=0
		for(ll in 1:length(sequ)){
		
		k=k+1
		syst[k,1,dnum,r]<-out[[k]]$syst[1]
		syst[k,2,dnum,r]<-out[[k]]$syst[2]
		syst[k,3,dnum,r]<-out[[k]]$syst[3]

		err[k,1,dnum,r]<-out[[k]]$err[1]
		err[k,2,dnum,r]<-out[[k]]$err[2]
		err[k,3,dnum,r]<-out[[k]]$err[3]

		roc[k,1,1,dnum,r]<-out[[k]]$roc[1,1]
		roc[k,1,2,dnum,r]<-out[[k]]$roc[1,2]
		roc[k,2,1,dnum,r]<-out[[k]]$roc[2,1]
		roc[k,2,2,dnum,r]<-out[[k]]$roc[2,2]
		roc[k,3,1,dnum,r]<-out[[k]]$roc[3,1]
		roc[k,3,2,dnum,r]<-out[[k]]$roc[3,2]
		sparsity[k,1,dnum,r]<-out[[k]]$sparsity[1]
		sparsity[k,2,dnum,r]<-out[[k]]$sparsity[2]
		sparsity[k,3,dnum,r]<-out[[k]]$sparsity[3]

		}
			}
	
}




rocout<-array(0,dim=c(3,2,length(dlist),reps))
sparout<-array(0,dim=c(3,length(dlist),reps))
for(i in 1:3){
	for(j in 1:length(dlist)){
		for(k in 1:reps){
			rocout[i,,j,k]<- roc[apply(err,2:4,which.min)[i,j,k],i,,j,k]
			sparout[i,j,k]<- sparsity[apply(err,2:4,which.min)[i,j,k],i,j,k]
		}
	}
}



roctab<-apply(rocout,1:3,mean)


errtab=apply(apply(err,2:4,min),1:2,mean)
colnames(errtab)<-dlist
rownames(errtab)<-c("quasr","glasso","clime")
write.table(
	x=errtab,
	sep=",",
	file="Error.csv"
) 



dplot=3

pdf("roc.pdf")
plot(apply(roc,1:4,mean)[,1,,dplot],type="l",xlim=c(0,1),ylim=c(0,1),col="red")
points(apply(roc,1:4,mean)[,2,,dplot],type="l",col="blue");
points(apply(roc,1:4,mean)[,3,,dplot],type="l",col="purple");

points(t(rocout[1,,dplot,]),col="red")
points(t(rocout[2,,dplot,]),col="blue")
points(t(rocout[3,,dplot,]),col="purple")

legend("bottomleft",c("quasr","glasso"),pch=1,col=c("red","blue"))
dev.off()

pdf("cpu.pdf")
plot(sequ,apply(syst,1:3,mean)[,1,dplot],pch=16,log="y",ylim=range(apply(syst,1:3,mean)[,1:3,dplot]),
ylab="CPU time (seconds (log))",xlab="# of included edges")
points(sequ,apply(syst,1:3,mean)[,2,dplot],pch=2,col="red")
points(sequ,apply(syst,1:3,mean)[,3,dplot],pch=2,col="red")

#points(sparsity[spind3,3],syst[spind3,3],pch=6,col="blue")
#abline(v=h$sparsity*d*(d-1)/2)
legend("topleft",c("new","glasso","clime"),pch=c(16,2,6))
dev.off()


pdf("regpath.pdf")
plot(sparsity[,1,dplot,1],err[,1,dplot,1],pch=10,
	ylim=range(err[,1:3,dplot,1]),
	xlim=range(sparsity[,1:3,dplot,1]+1),
	col="red",log="x",cex=.5)
points(sparsity[,2,dplot,1],err[,2,dplot,1],pch=13,col="blue",cex=.5)
points(sparsity[,3,dplot,1],err[,3,dplot,1],pch=5,col="purple",cex=.5)
abline(v=h$sparsity*dlist[dplot]*(dlist[dplot]-1)/2)
dev.off()

save.image("quasr.rdata")

# plot(sequ,err[,1],pch=10,ylim=range(err))
# points(sequ,err[,2],pch=13)
# abline(v=r)
# legend("topright",c("quasr","glasso"),pch=c(10,13))

