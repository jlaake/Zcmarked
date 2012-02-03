if(length(grep("RMark",.packages()))!=0)detach("package:RMark")
library(Zcmarked)
library(marked)
library(plyr)
library(RMark)
# define models function
crm_models=define.models(zcdata.proc,maxyear)
# Time bins for p to handle lack of sampling in 1988 & 1989 [1988-1990],1991,...,
dp1=list(p=list(time.bins=c(1987,1990,1991:maxyear)))
############################################################################################
# build set of p models with single Phi for parameter counts
model.list=as.matrix(expand.grid(p=paste("p",0:8,sep="."),Phi="Phi.00"))
# Each model is written out to a file spec_p.#.Phi.##.rda for p counts but not run
crm.wrapper(model.list,zcdata.proc,new.ddl,crm_models,design.parameters=dp1,method=c("nlminb"))
p.count=apply(as.matrix(model.list),1,function.wrapper,fx=fx.par.count,par="p")
# Get Phi counts from models that have already been run
model.list=as.matrix(expand.grid(p=paste("p",6,sep="."),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep=".")))
# fit all models -- each model is written out to a file spec_p.#.Phi.##.rda
Phi.count=apply(as.matrix(model.list),1,function.wrapper,fx=fx.par.count,par="Phi")
############################################################################################
# Fit all models -- commented out
# create list of all models 9*76=684 models 
# model.list=as.matrix(expand.grid(p=paste("p",0:8,sep="."),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep=".")))
# fit all models -- each model is written out to a file p.#.Phi.##.rda
# crm.wrapper(model.list,zcdata.proc,zc.ddl,crm_models,design.parameters=dp1,method=c("nlminb"))
############################################################################################
#    Fit models with -p.5,p.6,p.7 -- commented out
# create list of all models 3*76=228 models 
# model.list=as.matrix(expand.grid(p=paste("p",5:7,sep="."),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep=".")))
# fit models -- each model is written out to a file p.#.Phi.##.rda
# crm.wrapper(model.list,zcdata.proc,zc.ddl,crm_models,design.parameters=dp1,method=c("nlminb"))
############################################################################################
# get aic values of models with p.5,p.6,p.7 - chat = 1.11
mat=as.matrix(expand.grid(p=c("p.5","p.6","p.7"),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep=".")))
aic=apply(mat,1,function.wrapper,fx=fx.aic,chat=1.11)
aic_df=data.frame(mat)
aic_df$aic=aic
aic_df$delta_aic=aic_df$aic-min(aic_df$aic)


# compute v-c matrices for best models and save results
modelset=paste("p",best_models[["p Model"]],"Phi",best_models[["Phi Model"]],sep=".")
for(modelname in modelset)
{
	eval(parse(text=paste("load(file='",modelname,".rda')",sep="")))
	eval(parse(text=paste("model=",modelname,sep="")))
	vcv=cjs.hessian(model)
	model$vcv=vcv
	eval(parse(text=paste(modelname,"=model",sep="")))
	eval(parse(text=paste("save(",modelname,",file='",modelname,".rda')",sep="")))
}	



# Age and sex-specific predictions
# Set up data for survival predictions using average covariate values for weight, environment and
# disease
newdata=NULL
for(i in 1987:2007)
  newdata=rbind(newdata,as.character(zcdata$brand[zcdata$cohort==i&zcdata$sex=="F"][1]),
			as.character(zcdata$brand[zcdata$cohort==i&zcdata$sex=="M"][1]))
avgWtAnomaly.byyear=ldply(split(zcdata$weight,list(zcdata$sex,zcdata$cohort)),mean)
newPhi=zc.ddl$Phi[zc.ddl$Phi$brand%in%newdata,]
newPhi=newPhi[order(newPhi$cohort),]
newPhi$.id=paste(newPhi$sex,newPhi$cohort,sep=".")
newPhi$seq=1:nrow(newPhi)
newPhi=merge(newPhi,avgWtAnomaly.byyear)
newPhi$weight=newPhi$V1
newPhi=newPhi[order(newPhi$seq),]
newPhi$seq=NULL
newPhi$V1=NULL
avgEnv=sapply(envcovdf[envcovdf$time<=2008 & envcovdf$time>=1989,],mean)[-1]
nr=nrow(newPhi)
newPhi[,48:63]= matrix(rep(avgEnv,nr),byrow=T,nrow=nr)
newPhi$leptom=mean(leptostatus$leptom[leptostatus$time>=1989&leptostatus$time<=2008])/100
newPhi$DA=mean(dastatus$DA[dastatus$time>=1989&dastatus$time<=2008])/100
newPhi$repro=0
newPhi$repro[newPhi$cohort==1992&newPhi$Age>=3]=1
# Loop over models in top set and store in list
savePhi=newPhi

# Compute model weights
weight=best_models[["Delta QAIC"]]
weight=exp(-.5*weight)
weight=weight/sum(weight)
# Compute predictions for each model and their v-c matrix and
# store in lists
estimate.list=vector("list",length=length(modelset))
vcv.list=vector("list",length=length(modelset))
i=0
for(modelname in modelset)
{
	i=i+1
	eval(parse(text=paste("load(file='",modelname,".rda')",sep="")))
	eval(parse(text=paste("model=",modelname,sep="")))
	dml=marked:::create.dm(newPhi,model$model.parameters$Phi$formula)
	dml=dml[newPhi$Cohort>1,]
	newPhi=newPhi[newPhi$Cohort>1,]
	dml=dml[newPhi$Cohort%in%4:5|newPhi$Age<=1,]
	newPhi=newPhi[newPhi$Cohort%in%4:5|newPhi$Age<=1,]
	dml=dml[newPhi$Time>=newPhi$Cohort,]
	newPhi=newPhi[newPhi$Time>=newPhi$Cohort,]
	dml=dml[!newPhi$Age%in%c(4,6:7,9,11,13:14,16:25),]
	newPhi=newPhi[!newPhi$Age%in%c(4,6:7,9,11,13:14,16:25),]
	estimates=plogis(as.vector(dml%*%model$beta[1:ncol(dml)]))
	deriv=estimates*(1-estimates)*dml
	real.vcv=as.matrix(deriv%*%model$vcv[1:ncol(dml),1:ncol(dml)]%*%t(deriv))
    estimate.list[[i]]=estimates
	vcv.list[[i]]=real.vcv
	if(i<length(modelset)) newPhi=savePhi
}
# Call model average to get model averaged values
estimate.mat=do.call("rbind",estimate.list)
mavg.results=model.average(list(estimate=estimate.mat,weight=weight,vcv=vcv.list))
mavg.results$vcv[is.nan(mavg.results$vcv)]=0
# store results for age- and sex-specific means and std errors
newPhi$Phi=mavg.results$estimate
newPhi$Phi.se=mavg.results$se
real.vcv=mavg.results$vcv
# Use reml to compute means, se and process variance values
Fpup.index=which(newPhi$Age==0&newPhi$sex=="F")
Mpup.index=which(newPhi$Age==0&newPhi$sex=="M")
Fyearling.index=which(newPhi$Age==1&newPhi$sex=="F")
Myearling.index=which(newPhi$Age==1&newPhi$sex=="M")
num=length(Fpup.index)
Fpup=var.components.reml(newPhi$Phi[Fpup.index],matrix(1,nrow=num,ncol=1),matrix(real.vcv[as.matrix(expand.grid(Fpup.index,Fpup.index))],ncol=num,nrow=num))
Mpup=var.components.reml(newPhi$Phi[Mpup.index],matrix(1,nrow=num,ncol=1),matrix(real.vcv[as.matrix(expand.grid(Mpup.index,Mpup.index))],ncol=num,nrow=num))
Fyearling=var.components.reml(newPhi$Phi[Fyearling.index],matrix(1,nrow=num,ncol=1),matrix(real.vcv[as.matrix(expand.grid(Fyearling.index,Fyearling.index))],ncol=num,nrow=num))
Myearling=var.components.reml(newPhi$Phi[Myearling.index],matrix(1,nrow=num,ncol=1),matrix(real.vcv[as.matrix(expand.grid(Myearling.index,Myearling.index))],ncol=num,nrow=num))
# for 2+ the mean is a single value; the means for ages 0 and 1 are replaced with reml estimates
non.repro=with(newPhi[newPhi$cohort==1991,],tapply(Phi,list(sex,Age),mean))
non.repro[1,1]= Fpup$beta$Estimate
non.repro[2,1]= Mpup$beta$Estimate
non.repro[1,2]= Fyearling$beta$Estimate
non.repro[2,2]= Myearling$beta$Estimate
non.repro.se=with(newPhi[newPhi$cohort==1991,],tapply(Phi.se,list(sex,Age),mean))
non.repro.se[1,1]= Fpup$beta$SE
non.repro.se[2,1]= Mpup$beta$SE
non.repro.se[1,2]= Fyearling$beta$SE
non.repro.se[2,2]= Myearling$beta$SE
repro=with(newPhi[newPhi$cohort==1992,],tapply(Phi,list(sex,Age),mean))
repro[1:2,1:3]=NA
repro[2,4:5]=NA
repro.se=with(newPhi[newPhi$cohort==1992,],tapply(Phi.se,list(sex,Age),mean))
repro.se[1:2,1:3]=NA
repro.se[2,4:5]=NA


#Table
Female.table=rbind(paste(sprintf("%.2f",repro[1,]),"(",sprintf("%.3f",repro.se[1,]),")",sep=""),
           paste(sprintf("%.2f",non.repro[1,]),"(",sprintf("%.3f",non.repro.se[1,]),")",sep=""))
row.names(Female.table)=c("Repro","Non-repro")
colnames(Female.table)=colnames(repro)
Female.table[Female.table=="NA (NA)"]=""

Male.table=rbind(paste(sprintf("%.2f",repro[2,]),"(",sprintf("%.3f",repro.se[2,]),")",sep=""),
		paste(sprintf("%.2f",non.repro[2,]),"(",sprintf("%.3f",non.repro.se[2,]),")",sep=""))
row.names(Male.table)=c("Repro","Non-repro")
colnames(Male.table)=colnames(repro)
Male.table[Male.table=="NA (NA)"]=""


#Plotting
eps=.07
plot(1:9,non.repro[1,],type="b",ylim=c(.5,1),xaxt="n",ylab="Annual survival rate",xlab="Age (years)",xlim=c(1,9.7))
points((1:9)+eps,repro[1,],type="b",pch=1,lty=2)
points((1:9)+2*eps,non.repro[2,],type="b",pch=2)
points(c(1:9)+3*eps,repro[2,],type="b",pch=2,lty=2)
axis(1,1:9,labels=c("0","1","2","3-4","5-7","8-9","10-11","12-14","15+"))
legend(3,.7,legend=c("Female Non-reproductive","Female Reproductive", "Male Non-reproductive", "Male Reproductive"),pch=c(1,1,2,2),lty=c(1,2,1,2),bty="n")


best_models=aic_df[aic_df$delta_aic<6,] 
best_models$p=sapply(best_models$p,function(x)as.numeric(strsplit(as.character(x),"\\.")[[1]][2]))
best_models$Phi=sapply(best_models$Phi,function(x)as.numeric(strsplit(as.character(x),"\\.")[[1]][2]))
best_models$aic=NULL
best_models=best_models[order(best_models$delta_aic),]
best_models$npar=p.count[best_models$p]+Phi.count[best_models$Phi]
best_models=data.frame(best_models[,1:2],best_models$npar,sprintf("%.1f",best_models$delta_aic))
names(best_models)=c("p Model","Phi Model","Number of parameters","Delta QAIC")
rownames(best_models)=NULL


# split a vector of capture histories into a matrix when they are all numeric
chmat=t(sapply(strsplit(ch,""),function(x)as.numeric(x)))
# split a vector of capture histories into a matrix when they are alpha-numeric (as in MS designs)
chmat=do.call("rbind",strsplit(ch,""))
# reverse process to take matrix and create vector of ch strings
ch=apply(chmat,1,paste,collapse="")

xx=marked:::process.ch(zcdata$ch,all=TRUE)
ages=t(apply(xx$First,1,cumsum)-1)
aliveatage=(ages*xx$FtoL)[,6:23]
reprovalue=zcdata[,43:60]
reproage=matrix(0,nrow=23,ncol=2)
for (i in 1:18)
{
	xx=tapply(aliveatage[,i]>0,list(factor(aliveatage[,i],levels=0:22),reprovalue[,i]),sum)
	xx[is.na(xx)]=0
	reproage=reproage+xx
}
reproage=reproage/rowSums(reproage)
barplot(reproage,beside=T)





