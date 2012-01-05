library(Zcmarked)
library(marked)
# define models function
crm_models=define.models(zcdata.proc,maxyear)
# Time bins for p to handle lack of sampling in 1988 & 1989 [1988-1990],1991,...,
dp1=list(p=list(time.bins=c(1987,1990,1991:maxyear)))
############################################################################################
# build set of p models with single Phi for parameter counts
model.list=as.matrix(expand.grid(p=paste("p",0:8,sep="."),Phi="Phi.00"))
# Each model is written out to a file spec_p.#.Phi.##.rda for p counts but not run
crm.wrapper(model.list,zcdata.proc,zc.ddl,crm_models,design.parameters=dp1,method=c("nlminb"))
p.count=apply(as.matrix(model.list),1,function.wrapper,fx=fx.par.count,par="p")
# Get Phi counts from models that have already been run
model.list=as.matrix(expand.grid(p=paste("p",6,sep="."),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep=".")))
# fit all models -- each model is written out to a file spec_p.#.Phi.##.rda
Phi.count=apply(as.matrix(model.list),1,function.wrapper,fx=fx.par.count,par="Phi")
############################################################################################
# Fit all models -- commented out
# create list of all models 9*76=684 models 
#model.list=expand.grid(p=paste("p",0:8,sep="."),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep="."))
# fit all models -- each model is written out to a file p.#.Phi.##.rda
#crm.wrapper(model.list,zcdata.proc,zc.ddl,crm_models,design.parameters=dp1,method=c("nlminb"))
############################################################################################
#    Fit models with -p.5,p.6,p.7 -- commented out
# create list of all models 3*76=228 models 
#model.list=expand.grid(p=paste("p",5:7,sep="."),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep="."))
# fit models -- each model is written out to a file p.#.Phi.##.rda
#crm.wrapper(model.list,zcdata.proc,zc.ddl,crm_models,design.parameters=dp1,method=c("nlminb"))
############################################################################################
# get aic values of models with p.5,p.6,p.7 - chat = 1.11
mat=as.matrix(expand.grid(p=c("p.5","p.6","p.7"),Phi=paste("Phi",formatC(0:75,digits=1,flag="00"),sep=".")))
aic=apply(mat,1,function.wrapper,fx=fx.aic,chat=1.11)
aic_df=data.frame(mat)
aic_df$aic=aic
aic_df$delta_aic=aic_df$aic-min(aic_df$aic)

