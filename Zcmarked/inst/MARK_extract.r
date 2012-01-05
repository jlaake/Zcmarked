library(Zcmarked)
library(RMark)
# extract data from Access database
lastcohort=2007
maxyear=2009
zcdata=extract.Zc(lastyear=maxyear,lastcohort=lastcohort)
zcdata=zcdata[zcdata$cohort%in%1993:2007,]
zcdata$ch=substr(zcdata$ch,7,23)
zcdata=zcdata[,c(1,3)]

# process data for cjs model with beginning time of 1987
zcdata.proc=process.data(zcdata,begin.time=1987,groups="sex")
# create design data 
zc.ddl=make.design.data(zcdata.proc)

################### p Design data ######################################################
# Add age bins/covariates for p
# age bins for p; [1],[2],[3-4],[5-7],[8-14],[15+]
zc.ddl$p$agep=cut(zc.ddl$p$Age,c(1,2,3,5,8,13,24),right=FALSE)
################### Phi Design data ######################################################
# Add age bins/covariates for Phi
# age bins 0,1,2,[3-4],[5-7],[8-10],[11-14],[15-22]
zc.ddl$Phi$AgeS=cut(zc.ddl$Phi$Age,c(0,1,2,3,5,8,11,15,23),right=FALSE)

#####################  Model specifications ###################################################
p.global=list(formula=~-1+time:sex:agep)
Phi.global=list(formula=~-1+time:sex:AgeS)

# fit global model
global.model=mark(zcdata.proc,zc.ddl,model.parameters=list(Phi=Phi.global,p=p.global),invisible=FALSE)

# export to MARK
export.MARK(zcdata.proc,"GlobalZcSurvival",global.model)

# import into MARK and run median c-hat procedure
# resulting c-hat - saved in chat.out
