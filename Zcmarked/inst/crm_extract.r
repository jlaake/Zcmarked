# detach RMark and attach marked and Zcmarked 
if(length(grep("RMark",.packages()))!=0)detach("package:RMark")
library(Zcmarked)
library(marked)
library(CalcurData)
library(CIPinnipedAnalysis)
# read in Domoic acid strandings (da.txt) and lepto stranding (lepto.txt)
dastatus=read.delim("DA.txt")
leptostatus=read.delim("lepto.txt")
# extract data from Access database
lastcohort=2007
maxyear=2009
zcdata=extract.Zc(lastyear=maxyear,lastcohort=lastcohort)
# modify weight field to subtract off sex-specific mean
avg.byyear=tapply(zcdata$weight,list(zcdata$sex,zcdata$cohort),mean)
for(cohort in 1987:lastcohort)
{
	zcdata$weight.yr[zcdata$sex=="F"&zcdata$cohort==cohort]=zcdata$weight[zcdata$sex=="F"&zcdata$cohort==cohort]-avg.byyear[1,cohort-1986]
	zcdata$weight.yr[zcdata$sex=="M"&zcdata$cohort==cohort]=zcdata$weight[zcdata$sex=="M"&zcdata$cohort==cohort]-avg.byyear[2,cohort-1986]
}	
zcdata$weight[zcdata$sex=="F"]=zcdata$weight[zcdata$sex=="F"]-mean(zcdata$weight[zcdata$sex=="F"])
zcdata$weight[zcdata$sex=="M"]=zcdata$weight[zcdata$sex=="M"]-mean(zcdata$weight[zcdata$sex=="M"])
# add missing time-dependent covariate values
zcdata$repro1987=0
zcdata$repro1988=0
zcdata$repro1989=0
zcdata$repro1990=0
zcdata$repro1991=0
zcdata$RepAge1987=0
zcdata$RepAge1988=0
zcdata$RepAge1989=0
zcdata$RepAge1990=0
zcdata$RepAge1991=0
zcdata$area11987=0
zcdata$area11988=0
zcdata$area11989=0
zcdata$area11990=0
zcdata$area21987=0
zcdata$area21988=0
zcdata$area21989=0
zcdata$area21990=0
zcdata$area21991=0
zcdata$td1987=0
# create an area3 time dependent covariate; area1=1 if seen in area1; area2=1 if seen in area2 but
# has never been seen in area1; area3=1 if seen in area3 but not in either area1 or area2 ever.
for (i in 1987:maxyear)
	zcdata[,paste("area3",i,sep="")]=1-zcdata[,paste("area1",i,sep="")]-zcdata[,paste("area2",i,sep="")]
# create a prepro time dependent covariate which is the same as the repro covariate except
# that it is shifted in time (repro1987=prepro1988)
for (i in 1987:(maxyear-1))
	zcdata[,paste("prepro",i+1,sep="")]=zcdata[,paste("repro",i,sep="")]
zcdata$prepro1987=0
# process data for cjs model with beginning time of 1987
zcdata$brand=NULL
zcdata.proc=process.data(zcdata,model="cjs",begin.time=1987,accumulate=T)
# create design data with time-varying fields area1,area2,area3,td,repro,prepro
# if a field is not identified as time varying then all of those fields with the year value
# will be included in data unless you specify fields argument.  Thus, even though prepro isn't
# used for Phi and repro isn't used for p, they are included here and then removed below.
zc.ddl=make.design.data(zcdata.proc,parameters=list(Phi=list(time.varying=c("area1","area2","area3","repro","td","prepro","RepAge")),
				p=list(time.varying=c("area1","area2","area3","repro","td","prepro","RepAge"))))
# remove unused fields in Phi and p design data								
zc.ddl$Phi$repro1987=NULL
zc.ddl$Phi$RepAge1987=NULL
zc.ddl$Phi$RepAge2009=NULL
zc.ddl$Phi$prepro2009=NULL
zc.ddl$Phi$area12010=NULL
zc.ddl$Phi$area22010=NULL
zc.ddl$Phi$area32010=NULL
zc.ddl$Phi$repro2010=NULL
zc.ddl$Phi$RepAge2010=NULL
zc.ddl$Phi$prepro2010=NULL
zc.ddl$Phi$prepro=NULL
zc.ddl$Phi$td2010=NULL
zc.ddl$Phi$area1=NULL
zc.ddl$Phi$area2=NULL
zc.ddl$Phi$area3=NULL
zc.ddl$Phi$td=NULL
zc.ddl$p$area11987=NULL
zc.ddl$p$area21987=NULL
zc.ddl$p$area31987=NULL
zc.ddl$p$td1987=NULL
zc.ddl$p$weight=NULL
zc.ddl$p$recap=NULL
zc.ddl$p$prepro1987=NULL
zc.ddl$p$repro1987=NULL
zc.ddl$Phi$recap=NULL
zc.ddl$p$TotalTimesResighted=NULL
zc.ddl$Phi$TotalTimesResighted=NULL
zc.ddl$p$initial.age=NULL
zc.ddl$Phi$initial.age=NULL
zc.ddl$p$repro=NULL
zc.ddl$p$RepAge=NULL
zc.ddl$Phi$area12009=NULL
zc.ddl$Phi$area22009=NULL
zc.ddl$Phi$area32009=NULL
zc.ddl$Phi$repro2009=NULL
zc.ddl$Phi$td2009=NULL
################### p Design data ######################################################
# Add age bins/covariates for p
# 2+ indicator covariate for p
zc.ddl$p$twoplus = 0
zc.ddl$p$twoplus[zc.ddl$p$Age>=2]= 1
# age bins for p; [1],[2],[3-4],[5-7],[8-11],[12+]
zc.ddl$p$agep=cut(zc.ddl$p$Age,c(1,2,3,5,8,12,24),right=FALSE)
# These age bins are used for <94 because the max age is 6 in 1993 for 1987 cohort
zc.ddl$p$agep94=cut(zc.ddl$p$Age,c(1,2,3,5,24),right=FALSE)
# indicator variables for <1994 (pre94) and >=1994 (p94plus)
zc.ddl$p$pre94=0
zc.ddl$p$pre94[zc.ddl$p$Time<6]=1
zc.ddl$p$post=0
zc.ddl$p$post[zc.ddl$p$Time>=2]=1
zc.ddl$p$p94plus=1-zc.ddl$p$pre94
#  Add a time covariate to p to limit area effects to 1991 and after
#zc.ddl$p$atime=0
#zc.ddl$p$atime[zc.ddl$p$Time>=3]=1
#  Add a male covariate to p to adjust area1/area2 affects relative to area3
zc.ddl$p$male=0
zc.ddl$p$male[zc.ddl$p$sex=="M"]=1
################### Phi Design data ######################################################
# Add age bins/covariates for Phi
# age bins 0,1,2,[3-4],[5-7],[8-9],[10-11],[12-14],[15-22]
zc.ddl$Phi$AgeS=cut(zc.ddl$Phi$Age,c(0,1,2,3,5,8,10,12,15,23),right=FALSE)
# age bins for 1+ [0-1],2,[3-4],[5-7],[8-10],[11-14],[15-22]
zc.ddl$Phi$age1p=cut(zc.ddl$Phi$Age,c(0,2,3,5,8,10,12,15,23),right=FALSE)
# age bins for 2+ [0-2],[3-4],[5-7],[8-10],[11-14],[15-22]
zc.ddl$Phi$age2p=cut(zc.ddl$Phi$Age,c(0,3,5,8,11,15,23),right=FALSE)
# indicator variables for f/m 1987 and 1988 cohorts which are 3 and 2 year survivals
zc.ddl$Phi$f1987=ifelse(zc.ddl$Phi$cohort==1987&zc.ddl$Phi$sex=="F"&zc.ddl$Phi$Time<=2,1,0)
zc.ddl$Phi$f1988=ifelse(zc.ddl$Phi$cohort==1988&zc.ddl$Phi$sex=="F"&zc.ddl$Phi$Time<=2,1,0)
zc.ddl$Phi$m1987=ifelse(zc.ddl$Phi$cohort==1987&zc.ddl$Phi$sex=="M"&zc.ddl$Phi$Time<=2,1,0)
zc.ddl$Phi$m1988=ifelse(zc.ddl$Phi$cohort==1988&zc.ddl$Phi$sex=="M"&zc.ddl$Phi$Time<=2,1,0)
# indicator variables for pup/non-pup (for pups exclude 1987 and 1988, 
# yearling (exclude 1987-1989), 2 year-old(exclude 1987-1989), and 2+(exclude 1987-1989) 
zc.ddl$Phi$pup = 0
zc.ddl$Phi$yearling = 0
zc.ddl$Phi$two = 0
zc.ddl$Phi$twoplus = 0
zc.ddl$Phi$three = 0
zc.ddl$Phi$threeplus = 0
zc.ddl$Phi$pup[zc.ddl$Phi$Cohort==zc.ddl$Phi$Time&zc.ddl$Phi$Time>=2] = 1
zc.ddl$Phi$yearling[zc.ddl$Phi$age==1&zc.ddl$Phi$Time>=3]= 1
zc.ddl$Phi$two[zc.ddl$Phi$age==2&zc.ddl$Phi$Time>=3]= 1
zc.ddl$Phi$twoplus[zc.ddl$Phi$Age>=2&zc.ddl$Phi$Time>=3]= 1
zc.ddl$Phi$three[zc.ddl$Phi$age==3&zc.ddl$Phi$Time>=3]= 1
zc.ddl$Phi$threeplus[zc.ddl$Phi$Age>=3&zc.ddl$Phi$Time>=3]= 1
# indicator variables for <1993 (pre93) and >=1993 (p93plus); >=1989 post (excludes
# 1989 if non-pup).
zc.ddl$Phi$p89plus=ifelse(zc.ddl$Phi$Time>=2&zc.ddl$Phi$time!="1999",1,0)
zc.ddl$Phi$p90plus=ifelse(zc.ddl$Phi$Time>=3&zc.ddl$Phi$time!="1999",1,0)
zc.ddl$Phi$p93plus=ifelse(zc.ddl$Phi$Time>=6,1,0)
zc.ddl$Phi$p94plus=ifelse(zc.ddl$Phi$Time>=7,1,0)
zc.ddl$Phi$post = 1
zc.ddl$Phi$post[zc.ddl$Phi$Time<2 | (zc.ddl$Phi$Time==2&zc.ddl$Phi$Age>=1) ]= 0
zc.ddl$Phi$i99=ifelse(zc.ddl$Phi$time==1999,0,1)
# create additional 0/1 variables for Phi; male, female, adult female and ENSO events
zc.ddl$Phi$male=0
zc.ddl$Phi$male[zc.ddl$Phi$sex=="M"]=1
zc.ddl$Phi$female= 1-zc.ddl$Phi$male
zc.ddl$Phi$adfemale= zc.ddl$Phi$female
zc.ddl$Phi$adfemale[zc.ddl$Phi$Age<5]=0
zc.ddl$Phi$admale= zc.ddl$Phi$male
zc.ddl$Phi$admale[zc.ddl$Phi$Age<5]=0
zc.ddl$Phi$EN92=0
zc.ddl$Phi$EN97=0
zc.ddl$Phi$EN03=0
zc.ddl$Phi$EN97[zc.ddl$Phi$time==1997]=1
zc.ddl$Phi$EN92[zc.ddl$Phi$time==1992]=1
zc.ddl$Phi$EN03[zc.ddl$Phi$time==2003]=1
zc.ddl$Phi$EN= zc.ddl$Phi$EN92+zc.ddl$Phi$EN97
# Add time intervals for initial pup cohort; these are setup in cjs/js models
zc.ddl$Phi$time.interval=1
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1987"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.816
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1988"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.814
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1989"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.729
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1990"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.701
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1991"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.732
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1992"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.671
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1993"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.737
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1994"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.784
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1995"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.792
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1996"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.797
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1997"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.808
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1998"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.726
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="1999"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.797
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2000"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.792
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2001"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.797
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2002"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.795
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2003"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.797
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2004"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.759
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2005"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.745
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2006"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.814
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2007"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.814
zc.ddl$Phi$time.interval[zc.ddl$Phi$cohort=="2008"&zc.ddl$Phi$Cohort==zc.ddl$Phi$Time]=0.814
zc.ddl$Phi$timing=ifelse(zc.ddl$Phi$pup==0,0,zc.ddl$Phi$time.interval-0.671)
#   Add lepto and DA covariates	
zc.ddl$Phi=merge_design.covariates(zc.ddl$Phi,dastatus)
zc.ddl$Phi=merge_design.covariates(zc.ddl$Phi,leptostatus)
zc.ddl$Phi$leptom= zc.ddl$Phi$leptom/100
zc.ddl$Phi$leptof= zc.ddl$Phi$leptof/100
zc.ddl$Phi$DA= zc.ddl$Phi$DA/100
#  Add in environmental covariates
#envcovdf=EnvironCovariates()
zc.ddl$Phi=merge_design.covariates(zc.ddl$Phi,envcovdf,bytime=TRUE)
zc.ddl$Phi=droplevels(zc.ddl$Phi)
zc.ddl$p=droplevels(zc.ddl$p)

dp1=list(p=list(time.bins=c(1987,1990,1991:maxyear)))
fixed.0=fix.parameters(zcdata.proc,occasions=2:3,value=c(0))
p.6=list(formula=~pre94:(sex*agep94)  + pre94:time + p94plus:time:sex:agep + p94plus:time:prepro + p94plus:male:prepro + p94plus:sex:agep:area1 + p94plus:sex:agep:area2, fixed=fixed.0, remove.intercept=TRUE)
Phi.44=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
mark6.44b=crm(zcdata.proc,zc.ddl,design.parameters=dp1,model.parameters=list(p=p.6,Phi=Phi.44),accumulate=T,method=c("BFGS"))


# save files for passing to Linux machine
save(zcdata,file="zcdata.rda")
save(envcovdf,file="envcovdf.rda")
save(zc.ddl,file="zc.ddl.rda")
save(zcdata.proc,file="zcdata.proc.rda")

xx=.Fortran("cjs", as.double(model_data$imat$chmat), as.double(Phibeta), 
		as.double(pbeta), as.double(model_data$imat$first), as.double(model_data$imat$last), 
		as.double(model_data$imat$freq), as.integer(model_data$imat$loc), 
		as.double(model_data$Phi.fixed), as.double(model_data$p.fixed), 
		as.double(model_data$time.intervals), as.integer(nrow(model_data$imat$chmat)), 
		as.integer(ncol(model_data$imat$chmat)), as.integer(nrow(model_data$Phi.fixed)), 
		as.integer(nrow(model_data$p.fixed)), lnl = double(1), p0 = double(nrow(model_data$imat$chmat)), 
		PACKAGE = "marked")

detach("package:RMark")
library(marked)
data(dipper)
dp=process.data(dipper,model="cjs")
fixed.0=fix.parameters(dp,occasions=3,value=0)
crm(dp,model.parameters=list(p=list(formula=~1,fixed=fixed.0)))
 