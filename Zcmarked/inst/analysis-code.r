if(length(grep("RMark",.packages()))!=0)detach("package:RMark")
library(Zcmarked)
library(marked)
library(plyr)
library(RMark)
# define models function
#' Define model set used in Zc survival analysis
#' 
#' Creates and returns a function that defines p and Phi model specifications for 
#' survival analysis.
#' 
#' @param zcdata.proc processed dataframe for Zc analysis
#' @param maxyear Last year used in analysis
#' @return function that can be used to select a subset of models to be
#' used in the analysis
#' 
#' @export define.models
#' 
#' @author Jeff Laake
define.models=function(zcdata.proc,maxyear)
{
	########################## Fixed parameters ##############################################
# Fixed p=0 and Phi=1 parameters to handle lack of sampling in 1988 & 1989
	fixed.0=fix.parameters(zcdata.proc,occasions=2:3,value=c(0))
	fixed.1=fix.parameters(zcdata.proc,occasions=1:2,value=c(1))
	#####################  Model specifications ###################################################
#    p
# constant 
	p.0=list(formula=~1)
# additive time + sex + age class
	p.1=list(formula=~time+sex+agep,fixed=fixed.0)
# for time <1994, sex:agep and for >= 1994, time:sex:agep 
	p.2=list(formula=~pre94:(sex*agep94) + pre94:time + p94plus:time:sex:agep, fixed=fixed.0,remove.intercept=TRUE)
# for time <1994, sex:agep and for >= 1994, time:sex:agep + sex:reproductive 
	p.3=list(formula=~pre94:(sex*agep94) + pre94:time + p94plus:time:sex:agep + p94plus:prepro + p94plus:male:prepro, fixed=fixed.0,remove.intercept=TRUE)
# for time <1994, sex:agep and for >= 1994, time:sex:agep + sex:reproductive + sex*area (>1994)
	p.4=list(formula=~pre94:(sex*agep94) + pre94:time + p94plus:time:sex:agep + p94plus:prepro + p94plus:male:prepro + p94plus:area2+ p94plus:male:area2 + p94plus:area1 + p94plus:male:area1,fixed=fixed.0,remove.intercept=TRUE)
# for time <1994, sex:agep and for >= 1994, time:sex:agep + sex:reproductive + sex:agep for area 1 and 2
	p.5=list(formula=~pre94:(sex*agep94)  + pre94:time + p94plus:time:sex:agep + p94plus:prepro + p94plus:male:prepro + p94plus:sex:agep:area1 + p94plus:sex:agep:area2, fixed=fixed.0,remove.intercept=TRUE)
# for time <1994, sex:agep and for >= 1994, time:sex:agep + sex:reproductive by time due to ENSO events but additive for males + sex:agep for area 1 and 2
	p.6=list(formula=~pre94:(sex*agep94)  + pre94:time + p94plus:time:sex:agep + p94plus:time:prepro + p94plus:male:prepro + p94plus:sex:agep:area1 + p94plus:sex:agep:area2, fixed=fixed.0,remove.intercept=TRUE)
# for time <1994, sex:agep and for >= 1994, time:sex:agep + sex:reproductive:time due to ENSO events + sex:agep for area 1 and 2
	p.7=list(formula=~pre94:(sex*agep94)  + pre94:time + p94plus:time:sex:agep + p94plus:time:prepro + p94plus:male:time:prepro + p94plus:sex:agep:area1 + p94plus:sex:agep:area2, fixed=fixed.0,remove.intercept=TRUE)
# for time <1994, sex:agep and for >= 1994, time:sex:agep + sex:reproductive by time due to ENSO events but additive for males + sex:agep for area 1 and 2 + time:sex:agep for area 1
	p.8=list(formula=~pre94:(sex*agep94)  + pre94:time + p94plus:time:sex:agep + p94plus:time:prepro + p94plus:male:prepro + p94plus:sex:agep:area1:time + p94plus:sex:agep:area2, fixed=fixed.0,remove.intercept=TRUE)
#    Phi
	Phi.00=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex+AgeS))
	Phi.01=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS))
	Phi.02=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:time))
	Phi.03=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time))
	Phi.04=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time))
	Phi.05=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + two:time)))
	Phi.06=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + twoplus:time)))
	Phi.07=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + two:time + threeplus:time)))
	
	Phi.08=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + pup:weight))
	Phi.09=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight))
	Phi.10=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + two:time) + pup:weight))
	Phi.11=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + twoplus:time) + pup:weight))
	
	Phi.12=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + pup:weight + yearling:weight))
	Phi.13=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight))
	Phi.14=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + two:time) + pup:weight+ yearling:weight))
	Phi.15=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + twoplus:time) + pup:weight+ yearling:weight))
	
	Phi.16=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + pup:weight + yearling:weight + two:weight))
	Phi.17=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight + yearling:weight + two:weight))
	Phi.18=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + two:time) + pup:weight+ yearling:weight + two:weight))
	Phi.19=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:(yearling:time + twoplus:time) + pup:weight+ yearling:weight + two:weight))
	
	Phi.20=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA))
	Phi.21=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneSSTAnomalies))
	Phi.22=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneSSTAnomalies))
	
	Phi.23=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI))
	Phi.24=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI))
	
	Phi.25=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + repro + male:repro))
	Phi.26=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro))
	Phi.27=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro))
	
	Phi.28=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA +  male:repro + repro:AgeS))
	Phi.29=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.30=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + male:repro + repro:AgeS))
	
	Phi.31=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + repro:AgeS + male:repro:AgeS))
	Phi.32=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI +  repro:AgeS + male:repro:AgeS))
	Phi.33=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	
	Phi.34=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + two:male:JulytoJuneMEI +  repro:AgeS + male:repro:AgeS))
	Phi.35=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI +  repro:AgeS + male:repro:AgeS))
	Phi.36=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + two:male:JulytoJuneMEI +  male:repro + repro:AgeS))
	Phi.37=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI  + male:repro + repro:AgeS))
	
	
	Phi.38=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + two:male:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.39=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI +  repro:AgeS + male:repro:AgeS))
	Phi.40=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + two:male:JulytoJuneMEI +  male:repro + repro:AgeS))
	Phi.41=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI + male:repro + repro:AgeS))
	
	
	Phi.42=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + threeplus:JulytoJuneMEI +  male:repro + repro:AgeS))
	Phi.43=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + two:male:JulytoJuneMEI + threeplus:JulytoJuneMEI + threeplus:male:JulytoJuneMEI +  male:repro + repro:AgeS))
	Phi.44=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.45=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + two:male:JulytoJuneMEI + threeplus:JulytoJuneMEI + threeplus:male:JulytoJuneMEI + male:repro + repro:AgeS))
	
	
	Phi.46=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI  + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.47=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.48=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.49=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.50=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
	
	Phi.51=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI  + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.52=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.53=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.54=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.55=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
	
	Phi.56=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time  + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI +  male:repro + repro:AgeS))
	Phi.57=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI  + male:repro + repro:AgeS))
	Phi.58=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.59=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.60=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
	
	Phi.61=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies  + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.62=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI +  male:repro + repro:AgeS))
	Phi.63=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.64=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.65=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
	
	Phi.66=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies  + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.67=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.68=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI +  repro:AgeS + male:repro:AgeS))
	Phi.69=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.70=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoJuneSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
	
	Phi.71=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time  + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.72=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + male:repro + repro:AgeS))
	Phi.73=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.74=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro:AgeS + male:repro:AgeS))
	Phi.75=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + male:repro + repro:AgeS))
	model.list=create.model.list(c("Phi","p"))
	models=function(x)
	{
		model.parameters=list()
		for(j in 1:length(x))
		{
			if(!is.list(eval(parse(text=x[j]),envir=environment(models))[[1]]))
				model.parameters[[names(x)[j]]]=eval(parse(text=(as.character(x[j]))),envir=environment(models))
		}
		for(j in 1:length(x))
		{
			if(is.list(eval(parse(text=x[j]),envir=environment(models))[[1]]))
				model.parameters=c(model.parameters,eval(parse(text=(as.character(x[j]))),envir=environment(models)))
		}
		model.name=paste(x,collapse=".")
		return(list(model.name=model.name,model.parameters=model.parameters))
		
	}
	return(models)
}
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
	model$vcv=vcv*1.11
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





