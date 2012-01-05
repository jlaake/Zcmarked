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

Phi.28=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + repro + male:repro + repro:RepAge))
Phi.29=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.30=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))

Phi.31=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.32=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.33=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))

Phi.34=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + two:male:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.35=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.36=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + two:male:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.37=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI + repro + male:repro + repro:RepAge))


Phi.38=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + two:male:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.39=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.40=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + two:male:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.41=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + twoplus:JulytoJuneMEI + twoplus:male:JulytoJuneMEI + repro + male:repro + repro:RepAge))


Phi.42=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.43=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + two:JulytoJuneMEI + two:male:JulytoJuneMEI + threeplus:JulytoJuneMEI + threeplus:male:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.44=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.45=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + two:male:JulytoJuneMEI + threeplus:JulytoJuneMEI + threeplus:male:JulytoJuneMEI + repro + male:repro + repro:RepAge))


Phi.46=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI  + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.47=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.48=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.49=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.50=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))

Phi.51=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI  + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.52=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.53=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.54=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.55=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarMEI + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))

Phi.56=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time  + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.57=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.58=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.59=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.60=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneMEI + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))

Phi.61=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies  + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.62=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.63=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.64=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.65=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + p90plus:yearling:time + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))

Phi.66=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies  + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.67=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.68=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.69=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.70=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + pup:I(Time-12) + pup:OcttoMarSSTAnomalies + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))

Phi.71=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time  + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.72=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
Phi.73=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + two:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.74=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA + twoplus:JulytoJuneMEI + repro + male:repro + repro:RepAge + male:repro:RepAge))
Phi.75=list(formula=~-1+f1987+f1988+m1987+m1988+post:(sex*AgeS + p89plus:pup:time + yearling:I(Time-12) + yearling:JulytoJuneSSTAnomalies + pup:weight+ yearling:weight + male:twoplus:leptom + female:twoplus:DA+ two:JulytoJuneMEI + threeplus:JulytoJuneMEI + repro + male:repro + repro:RepAge))
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





