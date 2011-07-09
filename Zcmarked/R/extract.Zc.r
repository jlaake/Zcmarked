

#' Extracts Zc SMI capture histories and covariates
#' : from ACCESS database, it constructs the relevant capture histories from the
#' branding(ZCBrand) and resighting (Alive) tables. It constructs all the
#' queries that used to be done in ACCESS.
#' 
#' Attaches directly to ACCESS database and extracts initial brand and all
#' resighting information.  Excludes any sea lions with missing sex or initial
#' weight. Also, excludes the small number of animals branded in the spring and
#' only uses those branded past August (month 8).  Restricts resightings to
#' those seen between \code{begin} and \code{end} inclusive of those dates.
#' Also restricts resightings to SMI and ANI from 1990 onwards and to FAI/FAR
#' from 1992 onwards. SMI resightings are split into 2 areas based on field
#' \code{two} in table smiareacodes.  All other resightings are lumped into
#' area 3. A sea lion is only considered resighted for the year if either it
#' was seen twice in the year or a photo was taken.
#'
#' @import RODBC 
#' @export
#' @param file ACCESS database filename
#' @param dir Directory containing ACCESS database
#' @param begin month-day at beginning of resight period (515 = May 15)
#' @param end month-day at end of resight period (815 = August 15)
#' @param select either 0 or 1; if 1 it requires at least 2 resightings of a
#'   sea lion to be consider a recapture
#' @param lastyear final year to be included in resight history not final cohort which 
#'   is currently hard-coded at 2009
#' @return dataframe containing following fields \item{ch}{capture history;
#'   character string} \item{cohort}{year branded; factor variable}
#'   \item{sex}{either M or F; factor variable} \item{weight}{weight (kg) at
#'   time of branding} \item{brand}{brand identifier} \item{area1}{sequence of
#'   fields named area1yyyy with values 0/1; it is 1 if resighted in a year
#'   prior to yyyy in area1; can be used as a time-varying covariate}
#'   \item{area2}{sequence of fields named area2yyyy with values 0/1; it is 1
#'   if resighted in a year prior to yyyy in area2 but not seen previously in
#'   area 1; if not area1 or area2 it is either not been resighted or resighted
#'   in area 3 (not SMI)} \item{repro}{sequence of fields named reproyyyy with
#'   values 0/1; it is 1 if seen as reproductive in a year prior to yyyy}
#'   \item{td}{sequence of fields named tdyyyy with values 0/1; it is 1 if seen
#'   in year yyyy-1} \item{recap}{0 if never resighted and 1 otherwise}
#'   \item{TotalTimesResighted}{number of years resighted}
#' 
#' @author Jeff Laake
#' @examples 
#' zcdata=extract.Zc()
extract.Zc <-
		function(file="BrandMaster.mdb",dir="",begin=515,end=815,select=1,lastyear=2010)
{
	if(dir=="")dir=system.file(package="Zcmarked")
	fdir=file.path(dir,file)
	connection=odbcConnectAccess(fdir)
	Alive=sqlFetch(connection,"Alive")
	Alive$pupyear=as.numeric(as.POSIXlt(Alive$sitedate)$year)+1900
	Alive=droplevels(Alive[Alive$pupyear<=lastyear,])
	Brand=sqlFetch(connection,"ZcBrand")
	mon=as.POSIXlt(Brand$branddate)$mon+1
# Only use brands with known sex and weight and branded in Sept or later and cohort <= 2009
# This code will no longer work with inclusion of V brands in 2010 and will have to be changed; 2009 hard coded and 
# does not work for lesser value either but code could be modified to adjust cohort.count.table
	Brand=droplevels(Brand[Brand$sex!="U" & !is.na(Brand$weight) & Brand$weight>0 & mon>8 & Brand$cohort<=2009,])
	Brand$sex=Brand$adjsex
	AreaCodes=sqlFetch(connection,"smiareacodes")
# Limited Resights
	xx=merge(Alive,AreaCodes[,c("code","two")],all.x=TRUE,by="code")
	xx=xx[(xx$sitecode%in%c("ANI","SMI")&xx$pupyear>=1990) | (xx$sitecode%in%c("FAI","FAR")&xx$pupyear>=1992),]
	mday=as.numeric(as.POSIXlt(xx$sitedate)$mon+1)*100+as.POSIXlt(xx$sitedate)$mday
	LimitedResights=xx[mday<=end& mday>=begin,]
	LimitedResights$capregion=ifelse(is.na(LimitedResights$two),3,LimitedResights$two)
    LimitedResights=merge(LimitedResights,Brand[,c("brand","sex","cohort")],by="brand")
	LimitedResights$adjsex=LimitedResights$sex
	LimitedResights$sex=NULL
	LimitedResights=LimitedResights[LimitedResights$cohort<=2009,]
	LimitedResights$cohort=NULL
	LimitedResights$pregnant=as.character(LimitedResights$pregnant)
	LimitedResights$withpup=as.character(LimitedResights$withpup)	
	LimitedResights$pregnant[is.na(LimitedResights$pregnant)]=""
	LimitedResights$withpup[is.na(LimitedResights$withpup)]=""
	LimitedResights$pregnant=factor(LimitedResights$pregnant)
	LimitedResights$withpup=factor(LimitedResights$withpup)	
	LimitedResights$repro=ifelse((LimitedResights$adjsex=="F"&(LimitedResights$withpup%in%c("Y","L","J","y","l","j")|LimitedResights$pregnant%in%c("p","P")))|
					             (LimitedResights$adjsex=="M"&LimitedResights$withpup%in%c("T","t")),1,0)
	LimitedResights=droplevels(LimitedResights)
# BrandResightJoin
	BrandResightJoin=merge(Brand,LimitedResights,by="brand",all.x=TRUE)
	if(lastyear>2009)BrandResightJoin$cohort=factor(BrandResightJoin$cohort,levels=c(sort(unique(BrandResightJoin$cohort)), 2010:lastyear))
# LimitedResightswithCount
	brandyear=paste(BrandResightJoin$brand,BrandResightJoin$pupyear,sep="")
	count.table=table(brandyear)
	xx=data.frame(brandyear=names(count.table),count=as.numeric(count.table))
	LimitedResights$brandyear=paste(LimitedResights$brand,LimitedResights$pupyear,sep="")
	xx=merge(LimitedResights,xx,by="brandyear")
	LimitedResightswithCount=xx[xx$count>1,]
# CaptureHistory
	resight.count.table=with(BrandResightJoin,table(brand,pupyear))
	resight.count.table=cbind(matrix(0,ncol=3,nrow=nrow(resight.count.table)),resight.count.table)
	cohort.count.table=with(BrandResightJoin,table(brand,cohort))
	resight.count.table=ifelse(resight.count.table<=select,0,1)
	cohort.count.table=ifelse(cohort.count.table<=0,0,1)
	BrandResightJoin$nphotos=ifelse(BrandResightJoin$photos%in%c("","N")|is.na(BrandResightJoin$photos),0,1)
	photo.table=with(BrandResightJoin,tapply(nphotos,list(brand,pupyear),sum))
	photo.table[is.na(photo.table)]=0
	photo.table=cbind(matrix(0,ncol=3,nrow=nrow(photo.table)),photo.table)
	capture.history=cohort.count.table+resight.count.table+photo.table
	capture.history[capture.history>1]=1
	xx=Brand[,c("brand","sex","cohort","weight")]
	CaptureHistory=as.data.frame(capture.history)
	CaptureHistory$brand=as.numeric(row.names(CaptureHistory))
	CaptureHistory=merge(xx,CaptureHistory,by="brand",all.x=TRUE)
	capture.history=CaptureHistory[,-(1:4)]
	CaptureHistory$TotalTimesResighted=rowSums(capture.history)-1
	CaptureHistory$recap=ifelse(CaptureHistory$TotalTimesResighted>0,1,0)
# ReproCovariates
	xx=with(LimitedResightswithCount[LimitedResightswithCount$repro==1,],tapply(pupyear,brand,min))
	repro.table=data.frame(brand=as.numeric(names(xx)),repro.year=as.numeric(xx))
#	repro.table$repro.year=ifelse(repro.table$repro.year<100,1900+repro.table$repro.year,repro.table$repro.year)
	repro.table=with(repro.table,table(brand,repro.year))
	xx=as.matrix(repro.table)
	class(xx)="matrix"
	xx=cbind(data.frame(brand=as.numeric(row.names(repro.table))),xx)
	xx=merge(Brand[,"brand",drop=FALSE],xx,all.x=TRUE,by="brand")
	xx=xx[,-1]
	xx[is.na(xx)]=0
	xx[]=t(apply(xx,1,cumsum))
	ReproCovariates=xx
	names(ReproCovariates)=paste("repro",names(ReproCovariates),sep="")
# RegionCovariates
	xx=with(LimitedResightswithCount[LimitedResightswithCount$capregion==1,],tapply(pupyear,brand,min))
	area1.table=data.frame(brand=as.numeric(names(xx)),area1.year=as.numeric(xx))
#	area1.table$area1.year=ifelse(area1.table$area1.year<100,1900+area1.table$area1.year,area1.table$area1.year)
	area1.table=with(area1.table,table(brand,area1.year))
	xx=as.matrix(area1.table)
	class(xx)="matrix"
	xx=cbind(data.frame(brand=as.numeric(row.names(area1.table))),xx)
	xx=merge(Brand[,"brand",drop=FALSE],xx,all.x=TRUE,by="brand")
	xx=xx[,-1]
	xx[is.na(xx)]=0
	xx[]=t(apply(xx,1,cumsum))
	area1.table=xx
	xx=with(LimitedResightswithCount[LimitedResightswithCount$capregion==2,],tapply(pupyear,brand,min))
	area2.table=data.frame(brand=as.numeric(names(xx)),area2.year=as.numeric(xx))
#	area2.table$area2.year=ifelse(area2.table$area2.year<100,1900+area2.table$area2.year,area2.table$area2.year)
	area2.table=with(area2.table,table(brand,area2.year))
	xx=as.matrix(area2.table)
	class(xx)="matrix"
	xx=cbind(data.frame(brand=as.numeric(row.names(area2.table))),xx)
	xx=merge(Brand[,"brand",drop=FALSE],xx,all.x=TRUE,by="brand")
	xx=xx[,-1]
	xx[is.na(xx)]=0
	xx[]=t(apply(xx,1,cumsum))
	area2.table=xx
	area2.table=area2.table*(1-area1.table)
	names(area1.table)=paste("area1",as.numeric(names(area1.table))+1,sep="")
	names(area2.table)=paste("area2",as.numeric(names(area2.table))+1,sep="")
	nc=ncol(area1.table)
	RegionCovariates=cbind(area1.table[,-nc],area2.table[,-c(1,nc)])
# MarkData
	CaptureHistory$ch=apply(capture.history,1,paste,collapse="")
	td=as.data.frame(capture.history[-ncol(capture.history)])
	td=td-table(CaptureHistory$brand,CaptureHistory$cohort)
	names(td)=paste("td",as.numeric(names(td))+1,sep="")
	MarkData=CaptureHistory[,c("ch","cohort","sex","weight","brand")]
	MarkData=cbind(MarkData,RegionCovariates,ReproCovariates,td,CaptureHistory[,c("recap","TotalTimesResighted")])
	MarkData$cohort=factor(MarkData$cohort)
	return(MarkData)
}