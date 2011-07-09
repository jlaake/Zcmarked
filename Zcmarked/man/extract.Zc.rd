\name{extract.Zc}
\alias{extract.Zc}
\title{Extracts Zc SMI capture histories and covariates...}
\usage{extract.Zc(file="BrandMaster.mdb", dir="", begin=515,
    end=815, select=1, lastyear=2010)}

\description{
  Extracts Zc SMI capture histories and covariates : from
  ACCESS database, it constructs the relevant capture
  histories from the branding(ZCBrand) and resighting
  (Alive) tables. It constructs all the queries that used
  to be done in ACCESS.
}

\details{
  Attaches directly to ACCESS database and extracts initial
  brand and all resighting information.  Excludes any sea
  lions with missing sex or initial weight. Also, excludes
  the small number of animals branded in the spring and
  only uses those branded past August (month 8).  Restricts
  resightings to those seen between \code{begin} and
  \code{end} inclusive of those dates. Also restricts
  resightings to SMI and ANI from 1990 onwards and to
  FAI/FAR from 1992 onwards. SMI resightings are split into
  2 areas based on field \code{two} in table smiareacodes.
  All other resightings are lumped into area 3. A sea lion
  is only considered resighted for the year if either it
  was seen twice in the year or a photo was taken.
}
\value{dataframe containing following fields \item{ch}{capture history;
character string} \item{cohort}{year branded; factor variable}
\item{sex}{either M or F; factor variable} \item{weight}{weight (kg) at
time of branding} \item{brand}{brand identifier} \item{area1}{sequence of
fields named area1yyyy with values 0/1; it is 1 if resighted in a year
prior to yyyy in area1; can be used as a time-varying covariate}
\item{area2}{sequence of fields named area2yyyy with values 0/1; it is 1
if resighted in a year prior to yyyy in area2 but not seen previously in
area 1; if not area1 or area2 it is either not been resighted or resighted
in area 3 (not SMI)} \item{repro}{sequence of fields named reproyyyy with
values 0/1; it is 1 if seen as reproductive in a year prior to yyyy}
\item{td}{sequence of fields named tdyyyy with values 0/1; it is 1 if seen
in year yyyy-1} \item{recap}{0 if never resighted and 1 otherwise}
\item{TotalTimesResighted}{number of years resighted}}
\author{Jeff Laake}
\arguments{
  \item{file}{ACCESS database filename}
  \item{dir}{Directory containing ACCESS database}
  \item{begin}{month-day at beginning of resight period (515 = May 15)}
  \item{end}{month-day at end of resight period (815 = August 15)}
  \item{select}{either 0 or 1; if 1 it requires at least 2 resightings of a
sea lion to be consider a recapture}
  \item{lastyear}{final year to be included in resight history not final cohort which 
is currently hard-coded at 2009}
}
\examples{zcdata=extract.Zc()}
