\name{EnvironCovariates}
\alias{EnvironCovariates}
\title{Extracts environmental covariate data for SMI Zc analysis
: by extracting data from ACCESS database and creating a series of
anomaly summaries for sea surface temperature(SST), upwelling index (UWI) and
multivariate ENSO index (MEI).}
\usage{
  EnvironCovariates(average.years = c(1994:1996, 1998:2008),
    fdir = NULL, sites = 1:5)
}
\arguments{
  \item{average.years}{years to use for average for anomaly
  creation}

  \item{fdir}{directory that contains
  environmental.data.mdb; NULL is std location}

  \item{sites}{SST sites to be used in SST averaging (from
  create.SST.anomalies) 1) ESB, 2)WSB, 3)PtArg, 4)PtSM,
  5)PtSL, 6)CSM, 7)MB, 8)PtReyes}
}
\value{
  dataframe with rows as years from 1987 to last year in
  data and the columns are:
  time,OcttoJuneSSTAnomalies,ApriltoJuneSSTAnomalies,JulytoJuneSSTAnomalies,
  OcttoJuneUWI33Anomalies, ApriltoJuneUWI33Anomalies,
  JulytoJuneUWI33Anomalies, OcttoJuneUWI36Anomalies,
  ApriltoJuneUWI36Anomalies, JulytoJuneUWI36Anomalies,
  OcttoJuneMEI, ApriltoJuneMEI, JulytoJuneMEI
}
\description{
  Extracts environmental covariate data for SMI Zc analysis
  : by extracting data from ACCESS database and creating a
  series of anomaly summaries for sea surface
  temperature(SST), upwelling index (UWI) and multivariate
  ENSO index (MEI).
}
\author{
  Jeff Laake
}

