\name{define.models}
\alias{define.models}
\title{Define model set used in Zc survival analysis}
\usage{
  define.models(zcdata.proc, maxyear)
}
\arguments{
  \item{zcdata.proc}{processed dataframe for Zc analysis}

  \item{maxyear}{Last year used in analysis}
}
\value{
  function that can be used to select a subset of models to
  be used in the analysis
}
\description{
  Creates and returns a function that defines p and Phi
  model specifications for survival analysis.
}
\author{
  Jeff Laake
}

