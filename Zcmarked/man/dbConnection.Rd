\name{dbConnection}
\alias{dbConnection}
\title{Connection to brand database}
\usage{
  dbConnection(file = "BrandMaster.mdb", dir = "")
}
\arguments{
  \item{file}{ACCESS database filename}

  \item{dir}{Directory containing ACCESS database}
}
\value{
  connection connection to database
}
\description{
  Creates connection to database either in package or
  another directory
}

