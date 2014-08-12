#' Connection to brand database
#' 
#' Creates connection to database either in package or another directory
#' 
#' @import RODBC 
#' @export dbConnection
#' @param file ACCESS database filename
#' @param dir Directory containing ACCESS database
#' @return connection connection to database
dbConnection <-
		function(file="BrandMaster.mdb",dir="")
{
	if(dir=="")dir=system.file(package="CIPinnipedAnalysis")
	fdir=file.path(dir,file)
	connection=odbcConnectAccess(fdir)
	return(connection)
}
