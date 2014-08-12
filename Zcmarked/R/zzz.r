.onAttach<-function(library,pkgname)
{
	# this now conforms with new R conventions
	# uses packageStartupMessage which can then be
	# surpressed
	print.Zcmarked.version()
}

# taken from mgcv
print.Zcmarked.version <- function(){ 
	library(help=Zcmarked)$info[[1]] -> linfo
	
	version <- linfo[pmatch("Version",linfo)]
	built <- linfo[pmatch("Built",linfo)]
	
	
	hello <- paste("This is Zcmarked ",version,sep="")
	packageStartupMessage(hello)
}
