\name{TSEN}
\alias{TSEN}
\title{
Sentinel target peaks from two-dimensional chromatogram
}
\description{
This function assign and integrate target peaks in two-dimentional chromatogram. 
}
\usage{
TSEN(TSENdata, database, countMS, phase.shift, search.factor, HR.wide.order, HR.wide.order.fluc, MPeriod, MSthrethold,inicutRT, endcutRT )
}
\arguments{
	\item{TSENdata}{Name of the netCDF file which containing GC x GC-HRTOFMS data.}
	\item{database}{Name of the csv file which containing list of target compounds for search.}
	\item{countMS}{The number of MS counting for searching/assignment of target compound.}
	\item{phase.shift}{The degree of modification of start signal timing.}
	\item{search.factor}{The degree of search range of target compound in the chromatogram.}
	\item{HR.wide.order}{The decimal point position in accurate mass which is used for assingment.}
	\item{HR.wide.order.fluc}{The tolerable degree of m/z fluctuation in the last digit which decided by the setting of HR.wide.order.}
	\item{MPeriod}{The modulation time period.}
	\item{MSthrethold}{The minimum threshold of intensity of ions for searching/assignment of target compound.}
	\item{inicutRT}{The setting of initial time (min) for conducting data analysis. It should be smaller than the (RT of target compounds) but larger than minimum RT in the chromatogram.}
	\item{endcutRT}{The setting of cut time (min) from end RT for conducting data analysis. It should be smaller than (end RT) -  (the RT of target compounds).}
}
\value{
The values are summarised in csv file created in "Document" folder (working directory). The jpeg file of TIC and combined EIC are also created.
}

\author{
Yasuyuki ZUSHI
}


\examples{
library(ncdf)
data(database.dat)
cdf.dat <- 1

#Beware that file will be downloaded from the Internet, if executing the following 3 lines. Check that you are online.

#web_download1 <- "http://www.nies.go.jp/analysis/downloads/demo_data_stdPCDD_DFs.cdf"
#download.file(url=web_download1,destfile="demo_data_stdPCDD_DFs.cdf",mode="wb")
#cdf.dat <- open.ncdf("demo_data_stdPCDD_DFs.cdf",write=TRUE )


TSEN(TSENdata=cdf.dat,database=database.dat,countMS=1,HR.wide.order=1,HR.wide.order.fluc=1, inicutRT=0, endcutRT=0 ) #The wide window of m/z = 0.1 * 1 = 0.1

}
\keyword{TSEN}