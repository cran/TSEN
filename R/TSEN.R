TSEN <-
function(TSENdata=1,database=1, countMS=5, phase.shift=-1.3, search.factor=1.5, HR.wide.order=0.01, HR.wide.order.fluc=5, MPeriod=4, MSthrethold=0, inicutRT=0, endcutRT=0 ){

defaultSamRate <- 16
MSmode <- "HRscan"



if (class(TSENdata)=="ncdf"){
flag <- get.var.ncdf(TSENdata,"flag_count")
scantime <- get.var.ncdf(TSENdata,"scan_acquisition_time")
scannum <- get.var.ncdf(TSENdata,"actual_scan_number")
medmsmax <- get.var.ncdf(TSENdata,"mass_range_max")
medmsmin <- get.var.ncdf(TSENdata,"mass_range_min")
ionid <- get.var.ncdf(TSENdata,"scan_index")
eachscannum <- get.var.ncdf(TSENdata,"point_count")
MStotint <- get.var.ncdf(TSENdata,"total_intensity")
MSvalue <- get.var.ncdf(TSENdata,"mass_values") 
MSint <- get.var.ncdf(TSENdata,"intensity_values")

} else {


scantime <- c(10, 20, 30)
scannum <- c(10, 20, 30)
eachscannum <- c(10, 20, 30)
MStotint <- c(10, 20, 30)
MSvalue <- c(10, 20, 30)
MSint <- c(10, 20, 30)
}


Timepara <- scantime[which(abs(eachscannum)<.Machine$integer.max )]
RTini <- min( Timepara )/60
RTruntime <- max( Timepara )/60-RTini
SamRate <- 1/(mean( Timepara[2:length(Timepara)]-Timepara[1: (length(Timepara)-1) ] ))


pixelnum <- length(scannum)-round(endcutRT*SamRate*60)
maxscannum <- max(eachscannum)
minscannum <- min(eachscannum)
MSdatabox <- matrix(0,pixelnum, maxscannum+1 )
MSvaluebox <- MSdatabox
MSintbox <- MSdatabox

inum.range <- round((inicutRT)*SamRate*60+1):pixelnum
for (inum in inum.range) {
if (abs(eachscannum[inum])<.Machine$integer.max){
initial <- sum(eachscannum[1:inum])-eachscannum[inum]+1
acqrange <- initial:(initial+eachscannum[inum]-1)
if(eachscannum[inum]==0){acqrange <-0}
remainrep <- rep(0,maxscannum-eachscannum[inum]+1 )
MSvaluebox[inum,] <- c( MSvalue[ acqrange ], remainrep )
MSintbox[inum,] <- c( MSint[ acqrange ], remainrep )
}
}

rownames(MSvaluebox) <-1:(pixelnum)
rownames(MSintbox) <-1:(pixelnum)
MSintbox[which(MSintbox < MSthrethold)] <-0
MSvaluebox[which(MSintbox < MSthrethold)] <-0
MSintbox <- replace(MSintbox,is.na(MSintbox),1)
MSvaluebox <- replace(MSvaluebox,is.na(MSvaluebox),1)

  
Count.MS <- countMS
Count.MS.seq <- 17:(17+Count.MS-1)
any.or.all <- "&&" 
 
if (any.or.all=="&&"){
area.get <- "st.sum"
}
if (any.or.all=="||"){
area.get <- "st.sum/Sum.MS"
}

SumMSnum <- countMS
Sum.MS <- SumMSnum
Sum.MS.seq <- 17:(17+Sum.MS-1) 

if(Count.MS==0){
use.Count.MS <- "nonMS"
}else if (Count.MS==1){
use.Count.MS <- "SingleMS"
}else if (Count.MS>=2){
use.Count.MS <- paste(Count.MS,"MS",sep="")
}


database[,7] <-database[,7]-phase.shift

if( any(database[seq(1,length(database[,1]),2),7]>MPeriod) ){
database[which(database[,7]>MPeriod),7] <- database[which(database[,7]>MPeriod),7] - MPeriod
}
database.origin <- database

row <- round(SamRate*MPeriod)
col <- round(length(scantime)/row)
if (row==0){row<-2
col<-2}


wide.op <-0

if (MSmode=="HRscan"){
widefactor.op <- SamRate/defaultSamRate 
HRfactor <- 1
HRcell <- 1
LRcell <- 0
} else {
HRfactor <- 1
widefactor.op <-1
HRcell <- 0
LRcell <- 1
}

heading <- LRcell
tailing <- LRcell
headbuffer <- round(heading )  
tailbuffer <- round(tailing ) 
serch.range <- round(1 *widefactor.op + wide.op )
serch.range2 <- round(8 *widefactor.op )
noise.range <- round(((row-(serch.range2*2+1))/2))/2 
ColRange <- 3 *HRfactor +15*HRcell + wide.op
RowRange <- round(7 *HRfactor*widefactor.op +7*HRcell + wide.op ) 
addsearch <- 3 *HRfactor+15*HRcell + wide.op  
ColRange.f <- ColRange
RowRange.f <- RowRange
addsearch.f <- addsearch 
third.colrange <- heading + search.factor *HRcell  
third.rowrange <- round(2 *widefactor.op + search.factor*HRcell + wide.op )
third.addserch <- tailing + (third.colrange - heading ) 



MScalib <- 1


if (MSmode=="HRscan"){ 
MScalib <- MScalib*HR.wide.order*HR.wide.order.fluc 
}
if (MSmode=="HRscan"){ 
keta <- ceiling(abs(log10(MScalib))) 
} else { 
keta <-0 
}
if (MSmode=="HRscan") {
MScalib.seq <- seq((-MScalib),MScalib,HR.wide.order )
} else { 
MScalib.seq <- seq((-MScalib),MScalib,MScalib ) 
}


targetcol.ini.set <-0
targetrow.ini.set <-0

for (h in 1:(length(database[,1])/2) ){
targetcol.ini <- as.matrix(round((database[2*h-1,6]-RTini)*60/MPeriod)- round( ColRange.f ) )
if( targetcol.ini <0 ){ 
targetcol.ini <- replace(targetcol.ini,targetcol.ini <1,1) 
} 
if ((round((database[2*h-1,7])/MPeriod*row)-RowRange.f)>0){
targetrow.ini <- as.matrix(round((database[2*h-1,7])/MPeriod*row)-RowRange.f)  
if(targetrow.ini > row){
targetrow.ini <- as.matrix(row)
}
if(targetrow.ini < 1){
targetrow.ini <- as.matrix(1)
}
} else{
targetrow.ini <- as.matrix(1)}
rownames(targetcol.ini) <- database[2*h-1,2]
rownames(targetrow.ini) <- database[2*h-1,2]
targetcol.ini.set <- rbind(targetcol.ini.set,targetcol.ini)
targetrow.ini.set <- rbind(targetrow.ini.set,targetrow.ini)
}

targetcol.ini.set <- as.matrix(targetcol.ini.set[-c(1),])
targetrow.ini.set <- as.matrix(targetrow.ini.set[-c(1),])

MS.Fragment <- paste("any(frag.i[m,]==round(database[2*i-1,",Count.MS.seq[1], "]+(",MScalib.seq,"),digits=keta))||", sep="", collapse="")  
MS.Fragment <- paste("(",substr(MS.Fragment,1,(nchar(MS.Fragment)-2)),")")

if (length(Count.MS.seq) >1 ){
for (isn in 2:length(Count.MS.seq) ){
MS.Fragment.pre <- paste("any(frag.i[m,]==round(database[2*i-1,",Count.MS.seq[isn], "]+(",MScalib.seq,"),digits=keta))||", sep="", collapse="")  
MS.Fragment.pre <- substr(MS.Fragment.pre,1,(nchar(MS.Fragment.pre)-2))
MS.Fragment <- paste( MS.Fragment,any.or.all,"(",MS.Fragment.pre,")" )
}
}

MS.ext <- MS.Fragment


Sum.MS.ratio <- paste("replace( obj.mat[m,which(  frag.i[m,]==round(database[2*i-1,",Sum.MS.seq[1],"]+(",MScalib.seq,"),digits=keta) )][1],is.na(obj.mat[m,which(  frag.i[m,]==round(database[2*i-1,",Sum.MS.seq[1],"]+(",MScalib.seq,"),digits=keta) )][1]) ,0 ),  ",collapse="",sep="")
Sum.MS.ratio <- paste( "max(",substr(Sum.MS.ratio,1,(nchar(Sum.MS.ratio)-3)),")/database[2*i,Sum.MS.seq[1]]",sep="")

if (length(Sum.MS.seq)>1 ){
for (isn in 2:length(Sum.MS.seq) ){
Sum.MS.ratio.pre <- paste( "replace( obj.mat[m,which(  frag.i[m,]==round(database[2*i-1,",Sum.MS.seq[isn],"]+(",MScalib.seq,"),digits=keta) )][1],is.na(obj.mat[m,which(  frag.i[m,]==round(database[2*i-1,",Sum.MS.seq[isn],"]+(",MScalib.seq,"),digits=keta) )][1]) ,0 ),",collapse="",sep="")
Sum.MS.ratio.pre <- paste("max(",substr(Sum.MS.ratio.pre,1,(nchar(Sum.MS.ratio.pre)-1)),")/database[2*i,Sum.MS.seq[isn]]")
Sum.MS.ratio <- paste(  Sum.MS.ratio,",",Sum.MS.ratio.pre )
}
}

Sum.MS.ratio.ext <- Sum.MS.ratio
res.case <- matrix(0,pixelnum,1)

rownames(res.case) <- 1:pixelnum
res.box0 <- matrix(0,pixelnum,1)


result.box <- matrix(0,(length(database[,1])/2) ,10)
colnames(result.box) <- c("Compound","Volume","SN","Noise candidate","Noise number","RT1","RT2", "peak top ID","IS/Native","ID of IS")
oneline.box.num <- matrix(1:pixelnum,pixelnum,1)
icol.range <- round(ColRange+addsearch+1 + headbuffer + tailbuffer )


database.firstIS <- database
ColRange.firstIS <- third.colrange
RowRange.firstIS <- third.colrange
addsearch.firstIS <- third.addserch


redidu <- row-SamRate*MPeriod
redidu.accum <-0
blank.id <- round(phase.shift*SamRate+1)
cont.num<-0
tic.box <- as.matrix( MStotint  )
mat <- matrix(1,row,col)

try(for (u in round(phase.shift*SamRate+1):(length(tic.box[,1])+round(phase.shift*SamRate)) ){
if (u%%row==0){
redidu.accum <-redidu.accum + redidu
}
if ( redidu.accum > 1  ){
blank.id <- blank.id +1
mat[row*ceiling(blank.id/row)-(blank.id-1),ceiling( (u-round(phase.shift*SamRate))/row )] <- tic.box[round(u-phase.shift*SamRate),]
blank.id <- blank.id +1
mat[row*ceiling(blank.id/row)-(blank.id-1),ceiling( (u-round(phase.shift*SamRate))/row )] <- tic.box[round(u-phase.shift*SamRate),]
redidu.accum <- redidu.accum -1 
} else {
blank.id <- blank.id +1
mat[row*ceiling(blank.id/row)-(blank.id-1),ceiling( (u-round(phase.shift*SamRate))/row )] <- tic.box[round(u-phase.shift*SamRate),]
}
})

zz1 <- t(apply(as.matrix(mat[1:nrow(mat),1:ncol(mat)]),2,rev))
xx1 <- (1:nrow(zz1))/col*RTruntime+RTini
yy1 <-(((1:ncol(zz1))/row*MPeriod))
max.zz<-max(zz1)

dev.new(width=8,height=4)
par(mar=c(5,7,1,1))
filled.contour(xx1,yy1,zz1,xlab="GC1 (min)",ylab="GC2 (sec)",col=terrain.colors(30),levels=seq(0,max(tic.box)+0.1,max(tic.box)/30+0.1) )
text(20,3.8,"TIC")
dev.copy(jpeg,paste("TIC.jpeg",sep=""),width=960,height=480)
dev.off()
dev.off()


try(for (i in 1:(length(database[,1])/2)  ){ 

cand.col <- matrix(0, row,icol.range  )
oneline.box <- matrix(0,pixelnum,1)
rownames(oneline.box) <- 1:length(oneline.box)
original.max.box2 <- oneline.box

if ( database[2*i-1,4]==1 ){
ColRange <- ColRange.f
RowRange <- RowRange.f
addsearch <- addsearch.f

icol.range <- round( ColRange+addsearch+1 + headbuffer + tailbuffer ) 
cand.col <- matrix(0, row,icol.range  )
}


if ( (database[2*i-1,5]==0) & (database[which(database[,3]==database[2*i-1,3]),4][1]==1) ){
if(database[2*i-3,5]==1){ 
database <- database.firstIS
ColRange <- ColRange.firstIS 
RowRange <- RowRange.firstIS 
addsearch <- addsearch.firstIS 


targetcol.ini.set <-0
targetrow.ini.set <-0

for (h in 1:(length(database[,1])/2) ){
targetcol.ini <- as.matrix(round((database[2*h-1,6]-RTini)*60/MPeriod)-round( ColRange) )
if( targetcol.ini <0 ){ 
targetcol.ini <- replace(targetcol.ini,targetcol.ini <1,1) 
} 
if ((round((database[2*h-1,7])/MPeriod*row)-RowRange)>0){
targetrow.ini <- as.matrix(round((database[2*h-1,7])/MPeriod*row)-RowRange)
if(targetrow.ini > row){
targetrow.ini <- as.matrix(row)
}
if(targetrow.ini < 1){
targetrow.ini <- as.matrix(1)
}
} else{
targetrow.ini <- as.matrix(1)}
rownames(targetcol.ini) <- database[2*h-1,2]
rownames(targetrow.ini) <- database[2*h-1,2]
targetcol.ini.set <- rbind(targetcol.ini.set,targetcol.ini)
targetrow.ini.set <- rbind(targetrow.ini.set,targetrow.ini)
}

targetcol.ini.set <- as.matrix(targetcol.ini.set[-c(1),])
targetrow.ini.set <- as.matrix(targetrow.ini.set[-c(1),])

 
icol.range <- round( ColRange+addsearch+1 + headbuffer + tailbuffer  )
cand.col <- matrix(0, row,icol.range  ) 
}
}


for (j in 1:(icol.range) ){
cand.col[,j] <- as.matrix( ((targetcol.ini.set[i])*row+1):( ((targetcol.ini.set[i]+1))*row)  ) + row*(j-1) 
}
if(row==2){cand.col<-matrix(1,2,2)}

obj.mat <- MSintbox[cand.col,]
obj.mat.ms <- MSvaluebox[cand.col,]
if(row==2){rownames(obj.mat)<-1:4
rownames(obj.mat.ms)<-1:4
}


st.t <- matrix(0,nrow(cand.col),ncol(cand.col))
oneline <- matrix(0,length(cand.col),1)
rownames(oneline) <- rownames(obj.mat[1:(length(cand.col)),])

obj.mat.line <- 1:length(obj.mat[1,]) 


if(Count.MS==0){
for (mm in 1:length(cand.col)){
st<- sum(obj.mat[mm, ])
}
} else {
frag.i <- matrix(NA, length(cand.col),length(obj.mat.line))
frag.i <- round(obj.mat.ms[1:length(cand.col),],digits=keta)
for (m in 1:length(cand.col) ){
if ( eval(parse(text=paste(MS.ext,sep=""))) ){ 
deconv.ratio <- eval(parse(text=paste("min(",Sum.MS.ratio.ext,")",sep="") ))
st.sum <- max(as.matrix(deconv.ratio*database[2*i,Sum.MS.seq]))
st <- eval(parse(text=area.get))
} else{
st <- 0
}
oneline[m,] <- st
}
}


max.cand.box <- matrix(0, 2*RowRange+1,(icol.range))
max.cand.num.box <- matrix(0, 2*RowRange+1,(icol.range))

noise.sample.num <- 0
maxbox2remb <- rbind(0,0)
r.range <- 0


oneline.mat <- oneline[ (row*headbuffer+1):(length(oneline)-row*tailbuffer) ]
for( oneline.r in 0:( min(icol.range, round(2*ColRange)) ) ){
max.cand <- oneline.mat[((oneline.r*row)+targetrow.ini.set[i]):((oneline.r*row)+targetrow.ini.set[i]+2*RowRange)]
max.cand.num <- ((oneline.r*row)+targetrow.ini.set[i]):((oneline.r*row)+targetrow.ini.set[i]+2*RowRange)
max.cand.box[,oneline.r + 1] <- max.cand
max.cand.num.box[,oneline.r + 1] <- max.cand.num
}

max.cand.box <- replace(max.cand.box, is.na(max.cand.box), 0)
max.cand.box.row <- ceiling(which(max.cand.box==max(max.cand.box))[1]/nrow(max.cand.box))
peak.cand.range <- as.vector(max.cand.num.box)
peak.top.row <- eval(parse(text= rownames(as.matrix(oneline[ round(median(which(oneline==max(  replace(oneline[  peak.cand.range  ],is.na(oneline[  peak.cand.range  ]),0)  )))) ,]))   ))

if ( max.cand.box.row<=1 ){ 
tempnumber <- 1 
} else {
for (qq in 2:max.cand.box.row ) { 
if ( (row*qq) > (peak.top.row%%row-serch.range) ) {
tempnumber <- 1 
break 
} else {
if ( max( na.omit(oneline.mat[(-(row*qq)+peak.top.row%%row-serch.range ):(-(row*qq)+peak.top.row%%row+serch.range )])) >= max( na.omit(oneline.mat[(-row*(qq-1)+peak.top.row%%row-serch.range):(-row*(qq-1)+peak.top.row%%row+serch.range)] )) ){ 
start.row <- qq 
break  
} else {
 start.row <- qq-1 
}
}
}
tempnumber <- qq-1 #max.cand.box.row - qq  
}


try(repeat { if( (( r.range>3+tempnumber )&&  ( maxbox2remb[r.range+1] <= maxbox2remb[r.range+2] ) ) || (maxbox2remb[r.range+1]>0 && maxbox2remb[r.range+1] == maxbox2remb[r.range+2]) ) break

peak.top.row.grad <- peak.top.row + row*(r.range-tempnumber)
top.row.range <- (peak.top.row.grad -serch.range):(peak.top.row.grad +serch.range)

if ( eval(parse(text=rownames(oneline)[1] ))>top.row.range[1] ) {
candidate <- oneline[1+serch.range,]
}
if ( eval(parse(text=rownames(oneline) )) <top.row.range[1]+length(top.row.range)  ) {
candidate <- "nun" 
} 
if ( any( (eval(parse(text=rownames(oneline)[1])):eval(parse(text=rownames(oneline)))  )==top.row.range[1]+length(top.row.range) ) & any( (eval(parse(text=rownames(oneline)[1])):eval(parse(text=rownames(oneline)))  )==top.row.range[1] )  ){
candidate <- oneline[ (which(rownames(oneline)==top.row.range[1])):(which(rownames(oneline)==top.row.range[1]+ length(top.row.range) )),]
}
if ( class(candidate) == "numeric" ){
first.top <- eval(parse(text=rownames( as.matrix(which(candidate==max(candidate))) )))
} else {
first.top <- 0
}

first.range <- (first.top-serch.range):(first.top+serch.range)
if (  first.range[1] < eval( parse(text=rownames(oneline)[1] ))  ){
first.range <- eval( parse(text=rownames(oneline)[1] )):(eval( parse(text=rownames(oneline)[1] ))+2*serch.range) 
}
if (  tail(first.range,1) > eval( parse(text=rownames(oneline) ))  ){
first.range <- (eval( parse(text=rownames(oneline) ))-2*serch.range):eval( parse(text=rownames(oneline) )) 
}


candidate2 <- oneline[ (which(rownames(oneline)==first.range[1])):(which(rownames(oneline)==first.range[1])+length(first.range)-1),]

max.rowline <- (first.top-round(serch.range2/2)):(first.top+serch.range2)
if (  max.rowline[1]< eval( parse(text=rownames(oneline)[1] ))  ){
max.rowline <- round((eval( parse(text=rownames(oneline)[1] ))):round((eval( parse(text=rownames(oneline)[1] ))+serch.range2*1.5)))
}
if (  tail(max.rowline,1) > eval( parse(text=rownames(oneline) ))  ){
 max.rowline <- round((eval( parse(text=rownames(oneline) ))-serch.range2*1.5)):round(eval( parse(text=rownames(oneline) ))) 
}

if((which(rownames(oneline)==max.rowline[1])+(length(max.rowline))-1) > length(oneline)){
max.box2 <- as.matrix(oneline[ (which(rownames(oneline)==max.rowline[1])):length(oneline),])
original.max.box2[ (which(rownames(original.max.box2)==max.rowline[1])):(which(rownames(original.max.box2)==max.rowline[1])+length(max.box2)-1), ] <- max.box2
}else{
max.box2 <- as.matrix(oneline[ (which(rownames(oneline)==max.rowline[1])):(which(rownames(oneline)==max.rowline[1])+(length(max.rowline))-1),])
original.max.box2[ (which(rownames(original.max.box2)==max.rowline[1])):(which(rownames(original.max.box2)==max.rowline[1])+(length(max.rowline))-1), ] <- max.box2
}


if ( class(candidate) == "numeric" ){
order <- as.matrix(which(max.box2 <= max(max.box2) ))
max.id <- which(rownames(max.box2)==first.top )[1]
n <- length(order)
try(for (t in 1:n) {
if (order[t]< max.id) {
if (t < n){
if (max.box2[t,]>=max.box2[t+1,]){
max.box2[1:(t),] <- 0 
}
}

} else {
if (t < n){
if (max.box2[t,]<=max.box2[t+1,]){
max.box2[(t+1):n,] <- 0
}
} 
}
})
} else { 
max.box2[1:length(max.box2)] <- 0
}

maxbox2remb[r.range+3] <- max(max.box2) #sum(max.box2)
r.range <- r.range+1
noise.sample.num <- noise.sample.num + length(max.box2) - length(which(max.box2>0))
oneline.box[eval(parse(text=rownames(max.box2))[1]):(eval(parse(text=rownames(max.box2))[1])+length(max.box2)-1),] <- max.box2
})

oneline.box[eval(parse(text=rownames(max.box2))[1]):(eval(parse(text=rownames(max.box2))[1])+length(max.box2)-1),] <- oneline.box[eval(parse(text=rownames(max.box2))[1]):(eval(parse(text=rownames(max.box2))[1])+length(max.box2)-1),]-oneline.box[eval(parse(text=rownames(max.box2))[1]):(eval(parse(text=rownames(max.box2))[1])+length(max.box2)-1),]
res.box0 <- res.box0 + oneline.box


if ( database[2*i-1,4]==1 ){
if ( database[2*i-1,4]==1){
if ( sum(oneline.box) <= 0  ){
database.firstIS <- database.origin
} 
database.pre <- database.firstIS
ColRange.firstIS <- third.colrange
RowRange.firstIS <- third.colrange
addsearch.firstIS <- third.addserch
}

if ( sum(oneline.box) > 0 ){
GC1 <- database[2*i-1,6]
GC2 <- database[2*i-1,7]
RTcorr.cell <- which(oneline.box==max(oneline.box))[1]
RTcorr.col <- ceiling( RTcorr.cell /row )
RTcorr.row <- RTcorr.cell - ( row* (RTcorr.col-1)  )
GC1corr <- RTini  + (RTcorr.col-1)*MPeriod/60
if (GC1corr!=RTini){
GC2corr <- RTcorr.row/SamRate 
GC1.add <- GC1corr-GC1
GC2.add <- GC2corr-GC2
database[,6] <- database[,6]+GC1.add  
database[,7] <- database[,7]+GC2.add 
if (i < length(database[,1])/2){
if( any(database[seq(1,length(database[,1]),2),7] >MPeriod) ){
database[which(database[,7]>MPeriod),7] <- database[which(database[,7]>MPeriod),7] -MPeriod 
}
if( any(database[seq(1,length(database[,1]),2),7] <0) ){
database[which(database[,7]<0),7] <- database[which(database[,7]<0),7] + MPeriod 
}
}
}
ColRange <- third.colrange  
RowRange <- third.rowrange
addsearch <- third.addserch  
targetcol.ini.set <-0
targetrow.ini.set <-0

for (h in 1:(length(database[,1])/2) ){
targetcol.ini <- as.matrix(round((database[2*h-1,6]-RTini)*60/MPeriod)- round(ColRange) )
if( targetcol.ini <0 ){ 
targetcol.ini <- replace(targetcol.ini,targetcol.ini <1,1) 
} 
if ((round((database[2*h-1,7])/MPeriod*row)-RowRange)>0){
targetrow.ini <- as.matrix(round((database[2*h-1,7])/MPeriod*row)-RowRange)
if(targetrow.ini > row){
targetrow.ini <- as.matrix(row)
}
if(targetrow.ini < 1){targetrow.ini <- as.matrix(1)}
} else{
targetrow.ini <- as.matrix(1)
}
rownames(targetcol.ini) <- database[2*h-1,2]
rownames(targetrow.ini) <- database[2*h-1,2]
targetcol.ini.set <- rbind(targetcol.ini.set,targetcol.ini)
targetrow.ini.set <- rbind(targetrow.ini.set,targetrow.ini)
}

targetcol.ini.set <- as.matrix(targetcol.ini.set[-c(1),])
targetrow.ini.set <- as.matrix(targetrow.ini.set[-c(1),])

icol.range <- round( ColRange + addsearch + 1  + headbuffer + tailbuffer )

if ( database[2*i-1,4]==1){
database.firstIS <- database
database.pre <-database
ColRange.firstIS <- third.colrange
RowRange.firstIS <- third.colrange
addsearch.firstIS <- third.addserch
}


cand.col <- matrix(0, row,icol.range  )
}}


noise.chromat <- original.max.box2 - oneline.box
noise.sample <- length(which(noise.chromat>0))
noise.totalsample <- length(noise.chromat)


Area.value <- sum(oneline.box)
noise.chromat <- as.numeric(noise.chromat)
noise <- 2*3*sd(noise.chromat)
signal <- max(oneline.box) - mean(noise.chromat)
sn <- signal/noise


RT.cell <- which(oneline.box==max(oneline.box))[1]
RT.col <- ceiling( RT.cell /row )
RT.row <- RT.cell - ( row* (RT.col-1)  )

GC1time <- round(RTini  + (RT.col-1)*MPeriod/60,digits=2)
GC2time <- round(RT.row/SamRate,digits=2)
if (GC2time < 0) {
GC2time <- GC2time+MPeriod
}
if (round(GC1time,digits=1)==round(RTini,digits=1)){
GC1time <- 0
GC2time <- 0
}

Area.res <- cbind(paste(database[2*i-1,2]),Area.value,sn,noise.sample,noise.sample.num,GC1time,GC2time, names(oneline.box[which(oneline.box==max(oneline.box)),1][1]), paste(database[2*i-1,5]), paste(database[2*i-1,8]) )
result.box[i,] <- Area.res


redidu <- row-SamRate*MPeriod
redidu.accum <-0
blank.id <- round(phase.shift*SamRate+1)
cont.num<-0
res.box <- as.matrix( eval(parse(text=paste("res.box",cont.num,sep=""))) )
mat <- matrix(0,row,col)


try(for (u in round(phase.shift*SamRate+1):(length(res.box[,1])+round(phase.shift*SamRate)) ){
if (u%%row==0){
redidu.accum <-redidu.accum + redidu
}
if ( redidu.accum > 1  ){
blank.id <- blank.id +1
mat[row*ceiling(blank.id/row)-(blank.id-1),ceiling( (u-round(phase.shift*SamRate))/row )] <- res.box[round(u-phase.shift*SamRate),]
blank.id <- blank.id +1
mat[row*ceiling(blank.id/row)-(blank.id-1),ceiling( (u-round(phase.shift*SamRate))/row )] <- res.box[round(u-phase.shift*SamRate),]
redidu.accum <- redidu.accum -1 
} else {
blank.id <- blank.id +1
mat[row*ceiling(blank.id/row)-(blank.id-1),ceiling( (u-round(phase.shift*SamRate))/row )] <- res.box[round(u-phase.shift*SamRate),]
}
})

zz1 <- t(apply(as.matrix(mat[1:nrow(mat),1:ncol(mat)]),2,rev))
xx1 <- (1:nrow(zz1))/col*RTruntime+RTini
yy1 <-(((1:ncol(zz1))/row*MPeriod))
max.zz<-max(zz1)


dev.new(width=8,height=4)
par(mar=c(5,7,1,1))
if (i>1){
dev.off()
}
filled.contour(xx1,yy1,zz1,xlab="GC1 (min)",ylab="GC2 (sec)",col=terrain.colors(30),levels=seq(0,max(res.box/10)+0.1,max(res.box/10)/30+0.1) )
text(20,3.8,paste("Combined EIC","Current;",result.box[i,1],sep=","))


}) 

dev.copy(jpeg,paste("CombinedEIC.jpeg",sep=""),width=960,height=480)
dev.off()
dev.off()

ISref.result.box <- matrix(0,length(result.box[,1]),length(result.box[1,])+2 )
for (i in 1:length(result.box[,1]) ){
if ( eval(parse(text=paste(result.box[i,9])))==0 ){
ISref.result.box[i,1:length(result.box[1,])] <- result.box[i,]
ISref.result.box[i,length(result.box[1,])+1 ] <- eval(parse(text=paste(result.box[i,2]))) / eval(parse(text=paste(result.box[ which(database[seq(1,length(database[,1]),2),1]==database[(i*2-1),8]) ,2 ])))
}
}
ISref.result.box[ceiling(which(database[,5]==1)/2),1:length(result.box[1,])] <- result.box[ceiling(which(database[,5]==1)/2),]
colnames(ISref.result.box) <- c(colnames(result.box),"ISratio","ISratio_fitting") 
write.csv( ISref.result.box,paste("TSEN_result.csv",sep="") )

}
