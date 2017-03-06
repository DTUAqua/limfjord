## Script input:
###############
WEIGHT <- c("total", "commercial")[1]

###############

library(DATRAS)
d <- readICES("../data/Exchange_Data_all_years_v1.csv")

## Add 2016 data
d2 <- readICES("../data/Exchange_Data_2016.csv")
d <- c(d,d2)

## Remove invalid hauls (I/V)
d <- subset(d, HaulVal == "V")

d <- subset(d,lon<9.5 & lat>56)
d$sweptArea <- d$Distance * 1.20 ## Swept area in m^2
d <- subset(d, Year %in% 2004:2016)

## *Total* weight gram
w <- tapply(d[["HL"]]$CatCatchWgt,d[["HL"]]$haul.id,function(x)x[1])
w[is.na(w)] <- 0
d[["HH"]]$TotalWeight <- w[as.character(d$haul.id)]

## Select commercial size interval
## Fit length-weight relation and convert to 'target-catch'
d <- addSpectrum(d,by=.5)
spectrum2totalweight <- function(a,b){
    N <- d$N
    size <- .5 * (1:ncol(N))
    rowSums(  a * (N %*% (size^b) )  )
}
fit <- nls(TotalWeight ~ spectrum2totalweight(a,b), start=c(a=1,b=1)  , data=d[[2]])
## Got: a=0.2186047 b=2.44266826
## How good is it:
layout(t(1:2))
plot(rowSums(d$N),d$TotalWeight, main="Weigth ~ Total count")
plot(predict(fit),d$TotalWeight, main="Weight ~ sum a*count(l)^b")
## Convert to *commercial* weight
if(WEIGHT == "commercial"){
    above <- "[4.5,5)" ## From 4.5 cm and above
    keep <- as.logical(cumsum(colnames(d$N) == above))
    d[["HH"]]$N[,!keep] <- 0
    d[["HH"]]$Weight <- spectrum2totalweight(coef(fit)["a"],coef(fit)["b"])
} else if(WEIGHT == "total"){
    d[["HH"]]$Weight <- d[["HH"]]$TotalWeight
} else {
    stop()
}

## plot(d)
## d <- addSpectrum(d)
source("gridConstruct.R") ## Create objects: gr, map.pol

## Attach grid cell to each haulid:
d[[2]]$gf <- gridFactor(d,gr)

## Add spatial regions for aggregation:
## plot(as.data.frame(gr),col=unclass(spatialRegions))

## Avoid NAs in spatialRegions:
levels(spatialRegions) <- c(levels(spatialRegions), "NA")
spatialRegions[is.na(spatialRegions)] <- "NA"

## Get area of grid cells (km^2):
gridCellArea <- mean(summary(as.polygons(gr))$side.length)^2

## Change invalid TMB name:
d[[2]]$haulId <- d[[2]]$haul.id

## Response variables
response <- d$Weight
presence <- as.numeric(response!=0)
d[[2]]$response <- response
d[[2]]$presence <- presence

## Factor defining our time unit
time <- factor(d$Year)
d[[2]]$time <- time

## Sparse matrices for GMRF: Q = Q0+delta*I
Q0 <- -attr(gr,"pattern")
diag(Q0) <- 0
diag(Q0) <- -rowSums(Q0)
I <- .symDiagonal(nrow(Q0))

## Load commercial data
if(FALSE) {
load("/nobackup/kaskr/commercial_data/d_commercial.RData")
lon.tol <- 0.08*2 ## Lon = 0.08  <===> ca 5 km = 10 grid cells ! :
dist.km(data.frame(lon=8.5, lat=56.5), data.frame(lon=8.5 + lon.tol, lat=56.5))
lat.tol <- .045*2 ## Lat = 0.045  <===> ca 5 km = 10 grid cells ! :
dist.km(data.frame(lon=8.5, lat=56.5), data.frame(lon=8.5, lat=56.5 + lat.tol))
## Take subset:
dim(d_commercial)
d_commercial <- subset(d_commercial,
                       is.finite(lon) &
                       lon.range < lon.tol &
                       lat.range < lat.tol &
                       fangst >= 0 )
dim(d_commercial)
d_commercial$gf <- gridFactor(d_commercial,gr)
d_commercial$Year <- format(d_commercial$date, "%Y")
d_commercial$month <- as.numeric(format(d_commercial$date, "%m"))
dist_to_nearest_grid <- dist.km(d_commercial[c("lon","lat")], gr[d_commercial$gf,], outer=FALSE)
d_commercial <- subset(d_commercial, dist_to_nearest_grid < 1)
}

library(TMB)
##compile("mussel.cpp","-O0 -g")
compile("mussel.cpp")
dyn.load(dynlib("mussel"))

##fitModel <- function(){
  ## DATRAS package 'subset' removes empty levels from gridFactor object (!)
  ## So need to re-create...
  ##d[[2]]$gf <- gridFactor(d,gr)
data <- list(
             Q0=Q0  ,
             I=I  ,
             time=d[[2]]$time ,
             gf=d[[2]]$gf  ,
             haulId=d[[2]]$haulId  ,
             response=d[[2]]$response  ,
             presence=d[[2]]$presence,
    spatialRegions = spatialRegions,
    sweptArea = d[[2]]$sweptArea,
    gridCellArea = gridCellArea,
    reportLog = 1
             )
parameters <- list(
                   eta_presence = matrix(0,nrow(Q0),nlevels(time)) ,
                   eta_density = matrix(0,nrow(Q0),nlevels(time)) ,
                   logdelta=c(0,0),
                   logscale=matrix(0,nlevels(time),2),
                   logsd_nugget=rep(0,nlevels(time)),
                   mu=rep(0,nlevels(time))
                   )
map <- parameters[-(1:2)]
map <- lapply(map,factor) ## Collect all !!!
obj <- MakeADFun(data=data,
                 parameters=parameters,
                 random="eta",
                 regexp=TRUE,
                 DLL="mussel",
                 map=map
                 )
obj$par <- structure(c(-6.52786685094663, 0.511183138204372, 0.60252204243289,  8.02634171785195), .Names = c("logdelta", "logscale", "logsd_nugget",  "mu"))
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))

pl <- obj$env$parList()

###########################################
## Use as initial guess for next run:
map3 <- list()
map3$logscale <- factor(col(parameters$logscale))
map3$logsd_nugget<- factor(parameters$logsd_nugget* 0)
obj3 <- MakeADFun(data=data,
                  parameters=pl,
                  random=c("^eta"),
                  profile="mu",
                  regexp=TRUE,
                  DLL="mussel",
                  map=map3
                  )
obj3$par <- structure(c(-7.5724043055008, -4.46731216093986, 0.207982042852333,  0.907872821912759, 0.278968517094824), .Names = c("logdelta",  "logdelta", "logscale", "logscale", "logsd_nugget"))
system.time(opt3 <- nlminb(obj3$par,obj3$fn,obj3$gr))
hessian <- optimHess(opt3$par, obj3$fn, obj3$gr)
eigen(hessian)$val

system.time( sdrep <- sdreport(obj3, hessian = hessian, bias.correct=TRUE,
                               getReportCovariance=FALSE,
                               bias.correct.control=list(sd=FALSE, nsplit=10) ) ) ## <--- very memory intensive when reportLog = TRUE
pl <- obj3$env$parList(par=obj3$env$last.par.best)
rep <- obj3$report(obj3$env$last.par.best)
obj3$env$data$reportLog <- 0 ## Natural scale report
system.time( sdrep0 <- sdreport(obj3, hessian = hessian, bias.correct=TRUE) )

## Update with Natura 2000 areas:
stop()
data$spatialRegions <- spatialRegions2
levels(spatialRegions2) <- c(levels(spatialRegions2), "NA")
spatialRegions2[is.na(spatialRegions2)] <- "NA"
obj3 <- MakeADFun(data=data,
                  parameters=pl,
                  random=c("^eta"),
                  profile="mu",
                  regexp=TRUE,
                  DLL="mussel",
                  map=map3
                  )
obj3$fn(opt3$par)
obj3$env$data$reportLog <- 0 ## Natural scale report
system.time( sdrep00 <- sdreport(obj3, hessian = hessian, bias.correct=TRUE) )
tab <- as.data.frame(summary(sdrep00, "report")); row.names(tab) <- NULL
tab$Year <- rep(levels(data$time),each=3)
tab$Area <- levels(spatialRegions2)
tab <- tab[tab$Area != "NA", ]
tab <- tab[c(5,6,3,2)]
names(tab)[3] <- "Estimate"

rm(obj,obj3)
save.image(file="rapport2016.RData")
