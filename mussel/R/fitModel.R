fitModel <- function(d) {

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

## Data
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

## Parameters
parameters <- list(
    eta_presence = matrix(0,nrow(Q0),nlevels(time)) ,
    eta_density = matrix(0,nrow(Q0),nlevels(time)) ,
    logdelta=c(0,0),
    logscale=matrix(0,nlevels(time),2),
    logsd_nugget=rep(0,nlevels(time)),
    mu=rep(0,nlevels(time))
)


## map <- parameters[-(1:2)]
## map <- lapply(map,factor) ## Collect all !!!
## obj <- MakeADFun(data=data,
##                  parameters=parameters,
##                  random="eta",
##                  regexp=TRUE,
##                  DLL="mussel",
##                  map=map
##                  )
## obj$par <- structure(c(-6.52786685094663, 0.511183138204372, 0.60252204243289,  8.02634171785195), .Names = c("logdelta", "logscale", "logsd_nugget",  "mu"))
## system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
## pl <- obj$env$parList()


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

environment()
}
