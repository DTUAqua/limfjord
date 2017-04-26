fitModel <- function(d, getLog=FALSE) {

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
diag(Q0) <- -Matrix::rowSums(Q0)
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
                  parameters=parameters,
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

if(getLog){
    system.time( sdrep <- sdreport(obj3, hessian = hessian, bias.correct=TRUE,
                                   getReportCovariance=FALSE,
                                   bias.correct.control=list(sd=FALSE, nsplit=10) ) ) ## <--- very memory intensive when reportLog = TRUE
}
pl <- obj3$env$parList(par=obj3$env$last.par.best)
rep <- obj3$report(obj3$env$last.par.best)
obj3$env$data$reportLog <- 0 ## Natural scale report
system.time( sdrep0 <- sdreport(obj3, hessian = hessian, bias.correct=TRUE) )

environment()
}


#############################################################################



at10 <- 10^(-10:10)
label10 <- paste0(10,"^",-10:10)
label10 <- as.expression(lapply(label10,function(x)parse(text=x)[[1]]))
mypretty <- function(br, tol=1){
    rg <- range(log10(exp(br)))
    rg0 <- c(floor(rg[1]),ceiling(rg[2]))
    rg <- ifelse(abs(rg-rg0)<tol, rg0, rg)
    seq(log(10^rg[1]),
        log(10^rg[2]),
        length=length(br))
}

plotMap <- function(env,time="2016",...){
    env$showTime <- time
    local({
        layout(matrix(1:4,2))
##################### Plot results
        y <- which(levels(data$time)==showTime)
        presence.prob <- 1/(1+exp(-pl$eta_presence[,y]))
        plot(gr,type="n",las=1,xlab="Longitude",ylab="Latitude",main="Presence probability")
        br <- seq(0,1,length=17)
        image(gr,presence.prob, map=quote(plotMap(add=TRUE)), breaks=br, add=TRUE)
        par(usr=c(0,1,min(br),max(br)));axis(4, las=1, hadj=.3)
########################
        ## Map of density given presence
        plot(gr,type="n",las=1,xlab="Longitude",ylab="Latitude", main="Density given presence")
        resp <- pl$eta_density[,y]+pl$mu[y]    + log(1/(1000 * mean(d$sweptArea))) ## UNIT
        br <- seq(min(resp),max(resp), length=17)
        br <- mypretty(br)
        image(gr, resp , map=quote(plotMap(add=TRUE)),  add=TRUE, breaks=br)
        par(usr=c(0,1,min(br),max(br)));axis(4, at=log(at10), label=label10, las=1, hadj=.3)
        box()
########################
        ## Map of density
        log_density_survey <- pl$eta_density[,y] + log(1/(1+exp(-pl$eta_presence[,y]))) + pl$mu[y]
        resp <- log_density_survey    + log(1/(1000 * mean(d$sweptArea))) ## UNIT
        plot(gr,type="n",las=1,xlab="Longitude",ylab="Latitude", main="Density")
        br <- seq(min(resp),max(resp), length=17)
        br <- mypretty(br)
        image(gr, resp , map=quote(plotMap(add=TRUE)),  add=TRUE, breaks=br)
        par(usr=c(0,1,min(br),max(br)));axis(4, at=log(at10), label=label10, las=1, hadj=.3)
        box()
########################
## Map of biomass-density
## b = 2.703 * ((C/A)/1000)^0.29
        logC <- pl$eta_density[,y] + log(1/(1+exp(-pl$eta_presence[,y]))) + pl$mu[y]
        A <- mean(d$sweptArea)
        logb <- log(  2.703 * ((exp(logC)/A)/1000)^0.29  )
        resp <- logb
        plot(gr,type="n",las=1,xlab="Longitude",ylab="Latitude", main="Biomass-density")
        br <- seq(min(resp),max(resp), length=17)
        br <- mypretty(br, .5)
        image(gr, resp , map=quote(plotMap(add=TRUE)),  add=TRUE, breaks=br)
        par(usr=c(0,1,min(br),max(br)));axis(4, at=log(at10), label=label10, las=1, hadj=.3)
        box()
        ## End
        title(paste("Time:",showTime), outer=TRUE, line=-1)
    }, env)
}

########################
## Map of biomass-density
## b = 2.703 * ((C/A)/1000)^0.29
plotTimeSeriesLog <- function(env,selectRegions = c("Thisted Bredning, Sydvest", "Lovns Bredning", "Kås Bredning, Vest", "Skive Fjord"),...) {
    env$selectRegions <- selectRegions
    local({
        est <- sdrep$unbiased$value
        sd <- summary(sdrep,"report")[,2]
        ind <- (rbind(matrix(1:length(est),nlevels(spatialRegions),nlevels(time))))
        mat <- cbind(est,est-1.96*sd,est+1.96*sd)
        k <- match(selectRegions, levels(spatialRegions))
        for(i in 1:length(selectRegions)){
            ylim <- range(exp(mat[ind[k,],]))
            xlab <- "Year"
            ylab <- expression(`b  `(kg/m^2))
            matplot(as.numeric(levels(time)),exp(mat[ind[k[i],],]),type=c("b","l","l"),axes=FALSE, col=1, lty=c(1,2,2), pch=1, xlab="", ylab="", ylim=ylim)
            mtext(ylab,2,line=2.5)
            mtext(xlab,1,line=2.5)
            axis(1)
            axis(2, las=1)
            box()
            legend("topright", selectRegions[i])
        }
    }, env)
}

plotTimeSeries <- function(env, selectRegion = c("Lovns Bredning"),...) {
    env$selectRegions <- selectRegion
    local({
        est <- sdrep0$unbiased$value
        sd <- summary(sdrep0,"report")[,2]
        ind <- (rbind(matrix(1:length(est),nlevels(spatialRegions),nlevels(time))))
        mat <- cbind(est,est-1.96*sd,est+1.96*sd)
        k <- match(selectRegions, levels(spatialRegions))
        newmat <- mat[ind[k,],]
        rownames(newmat) <- levels(time)
        newmat <- newmat[-1,]
        matplot(as.numeric(rownames(newmat)),newmat,type="l", col=1, lty=c(1,2,2), xlab="", ylab="", las=1, ylim=c(0,max(newmat)))
        xlab <- "Year"
        ylab <- expression(`b  `(kg/m^2))
        mtext(ylab,2,line=2.5)
        mtext(xlab,1,line=2.5)
        ##points(as.numeric(names(b)),b)
    }, env)
}

