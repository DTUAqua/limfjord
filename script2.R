## Script for stockassessment.org
.libPaths(".")
devtools:::install_github("kaskr/gridConstruct/gridConstruct")
devtools:::install_github("DTUAqua/limfjord/mussel")

## Load packages and data
library(DATRAS)
library(mussel)
library(sp)
library(rgdal)
##d <- readExchange("data.zip")
load("blackbox/confidential/d_lag.RData")

## !!!! WARNING !!!! : May have to do this in original script as well !!!
## NA Species are mussels !!!
d <- subset(d, is.na(Species))

## Fit length-weight relation
d <- addSpectrum(d, by=.5)
w <- tapply(d[["HL"]]$CatCatchWgt,d[["HL"]]$haul.id,function(x)x[1])
w[is.na(w)] <- 0
d[["HH"]]$TotalWeight <- w[as.character(d$haul.id)]
spectrum2totalweight <- function(a,b,N){
    size <- .5 * (1:ncol(N))
    rowSums(  a * (N %*% (size^b) )  )
}
N <- d$N
fit <- nls(TotalWeight ~ spectrum2totalweight(a,b,N), start=c(a=1,b=1)  , data=d[[2]])
## Got: a=0.2186047 b=2.44266826
## How good is it:
layout(t(1:2))
plot(rowSums(d$N),d$TotalWeight, main="Weigth ~ Total count")
plot(predict(fit),d$TotalWeight, main="Weight ~ sum a*count(l)^b")
above <- "[4.5,5)" ## From 4.5 cm and above
remove <- as.logical(cumsum(colnames(d$N) == above)) ## Remove large individuals

## Add swept area in m^2 to dataset:
Width <- c("MSK" = 1.20, "MSKOES" = 1.00) ## FIXME: New gear 1.00 width ?
d$sweptArea <- d$Distance * Width[as.character(d$Gear)]

if (TRUE) {
    ## SMALL INDIVIDUALS
    ## Add response variable: Total weight (g) of selected size groups
    N <- d$N
    N[,remove] <- 0
    d[["HH"]]$WeightSmall <- spectrum2totalweight(coef(fit)["a"],coef(fit)["b"],N) /
        d$sweptArea
    ## LARGE INDIVIDUALS
    N <- d$N
    N[,!remove] <- 0
    d[["HH"]]$WeightLarge <- spectrum2totalweight(coef(fit)["a"],coef(fit)["b"],N) /
        d$sweptArea
    ## ALL INDIVIDUALS
    d[["HH"]]$WeightAll <- d[["HH"]]$TotalWeight /
        d$sweptArea
    ##
    d$sweptArea[] <- 1
}
## Remove invalid hauls (I/V)
d <- subset(d, HaulVal == "V")

## Space time subset:
d <- subset(d,lon<9.5 & lat>56)
latestYear <- max(levels(d$Year))
d <- subset(d, Year %in% 2004:as.numeric(latestYear))

## NO longer done here:
## Add response variable: Total weight (g)
## w <- tapply(d[["HL"]]$CatCatchWgt,d[["HL"]]$haul.id,function(x)x[1])
## w[is.na(w)] <- 0
## d[["HH"]]$Weight <- w[as.character(d$haul.id)]

## Standardize response to kg/m^2
## d[["HH"]]$Weight <- d[["HH"]]$Weight / d[["HH"]]$sweptArea
## d$sweptArea[] <- 1

#############################################################
## 3. Add spatial grid
#############################################################

data(grid, package="mussel") ## Get default grid

#############################################################
## 4. Fit model
#############################################################


## Model without effects
d$Weight <- d$WeightLarge
d[[2]]$X <- matrix(,length(d),0)
env <- fitModel(d,FALSE,FALSE,FALSE)
resid <- log(env$data$response)-env$rep$muvec
i <- is.finite(resid)
plot(d$lag1[i]^(1/4),resid[i])
plot(d$lag2[i],resid[i],log="x")
plot(cut(d$lag1[i],10),resid[i],log="x")
plot(factor(d$lag1[i]!=0),resid[i])


sdr <- sdreport(env$obj3)


d$Weight <- d$WeightLarge
d[[2]]$X <- model.matrix(~poly(lag1,2)*poly(lag2,2),data=d[[2]])
d[[2]]$X <- d[[2]]$X[,-1]
env <- fitModel(d,FALSE,FALSE,FALSE)
sdr <- sdreport(env$obj3,hessian=env$hessian)
## beta          0.08783800 0.08893661
## beta          0.01808333 0.10154617


d[[2]]$X <- log( cbind(d$lag1,d$lag2) + 1)
env <- fitModel(d,FALSE,FALSE,FALSE)
sdr <- sdreport(env$obj3)
## beta          0.50775918 0.32048733
## beta          0.25625001 0.35093103


d[[2]]$X <- cbind(d$lag1>0,d$lag2>0)[,-1]
env <- fitModel(d,FALSE,FALSE,FALSE)
sdr <- sdreport(env$obj3)
## beta          0.86238394  0.2723569 <--- !!!
## beta          0.09693322  0.2975847

d$flag1 <- cut(d$lag1, breaks=c(-1,0,1,100))
d$flag2 <- cut(d$lag2, breaks=c(-1,0,1,100))
d[[2]]$X <- model.matrix(~flag1:flag2-1,data=d[[2]])[,-1]
env <- fitModel(d,FALSE,FALSE,FALSE)
sdr <- sdreport(env$obj3)
## NaN


d$Weight <- d$WeightLarge
d[[2]]$X <- model.matrix(~flag1:flag2-1,data=d[[2]])[,-1]
env <- fitModel(d,FALSE,FALSE,FALSE)
sdrLarge <- sdreport(env$obj3)
## ==== ADULT
## beta          0.88051777 0.14803052
## beta          0.56786924 0.15545536
## beta          1.27682386 0.17117546


##################### FOR REPORT

##################### SIMPLE
d$flag1 <- factor(d$lag1>0)
d$flag2 <- factor(d$lag2>0)

d$Weight <- d$WeightSmall
d[[2]]$X <- model.matrix(~flag1:flag2-1,data=d[[2]])[,-1]
env <- fitModel(d,FALSE,FALSE,FALSE)
sdrSmall <- sdreport(env$obj3,hessian=env$hessian)
d[[2]]$X <- model.matrix(~flag1-1,data=d[[2]])[,-1,drop=FALSE]
env2 <- fitModel(d,FALSE,FALSE,FALSE)
sdrSmall <- sdreport(env$obj3,hessian=env$hessian)
1-pchisq(2*(env2$opt3$objective-env$opt3$objective), df=3)
sdrSmall2 <- sdreport(env2$obj3,hessian=env2$hessian)
## ==== RECRUITMENT
## beta          0.88310599 0.33502158  (lag1 && !lag2)
## beta          0.11997373 0.36829246  (lag2 && !lag1)
## beta          0.94340745 0.36855030  (lag1 &&  lag2)

## Model without effects
d[[2]]$X <- matrix(,length(d),0)
env3 <- fitModel(d,FALSE,FALSE,FALSE)
resid <- log(env3$data$response)-env3$rep$muvec
i <- is.finite(resid)
plot(log(d$lag1[i]),resid[i])
plot(d$lag2[i],resid[i],log="x")
save(resid,file="resid.RData")

#############################################################
## 5. Save results
#############################################################

## d[[2]]$X <- model.matrix(~flag2-1,data=d[[2]])[,-1,drop=FALSE]
## env4 <- fitModel(d,FALSE,FALSE,FALSE)
## sdrSmall4 <- sdreport(env4$obj3,hessian=env4$hessian)


save.image("res.RData.exe") ## .exe == don't check in :)
