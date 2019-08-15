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
above <- "[4.5,5)" ## From 4.5 cm and above
remove <- as.logical(cumsum(colnames(d$N) == above))
d[["HH"]]$N[,remove] <- 0
if (TRUE) {
    ## Add response variable: Total weight (g) of selected size groups
    d[["HH"]]$Weight <- spectrum2totalweight(coef(fit)["a"],coef(fit)["b"])
} else {
    ## OR, use the total weight
    d[["HH"]]$Weight <- d[["HH"]]$TotalWeight
}

## Remove invalid hauls (I/V)
d <- subset(d, HaulVal == "V")

## Space time subset:
d <- subset(d,lon<9.5 & lat>56)
latestYear <- max(levels(d$Year))
d <- subset(d, Year %in% 2004:as.numeric(latestYear))

## Add swept area in m^2 to dataset:
Width <- c("MSK" = 1.20, "MSKOES" = 1.00) ## FIXME: New gear 1.00 width ?
d$sweptArea <- d$Distance * Width[as.character(d$Gear)]

## NO longer done here:
## Add response variable: Total weight (g)
## w <- tapply(d[["HL"]]$CatCatchWgt,d[["HL"]]$haul.id,function(x)x[1])
## w[is.na(w)] <- 0
## d[["HH"]]$Weight <- w[as.character(d$haul.id)]

## Standardize response to kg/m^2
d[["HH"]]$Weight <- d[["HH"]]$Weight / d[["HH"]]$sweptArea
d$sweptArea[] <- 1

#############################################################
## 3. Add spatial grid
#############################################################

data(grid, package="mussel") ## Get default grid

#############################################################
## 4. Fit model
#############################################################

env <- fitModel(d)

#############################################################
## 5. Save results
#############################################################

save.image("res.RData.exe") ## .exe == don't check in :)
