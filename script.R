## Script for stockassessment.org
.libPaths(".")
devtools:::install_github("kaskr/gridConstruct/gridConstruct")
devtools:::install_github("DTUAqua/limfjord/mussel")

## Load packages and data
library(DATRAS)
library(mussel)
d <- readExchange("data.zip")

## Remove invalid hauls (I/V)
d <- subset(d, HaulVal == "V")

## Space time subset:
d <- subset(d,lon<9.5 & lat>56)
latestYear <- max(levels(d$Year))
d <- subset(d, Year %in% 2004:as.numeric(latestYear))

## Add swept area in m^2 to dataset:
Width <- c("MSK" = 1.20, "MSKOES" = 1.00) ## FIXME: New gear 1.00 width ?
d$sweptArea <- d$Distance * Width[as.character(d$Gear)]

## Add response variable: Total weight (g)
w <- tapply(d[["HL"]]$CatCatchWgt,d[["HL"]]$haul.id,function(x)x[1])
w[is.na(w)] <- 0
d[["HH"]]$Weight <- w[as.character(d$haul.id)]

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
