library(DATRAS)

#############################################################
## 1. Read the builtin data. Alternatively use 'readICES' or
##   'readExchange' to read DATRASraw object from a file.
#############################################################

data(mussel) ## Get data object 'd'

#############################################################
## 2. Data pre-processing
#############################################################

## Remove invalid hauls (I/V)
d <- subset(d, HaulVal == "V")

## Space time subset:
d <- subset(d,lon<9.5 & lat>56)
d <- subset(d, Year %in% 2004:2016)

## Add swept area in m^2 to dataset:
d$sweptArea <- d$Distance * 1.20

## Add response variable: Total weight (g)
w <- tapply(d[["HL"]]$CatCatchWgt,d[["HL"]]$haul.id,function(x)x[1])
w[is.na(w)] <- 0
d[["HH"]]$TotalWeight <- w[as.character(d$haul.id)]

#############################################################
## 3. Add spatial grid
#############################################################
