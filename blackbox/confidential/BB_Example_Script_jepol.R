#Load packages
library(data.table)
library(sf)
library(mapview)
library(dplyr)

# lst <- list.files("Q:/dfad/users/kaskr/data/18-12-20_Blackbox/Data", pattern = "2016_1|2017_2", full.names = T)
# 
# dat <- rbindlist(lapply(lst, readRDS), fill=T)
#Read data
##dat <- readRDS("Q:/dfad/users/kaskr/data/18-12-20_Blackbox/Data/Polygons_bb_2018.rds")
## dat <- readRDS("Polygons2018.rds")
## dat <- readRDS("Polygons2014.rds")
## dat <- readRDS("Polygons2015.rds")

dat <- lapply(dir(pattern="Polygons2"), readRDS)
Years <- 2012:2018
names(dat) <- Years

plot(table(substring((dat[[7]]$date),6,7)))
points(table(d$Month),col="red")

#Subset data
##dat2 <- dat[dat$date %in% as.Date(c("2018-03-12", "2018-03-13")),]

##Create example points:
##pts <- data.table("Site.ID."=1, y=c(56.759921), x=c(8.845))
##pts <- data.table("Site.ID."=1:2, y=c(56.759921, 56.791358), x=c(8.845, 8.887))

library(DATRAS)
load("../../d.RData")
d <- subset(d, Year %in% 2014:2018) ## Removing year 2012 and 2013

f <- function(i, lag=0) {
    print(i)
    lon <- d$lon[i]
    lat <- d$lat[i]
    k <- which( d$Year[i] == Years )
    pts <- data.table(x=lon, y=lat)
    pts <- pts %>% 
        sf::st_as_sf(coords = c("x","y")) %>% 
        sf::st_set_crs(4326)
    ##Turn into UTM32N
    pts <- st_transform(pts, 32632)
    ##Make buffer 50 m around each points
    buf <- st_buffer(pts, 50)
    ## Which commercial year to consider
    j <- k - lag
    sum(st_area(st_intersection(dat[[j]], buf)))
}

f(2)
library(parallel)
options(mc.cores=4)
##system.time(m <- mclapply(1:4, f))

## lag 0
system.time(lag0 <- mclapply(1:length(d), f, lag=0)) ## 18 minutes
lag0 <- unlist(lag0)
save(lag0, file="lag0.RData")

## lag 1
system.time(lag1 <- mclapply(1:length(d), f, lag=1)) ## 18 minutes
lag1 <- unlist(lag1)
save(lag1, file="lag1.RData")

## lag 2
system.time(lag2 <- mclapply(1:length(d), f, lag=2)) ## 18 minutes
lag2 <- unlist(lag2)
save(lag2, file="lag2.RData")



pts <- data.table("haul.id"=levels(d$haul.id), x=d$lon, y=d$lat)

pts <- pts %>% 
  sf::st_as_sf(coords = c("x","y")) %>% 
  sf::st_set_crs(4326)

#Turn into UTM32N
pts <- st_transform(pts, 32632)

#Make buffer 50 m around each points
buf <- st_buffer(pts, 50)

#Visualize data (doesnt work with large amounts of data)
mapview(dat2) + mapview(pts) + mapview(buf)

#Clip polygons to buffer
int <- st_intersection(dat[[2]], buf)
mapview(int)

#Calculate overlapping area for each site id and each date
ss <-   int %>%
  group_by(date, haul.id) %>% 
  summarise(area_overlap = sum(st_area(geometry)))

#print output
data.table(ss)

