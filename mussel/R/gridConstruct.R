if(FALSE){
## Script input: DATRASraw object d

## Grid
library(gridConstruct)
library(rgdal)

## Objects to construct
gr <- NULL
map.pol <- NULL
spatialRegions <- NULL

gridConstruct <- function(d){
    ## gr <- gridConstruct(d,km=.5)
    ## points(d$lon,d$lat)
    ## map("worldHires",add=TRUE)
    gr <- gridConstruct::gridConstruct(d,km=.5,filter=FALSE)
    ## Read shape file data
    ## Depth data:
    ##  shape <- readOGR("../shpfiles/","Bathy_Lim_1m_ploy")
    file <- system.file("shp/Denmark",package="mussel")
    shape <- readOGR(file,"Kystlinie")
    proj4 <- proj4string(shape)
    map <- spTransform(shape,CRS("+proj=longlat"))
    gr2 <- as.data.frame(gr)
    coordinates(gr2) <- ~lon + lat
    proj4string(gr2) <- CRS("+proj=longlat")
    xtra <- over(gr2, map)
    gr3 <- gr[is.na(xtra$ID),]
    cc <- connectedComponents(gr3)
    gr3 <- gr3[cc[[which.max(sapply(cc,length))]],]
    while(any(rowSums(attr(gr3,"pattern")) <=2 ) ){
        gr3 <- gr3[rowSums(attr(gr3,"pattern"))>2,]
    }
    ##plot(map,add=TRUE)
    f <- function(i){
        x <- map@polygons[[i]]
        do.call("rbind",lapply(x@Polygons,function(x)rbind(x@coords,NA)))
    }
    qw <- do.call("rbind",lapply(1:length(map@polygons),f))
    plot(gr3)
    points(qw[,1],qw[,2],type="l")
    polygon(qw[,1],qw[,2],col="grey")
    ##gr2 <- spTransform(gr2,CRS(proj4))
    ##xtra <- over(gr2, shape)
    ##image(gr,xtra$Lower_,map=NULL)
    gr <<- gr3
    map.pol <<- qw
    ## Additional regions
    file <- system.file("shp/Regions",package="mussel")
    shape <- readOGR(file,"Limfjord_omraader_nov_2005_area")
    proj4 <- proj4string(shape)
    regions <- spTransform(shape,CRS("+proj=longlat"))
    gr2 <- as.data.frame(gr3)
    coordinates(gr2) <- ~lon + lat
    proj4string(gr2) <- CRS("+proj=longlat")
    xtra <- over(gr2, regions)
    fac <- factor(xtra[[2]])
    remap <- c("Lovns Bredning, Øst", "Lovns Bredning, Vest")
    levels(fac)[levels(fac) %in% remap] <- "Lovns Bredning"
    spatialRegions <<- fac
    NULL
}


plotMap <- function(..., add=FALSE){
    if(!add)plot(gr,...)
    points(map.pol[,1],map.pol[,2],type="l")
    polygon(map.pol[,1],map.pol[,2],col="grey")
}

gridConstruct(d)


## New spatial region that approximate natura 2000 areas:
spatialRegions2 <- spatialRegions
area1 <- c("Bjørnsholm Bugt",
           "Løgstør Bredning",
           "Løgstør Bredning, Øst",
           "Løgstør Bredning, Vest",
           "Løgstør Grunde",
           "Livø Bredning, Øst")
area2 <- c("Lovns Bredning")
levels(spatialRegions2)[!(levels(spatialRegions2) %in% c(area1,area2))] <- NA
levels(spatialRegions2)[levels(spatialRegions2) %in% area1] <- "Natura Loegstoer"
levels(spatialRegions2)[levels(spatialRegions2) %in% area2] <- "Lovns Bredning"
}
