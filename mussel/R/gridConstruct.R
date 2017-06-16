## Script input: DATRASraw object d


## Objects to construct
## gr <- NULL
## map.pol <- NULL
## spatialRegions <- NULL
## spatialRegions2 <- NULL

## Optionally override:
gridConstruct <- function(d, km=.5){
    ## Grid
    library(gridConstruct)
    library(rgdal)
    ## gr <- gridConstruct(d,km=km)
    ## points(d$lon,d$lat)
    ## map("worldHires",add=TRUE)
    gr <- gridConstruct::gridConstruct(d,km=km,filter=FALSE)
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
    ## Alle_muslingeomraader_2011_region
    ## Forbudsomraader_musling
    ## Natura2000 = Habitatomraade16_intersect_hav
    ## Habitatomraade30_Intersect_hav_explode_kun_lovns
    ## Limfjord_omraader_nov_2005_area
    lookupRegions <- function(filename = "Limfjord_omraader_nov_2005_area",
                              getShape=FALSE) {
        folder <- system.file("shp/Regions",package="mussel")
        shape <- readOGR(folder, filename)
        proj4 <- proj4string(shape)
        regions <- spTransform(shape,CRS("+proj=longlat"))
        if(getShape)return(regions)
        gr2 <- as.data.frame(gr3)
        coordinates(gr2) <- ~lon + lat
        proj4string(gr2) <- CRS("+proj=longlat")
        xtra <- over(gr2, regions)
        fac <- factor(xtra[[2]])
        fac
    }
    ## Test content
    if(FALSE) {
        files <- sub("\\..*","",dir(folder, pattern="dbf"))
        myf <- function(x)levels(lookupRegions(x))
        out <- lapply(files, myf)
        names(out) <- files
    }

    ## Produktionsomraåder løgstør
    prodomr <- c("Feggesund / Hovsør Havn" = 32,
                 "Bjørnsholm Bugt" = 37,
                 "Løgstør Bredning" = 34,
                 "Løgstør Bredning, Øst" = 38,
                 "Løgstør Bredning, Vest" = 33,
                 "Løgstør Grunde" = 39,
                 "Livø Bredning, Øst" = 36)

    shp_natura2000 <<- lookupRegions("Habitatomraade16_intersect_hav",TRUE)
    shp_lovns <<- lookupRegions("Habitatomraade30_Intersect_hav_explode_kun_lovns",TRUE)
    shp_prod <<- lookupRegions("Alle_muslingeomraader_2011_region",TRUE)

    fac <- lookupRegions()
    remap <- c("Lovns Bredning, Øst", "Lovns Bredning, Vest")
    levels(fac)[levels(fac) %in% remap] <- "Lovns Bredning"
    spatialRegions <<- fac
    NULL
}

lookupShape <- function(gr, shape) {
    proj4 <- proj4string(shape)
    regions <- spTransform(shape,CRS("+proj=longlat"))
    gr2 <- as.data.frame(gr)
    coordinates(gr2) <- ~lon + lat
    proj4string(gr2) <- CRS("+proj=longlat")
    xtra <- over(gr2, regions)
    xtra
}

spatialRegionIndicator <- function(gr) {
    data(shp, package="mussel")
    lovns <- as.numeric(!is.na(lookupShape(gr, shp_lovns)[[1]]))
    natura2000 <- as.numeric(!is.na(lookupShape(gr, shp_natura2000)[[1]]))
    prodomr <- factor(lookupShape(gr, shp_prod)[[1]], exclude=NULL)
    A <- sparse.model.matrix( ~ prodomr + natura2000 + lovns - 1   )
    A
}

plotMap <- function(..., add=FALSE){
    if(!add)plot(gr,...)
    points(map.pol[,1],map.pol[,2],type="l")
    polygon(map.pol[,1],map.pol[,2],col="grey")
}

## Plot regions on top
plotRegions <- function(...) {
    data(shp, package="mussel")
    plot(shp_natura2000, add=TRUE, ...)
    plot(shp_lovns, add=TRUE, ...)
    omr <- c("32", "33", "34", "36", "37", "38", "39")
    for(x in omr) {
        sub <- subset(shp_prod, OMRådENUMM == x)
        plot(sub, add=TRUE, ...)
        co <- coordinates(sub)
        text(co[1], co[2], x)
    }
}

if(FALSE){ ## Create cache
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
    
    ## Cache default grid, regions etc
    map.pol <- map.pol[8<map.pol[,1] & map.pol[,1]<10 & 56<map.pol[,2] & map.pol[,2]<58  , ]
    save(gr, map.pol, spatialRegions, spatialRegions2, file="data/grid.RData")
    save(shp_lovns,
         shp_natura2000,
         shp_prod, file="data/shp.RData")
}


