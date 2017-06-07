plotBiomassMap <- function(env, time="2016", xlim=c(9.15,9.35), ylim=c(56.6,56.7) ) {
    env$showTime <- time
    env$xlim <- xlim
    env$ylim <- ylim
    local({
        y <- which(levels(data$time)==showTime)
        logb <- logB[,y]
        resp <- exp(logb)
        plot(gr,type="n",las=1,xlab="Longitude",ylab="Latitude",xlim=xlim , ylim=ylim)
        subresp <- resp[xlim[1]<gr$lon & gr$lon<xlim[2] & ylim[1]<gr$lat & gr$lat<ylim[2]]
        probs <- seq(0,1,length=26)
        br <- quantile(subresp,probs)
        col <- rev(topo.colors(25))
        image(gr, resp , map=quote(plotMap(add=TRUE)),  add=TRUE, breaks=br, col=col)
        i <- seq(1,26,by=5)
        qw <- subset(d,unclass(time) == y)
        if (FALSE) text(qw$lon,qw$lat+c(-.01,0,.01)/3,qw$TotalWeight,col="red")
        points(qw$lon,qw$lat)
        par(usr=c(0,1,0,1));axis(4, at=probs[i], label=round(br,2)[i], las=1, hadj=.3)
        box()
        title(paste("Time:",showTime), outer=TRUE, line=-1)
    }, env)
}
