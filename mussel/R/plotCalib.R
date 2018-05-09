plotCalib <- function(type=c("data", "efficiency", "efficiency2"), ...) {
    type <- match.arg(type)
    data(calib, package="mussel")
    ## Rough inside/outside estimates based on total catches
    y <- calib$Total / calib$Area
    tab <- tapply(y, list(calib$Dredge, calib$Station), mean, na.rm=TRUE)
    df1 <- as.data.frame(t(tab))
    ## Rough inside/outside estimates based on BMS catches
    y <- calib$BMS / calib$Area
    tab <- tapply(y, list(calib$Dredge, calib$Station), mean, na.rm=TRUE)
    df2 <- as.data.frame(t(tab))
    ## Join them
    df <- rbind(df1, df2)
    names(df) <- paste(names(df),"(kg/m^2)")
    df$Gear <- factor(c("old","new")[grepl(" ", rownames(df))+1])
    ## Fit
    a <- 0.05470782; b <- 1.81567366
    psi <- function(x) a*x^b + x  ## Maps 'Catch inside' to 'density outside'
    invpsi <- Vectorize(function(y) {
        uniroot(function(x)psi(x)-y, c(0,100))$root
    })
    fnew <- function(x) x - invpsi(x) ## Maps 'outside' to 'inside'
    ## Plot
    if(type=="data"){
        plot(df[1:2], log="xy", pch=c(16,1)[unclass(df$Gear)],
             main="Gear calibration", ...)
        legend("topleft",
               c("Samples new gear", "Samples old gear", "Model fit"),
               pch=c(16, 1, NA), lwd=c(NA,NA,2), col=c(1,1,1) )
        ygrid <- seq(1e-3,30,length=101)
        xval <- fnew(ygrid)
        lines(xval, ygrid, lwd=2)
    }
    ## Efficiency
    if(type=="efficiency") {
        eff <- function(C) C / psi(C)
        plot(eff, 0, 20, ylim=c(0,1), lwd=2,
             ylab="Gear efficiency", xlab="Catch (kg/m^2)", ...)
    }
    ## As fct of 'Density outside'
    if(type=="efficiency2") {
        plot(df$outside, 1-df$inside/df$outside, ylim=c(0,1),
             ylab="Gear efficiency", xlab="Density outside (kg/m^2)",
             ...)
        eff2 <- function(Dout) invpsi(Dout) / Dout
        plot(eff2, 0, 25, add=TRUE, lwd=2)
    }
    NULL
}
