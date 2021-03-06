## ========================== OLD GEAR CALIB DATA

d <- read.table("calib_old_gear.csv",sep=",",header=TRUE, check.names = FALSE)
d <- d[-(9:10)] ## Remove subsampling info

## Reshape dataset

## Extract diver data
head(d)
d$is_kg <- as.character(d$Station) %in% LETTERS[1:10]
d1 <- with(d, data.frame(
                  Total = `Total (g)` / ifelse(is_kg, 1, 1000),
                  BMS = `BMS (g)` / ifelse(is_kg, 1, 1000),
                  Station = Station,
                  Area = `Ring area`
              )
           )
d1$Dredge <- "outside"

## Extract dredge data
head(d)
d2 <- with(d, data.frame(
                  Total = `Total fangst (kg)`,
                  BMS = `BMS (kg)`,
                  Station = Station,
                  Area = `Area`
              )
           )
d2$Dredge <- "inside"
d2 <- subset(d2, !is.na(Total))

dd <- rbind(d1, d2)
dd$Gear <- "old"

## ========================== NEW GEAR CALIB DATA

d <- read.csv("calib_new_gear.csv", header=TRUE, check.names=FALSE, na.strings = "na")
names(d)[3] <- "Pos_within_dredge"
d$Rectangle <- factor(1:nrow(d))
dnew <- do.call("rbind", rep(list(d[-(4:7)]), 4))
dnew$haul.id <- factor(paste(dnew$Area, dnew$Dredge))
v <- do.call("c", d[4:7])
dnew$box <- substring(names(v), 1,1)
dnew$inside <- sub("\\).*","",sub(".*\\(","",names(v))) == "inside"
dnew$outside <- sub("\\).*","",sub(".*\\(","",names(v))) == "outside"
dnew$weight <- v
dnew <- na.omit(dnew)
names(dnew) <- sub(" .*", "", names(dnew))
dnew[1:6] <- lapply(dnew[1:6], factor)

d3 <- with(dnew, data.frame(
                     Total = NA, ## FIXME !!!!
                     BMS = `weight` / 1000,
                     Station = haul.id,
                     Area = 0.25,
                     Dredge = ifelse(inside, "inside", "outside")
                 )
           )

d3$Gear <- "new"

dd <- rbind(dd, d3)

## =========================== NEW GEAR TOTAL CATCHES

## Area 35
tot <- c(0, 61, 93, 89, 121)
bms <- c(0, 3.65, 4.2, 3.9, 14.15)
shells <- c(0, 0.8, 0.75, 1.85, 2.57)
bms35 <- tot * bms / (bms + shells)

## Area 39
tot <- c(253, 154, 306, 354)
bms <- c(7.70, 7.80, 7.75, 7.80)
shells <- c(0.200, 0.150, 0.200, 0.200)
bms39 <- tot * bms / (bms + shells)

d4 <- data.frame(
    Total = NA, ## FIXME !!!!
    BMS = c(bms35[-1], bms39),
    Station = c(paste(35, 1:4), paste(39, 1:4)),
    Area = 100,
    Dredge = "inside"
)
d4$Gear <- "new"
dd <- rbind(dd, d4)

## library(glmmTMB)
## fit <- glmmTMB( Total ~ Station + offset(log(Area)) + Gear - 1, data=dd, family=poisson(link="log") )

y <- dd$Total / dd$Area
tab <- tapply(y, list(dd$Dredge, dd$Station), mean, na.rm=TRUE)
df1 <- as.data.frame(t(tab))
plot(df1)

y <- dd$BMS / dd$Area
tab <- tapply(y, list(dd$Dredge, dd$Station), mean, na.rm=TRUE)
df2 <- as.data.frame(t(tab))
points(df2,col="red")

df <- rbind(df1,df2)

plot(sqrt(df), pch=16)
points(sqrt(df2), col="red" , pch=16)


plot(log(df), pch=16)
points(log(df2), col="red" , pch=16)


summary(fm0 <- lm(outside~inside-1,data=df))
summary(fm1 <- lm(log(outside)~log(inside),data=df))
summary(fm2 <- lm(sqrt(outside)~sqrt(inside),data=df))
fm4 <- nls(log(outside) ~ nlb + log(inside) - log(1-exp(-a*inside)) ,start=c(nlb=0, a=1),data=na.omit(df))

f0 <- function(x)coef(fm0) * x
f1 <- function(x)exp(coef(fm1)[1]) * x^coef(fm1)[2]
f2 <- function(x)coef(fm2)[2]*x+ coef(fm2)[1] +2*prod(coef(fm2))*sqrt(x)
f3 <- function(x)2.703*x^0.29
f4 <- function(x){b <- exp(-coef(fm4)["nlb"]); a <- coef(fm4)["a"]; x/(b*(1-exp(-a*x)))}
plot(f0,0,15)
plot(f1,0,15,add=TRUE,col=2)
plot(f2,0,15,add=TRUE,col=3)
plot(f3,0,15,add=TRUE,col=4)
plot(f4,0,15,add=TRUE,col=5)
## abline(v=range(df$inside, na.rm=TRUE))
## points(df, pch=16)
points(df2, pch=16, col=grepl(" ", rownames(df))+1)


## Plot in sqrt-transformed domain
tplot <- function(f, a, b, ...) {
    plot(function(x)sqrt(f(x^2)), sqrt(a), sqrt(b), ..., xlab="Density inside (kg/m^2)", ylab="Density outside (kg/m^2)", lwd=2, axes=FALSE)
    x <- 0:10
    axis(1, at=x, labels=(x)^2)
    axis(2, at=x, labels=(x)^2, las=1)
    box()
}
tplot(f0,0,15)
tplot(f1,0,15,add=TRUE,col=2)
tplot(f2,0,15,add=TRUE,col=3)
tplot(f3,0,15,add=TRUE,col=4)
tplot(f4,0,15,add=TRUE,col=5)
## abline(v=range(df$inside, na.rm=TRUE))
##points(sqrt(df), pch=16)
##points(sqrt(df2), pch=16, col="red")
points(sqrt(df), pch=16, col=grepl(" ", rownames(df))+1)

###################### TMB
refac <- function(data) data.frame(lapply(na.omit(data), function(x) if(is.factor(x)) factor(x) else x))

data1 <- with(dd, data.frame(
                      Density = Total / Area,
                      Station = Station,
                      AreaFac = cut(Area, breaks=c(0,1,200)),
                      Inside = as.integer(Dredge == "inside"),
                      Gear = Gear
                  ))
levels(data1$Station) <- paste(levels(data1$Station), "Total")
levels(data1$Gear)    <- paste(levels(data1$Gear),    "Total")
levels(data1$AreaFac) <- paste(levels(data1$AreaFac), "Total")
data1 <- refac(data1)

data2 <- with(dd, data.frame(
                      Density = BMS / Area,
                      Station = Station,
                      AreaFac = cut(Area, breaks=c(0,1,200)),
                      Inside = as.integer(Dredge == "inside"),
                      Gear = Gear
                  ))
levels(data2$Station) <- paste(levels(data2$Station), "BMS")
levels(data2$Gear)    <- paste(levels(data2$Gear),    "BMS")
levels(data2$AreaFac) <- paste(levels(data2$AreaFac), "BMS")
data2 <- refac(data2)


data <- refac(rbind(data1, data2))

parameters <- function(data)
    with(data,
         list(
             logmu  = numeric(nlevels(Station)),
             logphi = numeric(nlevels(AreaFac)),
             power = 1.5,
             a = 1 + numeric(nlevels(Gear)),
             b = 1 + numeric(nlevels(Gear)),
             c = 0 + numeric(nlevels(Gear))
         ))

require(TMB)
compile("calib.cpp")
dyn.load(dynlib("calib"))

################################################################################

data <- as.list(data); data$x <- numeric(0)
obj <- MakeADFun(data, parameters(data))
fit <- nlminb(obj$par + .01, obj$fn, obj$gr)
rep <- sdreport(obj)

print(rep)
cat(paste0(fit$objective, " (df=",length(fit$par),")\n") )
fit1 <- fit

## Joint: 686.877140102969 (df=36)
## Data1: 391.7284 (df=14)
## Data2: 294.97863087435 (df=23)

## Collect Gears
levels(data$Gear)[] <- "collect"
map <- list(c=factor(NA))
obj <- MakeADFun(data, parameters(data), map=map)
fit <- nlminb(obj$par + .01, obj$fn, obj$gr, control=list(iter.max=1e4,eval.max=1e4))
rep <- sdreport(obj)

print(rep)
cat(paste0(fit$objective, " (df=",length(fit$par),")\n") )
fit2 <- fit

## p value for H0: identical gears
1-pchisq(2*(fit2$objective-fit1$objective), df=length(fit1$par)-length(fit2$par))

tplot(function(x)(1/fit2$par["a"])^(1/(fit2$par["b"]+1))*x^(1/(fit2$par["b"]+1)), 0, 36, add=!TRUE, col="black")
tplot(f1,0,36,add=TRUE,col=2)
tplot(f3,0,36,add=TRUE,col=4)
##points(sqrt(df), pch=16, col=grepl(" ", rownames(df))+1)
points(sqrt(df), pch=c(1,16)[grepl(" ", rownames(df))+1])
legend("topleft",c("Tweedie", "Gaussian", "Dolmer"), lwd=2, col=c(1,2,4))
abline(0,1,lty="dashed")

sdr <- sdreport(obj)
summary(sdr,"report")
nm <- 1:2
quadform <- list(
    mu = sdr$value[nm],
    hessian = solve(sdr$cov[nm, nm])
)
