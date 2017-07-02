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

f0 <- function(x)coef(fm0) * x
f1 <- function(x)exp(coef(fm1)[1]) * x^coef(fm1)[2]
f2 <- function(x)coef(fm2)[2]*x+ coef(fm2)[1] +2*prod(coef(fm2))*sqrt(x)
f3 <- function(x)2.703*x^0.29
plot(f0,0,15)
plot(f1,0,15,add=TRUE,col=2)
plot(f2,0,15,add=TRUE,col=3)
plot(f3,0,15,add=TRUE,col=4)
## abline(v=range(df$inside, na.rm=TRUE))
## points(df, pch=16)
points(df2, pch=16, col=grepl(" ", rownames(df))+1)


## Plot in sqrt-transformed domain
tplot <- function(f,a,b,...){plot(function(x)sqrt(f(x^2)), sqrt(a), sqrt(b), ..., xlab="sqrt(x)", ylab="sqrt(y)", lwd=2)}
tplot(f0,0,15)
tplot(f1,0,15,add=TRUE,col=2)
tplot(f2,0,15,add=TRUE,col=3)
tplot(f3,0,15,add=TRUE,col=4)
## abline(v=range(df$inside, na.rm=TRUE))
##points(sqrt(df), pch=16)
##points(sqrt(df2), pch=16, col="red")
points(sqrt(df), pch=16, col=grepl(" ", rownames(df))+1)

###################### TMB
data1 <- with(dd, data.frame(
                      Density = Total / Area,
                      Station = Station,
                      AreaFac = cut(Area, breaks=c(0,1,200)),
                      Inside = as.integer(Dredge == "inside"),
                      Gear = Gear
                  ))
levels(data1$Station) <- paste(levels(data1$Station), "Total")
data1 <- na.omit(data1)

data2 <- with(dd, data.frame(
                      Density = BMS / Area,
                      Station = Station,
                      AreaFac = cut(Area, breaks=c(0,1,200)),
                      Inside = as.integer(Dredge == "inside"),
                      Gear = Gear
                  ))
levels(data2$Station) <- paste(levels(data2$Station), "BMS")
data2 <- na.omit(data2)


data <- rbind(data1, data2)
data <- lapply(data, function(x) if(is.factor(x)) factor(x) else x)

parameters <- with(data,
                   list(
                       logmu  = numeric(nlevels(Station)),
                       logphi = numeric(nlevels(AreaFac)),
                       power = 1.5,
                       a = 1 + numeric(nlevels(Gear)),
                       b = 1 + numeric(nlevels(Gear))
                   ))

require(TMB)
compile("calib.cpp")
dyn.load(dynlib("calib"))

################################################################################

obj <- MakeADFun(data, parameters)
fit <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)

print(rep)
cat(paste0(fit$objective, " (df=",length(fit$par),")\n") )

## Joint: 701.223884983253 (df=32)
## Data1: 391.7284 (df=14)
## Data2: 294.97863087435 (df=23)
