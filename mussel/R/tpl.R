cmb <- function(f, d) function(p) f(p, d)
tpl <- function(parms, data) {
  getAll(data, parms)

  ans <- 0

  ## Hardcode Gear calibration parameters and uncertainties
  calib_mean <- c(0.05470782, 1.81567366)
  ##matrix<Type> calib_hessian(2, 2);
  calib_hessian <- matrix(c(
    12513.5960, 923.45777,
    923.45777,  87.71207),
    2, 2)
  ans <- ans + .5 * sum( calib * (calib_hessian %*% calib) );
  calib <- calib + calib_mean;

  sd_nugget <- exp(logsd_nugget)
  nhaul <- length(response)
  ## Scale parameters for fields
  scale0 <- exp(logscale[,1])
  scale1 <- exp(logscale[,2])
  ## Time covariance
  ## N01<Type> nldens_time;
  ldens_time0 <- function(x) sum(dnorm(x, sd=scale0, log=TRUE))
  ## GMRF: Sparse precision matrix
  ## Eigen::SparseMatrix<Type> Q = Q0+exp(logdelta[0])*I;
  Q <- Q0 + exp(logdelta[1]) * I
  ##GMRF_t<Type> nldens0=GMRF(Q);
  ldens0 <- function(x) dgmrf(x, 0, Q, log=TRUE)
  ##ans += SEPARABLE(VECSCALE(nldens_time,scale0),nldens0)(eta_presence);
  ans <- ans - dseparable(ldens0, ldens_time0)(eta_presence, log=TRUE)
  
  ## SIMULATE {
  ##   SEPARABLE(VECSCALE(nldens_time,scale0),nldens0).simulate(eta_presence);
  ## }
  ##Eigen::SparseMatrix<Type> QQ= Q0+exp(logdelta[1])*I;
  QQ = Q0 + exp(logdelta[2]) * I
  ##GMRF_t<Type> nldens1=GMRF(QQ);
  ldens1 <- function(x) dgmrf(x, 0, QQ, log=TRUE)
  ldens_time1 <- function(x) sum(dnorm(x, sd=scale1, log=TRUE))
  ##ans += SEPARABLE(VECSCALE(nldens_time,scale1),nldens1)(eta_density);
  ans <- ans - dseparable(ldens1, ldens_time1)(eta_density, log=TRUE)
  ## SIMULATE {
  ##   SEPARABLE(VECSCALE(nldens_time,scale1),nldens1).simulate(eta_density);
  ## }
  
  ## Data for presence
  ##array<Type> prob(eta_presence.dim);
  ##prob=Type(1)/(Type(1)+exp(Type(-1)*eta_presence));

  prob <- 1 / (1 + exp(-eta_presence)) ## dim = space x time

  ans <- ans - sum(dbinom(presence, 1, prob[cbind(gf,time)], log=TRUE))

  i <- (presence > 0)
  ## Presence
  mu_i <- eta_density[cbind(gf[i],time[i])] + mu[time[i]] - .5 * sd_nugget[time[i]]^2
  sd_i <- sd_nugget[time[i]]
  ans <- ans - sum(dnorm(log(response[i]), mu_i, sd_i, log=TRUE))
  ## Absence
  response[i] <- 0
  
  ## for(int i=0;i<nhaul;i++){
  ##   ans -= dbinom(presence[i],Type(1),prob(gf[i],time[i]),true);
  ##   SIMULATE {
  ##     presence[i] = rbinom(Type(1), prob(gf[i],time[i]));
  ##   }
  ##   if(presence[i]>0){
  ##     Type mu_i = eta_density(gf[i],time[i]) + mu[time[i]] - .5*pow(sd_nugget[time[i]], 2);
  ##     Type sd_i = sd_nugget[time[i]];
  ##     ans -= dnorm(log(response[i]), mu_i, sd_i, true);
  ##     SIMULATE {
  ##       response[i] = exp( rnorm(mu_i, sd_i) );
  ##     }
  ##   } else {
  ##     SIMULATE {
  ##       response[i] = 0;
  ##     }
  ##   }
  ## }

  ## /* Report simulated data */
  ## SIMULATE {
  ##   REPORT(eta_presence);
  ##   REPORT(eta_density);
  ##   REPORT(response);
  ##   REPORT(presence);
  ## }

  ## report biomass density - notation as in report:
  ##    b = 2.703 * (C/A)^0.29
  ##    Units: C (kg) and A (m^2)
  ##    New formulation
  ##    ^^^^^^^^^^^^^^^
  ##    b = 0.05 * (C/A)^1.82 + (C/A)^1
  ##
  
  ##Type C;
  ##Type A = sweptArea.mean(); REPORT(A);
  A <- mean(sweptArea); REPORT(A)
  ## matrix<Type> b_full(NLEVELS(gf), NLEVELS(time));
  ## // Loop through *all* cells and times
  ## for(int i=0; i < NLEVELS(gf); i++){
  ##   for(int j=0; j < NLEVELS(time); j++){
  ##     C = exp(eta_density(i,j) + log(1.0/(1.0+exp(-eta_presence(i,j)))) + mu(j));
  ##     C = C / 1000.0; // Kg
  ##     b_full(i, j) = calib(0) * pow(C/A, calib(1)) + C/A;
  ##     //cellcount(k, j) = cellcount(k, j) + 1.0;
  ##   }
  ## }
  C <- exp(eta_density + log( 1.0 / (1.0 + exp(-eta_presence ))) + rep(mu,each=nlevels(gf)) )
  C <- C / 1000.0 ## Kg
  b_full <- calib[1] * (C/A)^calib[2] + C/A
  
  ## Aggregate biomass across regions
  ## int nreg = regionIndicator.cols();
  ## matrix<Type> b(nreg, NLEVELS(time));
  ## b = regionIndicator.transpose() * b_full;
  b <- t(regionIndicator) %*% b_full
  
  ## Divide with cell count withing region
  ## matrix<Type> one(NLEVELS(gf), NLEVELS(time));
  ## one.fill(Type(1));
  ## matrix<Type> cellcount(nreg, NLEVELS(time));
  ## cellcount = regionIndicator.transpose() * one;
  ## b = b.array() / cellcount.array();
  b <- b / colSums(regionIndicator)

  if(reportLog){
    ## matrix<Type> logb(b);
    ## logb = log(b.array());
    logb <- log(b)
    REPORT(logb);
    ADREPORT(logb);
  } else {
    REPORT(b);
    ADREPORT(b);
  }

  ans

}
