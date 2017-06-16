#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SPARSE_MATRIX(Q0);
  DATA_SPARSE_MATRIX(I);
  DATA_FACTOR(time);
  DATA_FACTOR(gf);
  DATA_VECTOR(response);
  DATA_VECTOR(presence);
  DATA_VECTOR(sweptArea);      // by haul
  //DATA_FACTOR(spatialRegions); // by grid index
  DATA_SPARSE_MATRIX(regionIndicator);
  //DATA_SCALAR(gridCellArea);   // km^2
  DATA_VECTOR(yearLevels);     // as.numeric(levels(Year))
  DATA_INTEGER(reportLog);
  /* Random fields */
  PARAMETER_ARRAY(eta_presence);
  PARAMETER_ARRAY(eta_density);
  /* 6 fixed effects x number of times */
  PARAMETER_VECTOR(logdelta);        // Length=2
  PARAMETER_ARRAY(logscale);         // Dim = c(ntime,2)
  PARAMETER_VECTOR(logsd_nugget);    // Length = ntime
  PARAMETER_VECTOR(mu);              // Length = ntime

  vector<Type> sd_nugget=exp(logsd_nugget);
  int nhaul=response.size();
  Type ans=0;
  using namespace density;
  /* Time covariance */
  N01<Type> nldens_time;
  /* Scale parameters for fields */
  vector<Type> scale0=exp(logscale.col(0));
  vector<Type> scale1=exp(logscale.col(1));
  /* GMRF: Sparse precision matrix */
  Eigen::SparseMatrix<Type> Q = Q0+exp(logdelta[0])*I;
  GMRF_t<Type> nldens0=GMRF(Q);
  ans += SEPARABLE(VECSCALE(nldens_time,scale0),nldens0)(eta_presence);
  SIMULATE {
    SEPARABLE(VECSCALE(nldens_time,scale0),nldens0).simulate(eta_presence);
  }
  Eigen::SparseMatrix<Type> QQ= Q0+exp(logdelta[1])*I;
  GMRF_t<Type> nldens1=GMRF(QQ);
  ans += SEPARABLE(VECSCALE(nldens_time,scale1),nldens1)(eta_density);
  SIMULATE {
    SEPARABLE(VECSCALE(nldens_time,scale1),nldens1).simulate(eta_density);
  }
  
  /* Data for presence */
  array<Type> prob(eta_presence.dim);
  prob=Type(1)/(Type(1)+exp(Type(-1)*eta_presence));
  for(int i=0;i<nhaul;i++){
    ans -= dbinom(presence[i],Type(1),prob(gf[i],time[i]),true);
    SIMULATE {
      presence[i] = rbinom(Type(1), prob(gf[i],time[i]));
    }
    if(presence[i]>0){
      Type mu_i = eta_density(gf[i],time[i]) + mu[time[i]] - .5*pow(sd_nugget[time[i]], 2);
      Type sd_i = sd_nugget[time[i]];
      ans -= dnorm(log(response[i]), mu_i, sd_i, true);
      SIMULATE {
        response[i] = exp( rnorm(mu_i, sd_i) );
      }
    } else {
      SIMULATE {
        response[i] = 0;
      }
    }
  }

  /* Report simulated data */
  SIMULATE {
    REPORT(eta_presence);
    REPORT(eta_density);
    REPORT(response);
    REPORT(presence);
  }

  /* report biomass density - notation as in report:

     b = 2.703 * (C/A)^0.29

     Units: C (kg) and A (m^2)
  */
  Type C;
  Type A = sweptArea.mean(); REPORT(A);
  array<Type> b(NLEVELS(spatialRegions), NLEVELS(time));
  b.setZero();
  array<Type> cellcount(NLEVELS(spatialRegions), NLEVELS(time));
  cellcount.setZero();
  // Loop through *all* cells and times
  for(int i=0; i < NLEVELS(gf); i++){
    for(int j=0; j < NLEVELS(time); j++){
      C = exp(eta_density(i,j) + log(1.0/(1.0+exp(-eta_presence(i,j)))) + mu(j));
      C = C / 1000.0; // Kg
      if (yearLevels(j) < 2017 ) {
        // Old Gear (until 2016)
        b(spatialRegions(i), j) += 2.703 * pow(C/A, 0.29);
      } else {
        // New Gear (from 2017)
        b(spatialRegions(i), j) += 2.470 * C/A;
      }
      cellcount(spatialRegions(i), j) = cellcount(spatialRegions(i), j) + 1.0;
    }
  }
  b = b / cellcount;

  if(reportLog){
    array<Type> logb(b);
    logb = log(b);
    REPORT(logb);
    ADREPORT(logb);
  } else {
    REPORT(b);
    ADREPORT(b);
  }

  return ans;

}

