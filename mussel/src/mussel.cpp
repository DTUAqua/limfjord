#define TMB_LIB_INIT R_init_mussel
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
  /* Gear calibration parameters (mle=zero) */
  PARAMETER_VECTOR(calib);

  Type ans=0;

  // Hardcode Gear calibration parameters and uncertainties
  vector<Type> calib_mean(2);
  calib_mean << 0.05470782, 1.81567366;
  matrix<Type> calib_hessian(2, 2);
  calib_hessian <<
    12513.5960, 923.45777,
    923.45777,  87.71207;
  ans += .5 * ( calib * (calib_hessian * calib) ).sum();
  calib = calib + calib_mean;

  vector<Type> sd_nugget=exp(logsd_nugget);
  int nhaul=response.size();
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

     New formulation
     ^^^^^^^^^^^^^^^

     b = 0.05 * (C/A)^1.82 + (C/A)^1

  */
  Type C;
  Type A = sweptArea.mean(); REPORT(A);
  matrix<Type> b_full(NLEVELS(gf), NLEVELS(time));
  // Loop through *all* cells and times
  for(int i=0; i < NLEVELS(gf); i++){
    for(int j=0; j < NLEVELS(time); j++){
      C = exp(eta_density(i,j) + log(1.0/(1.0+exp(-eta_presence(i,j)))) + mu(j));
      C = C / 1000.0; // Kg
      b_full(i, j) = calib(0) * pow(C/A, calib(1)) + C/A;
      //cellcount(k, j) = cellcount(k, j) + 1.0;
    }
  }

  // Aggregate biomass across regions
  int nreg = regionIndicator.cols();
  matrix<Type> b(nreg, NLEVELS(time));
  b = regionIndicator.transpose() * b_full;

  // Divide with cell count withing region
  matrix<Type> one(NLEVELS(gf), NLEVELS(time));
  one.fill(Type(1));
  matrix<Type> cellcount(nreg, NLEVELS(time));
  cellcount = regionIndicator.transpose() * one;
  b = b.array() / cellcount.array();

  if(reportLog){
    matrix<Type> logb(b);
    logb = log(b.array());
    REPORT(logb);
    ADREPORT(logb);
  } else {
    REPORT(b);
    ADREPORT(b);
  }

  return ans;

}

