#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Density);  // Density kg/m^2
  DATA_IVECTOR(Inside);  // Logical; Inside/Outside dredge ?
  DATA_FACTOR(Station);  // Mean constant within station. If we have
                         // multiple categories
                         // (BMS/Other/Total/Oyster) we paste it with
                         // the station code.
  DATA_FACTOR(AreaFac);  // Variance depends on area
  DATA_FACTOR(Gear);     // Which gear is it ?
  
  PARAMETER_VECTOR(logmu); // Mean 
  PARAMETER_VECTOR(logphi); // Tweedie
  PARAMETER(power); // Tweedie
  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(b);
  Type ans = 0;
  for (int i=0; i < Density.size(); i++) {
    Type mu = exp(logmu(Station(i)));
    Type eff = a(Gear(i)) * pow(mu, b(Gear(i)));
    Type mueff = (Inside(i) ? mu * eff : mu);
    Type phi = exp(logphi(AreaFac(i)));
    ans -= dtweedie(Density(i), mueff, phi, power, true);
  }

  b = Type(1) / (Type(1) + b);
  a = Type(1) / a;
  a = pow(a, b);
  ADREPORT(a);
  ADREPORT(b);

  return ans;
}
