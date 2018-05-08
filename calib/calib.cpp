#include <TMB.hpp>

/*
  Properties:  psi(x) >= x for all x
               Natural to require phi(0)=0
*/
template<class Type>
struct psi_t {
  Type a, b, c;
  Type eval(Type x) {
    return x + a * pow(x, b) + c;
  }
  Type deriv(Type x) {
    return 1 + a * b * pow(x, b - 1);
  }
  Type inverse(Type y0) {
    Type y = y0; // start
    for (int i=0; i<100; i++) {
      y = y - (eval(y) - y0) / deriv(y);
    }
    return y;
  }
};

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
  DATA_VECTOR(x);        // Optional grid to evaluate power function

  PARAMETER_VECTOR(logmu); // Mean 
  PARAMETER_VECTOR(logphi); // Tweedie
  PARAMETER(power); // Tweedie
  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(c);
  /* Notation:

     Din  = Density observation inside trawl path
     Dout = Density observation outside trawl path
     C = Catch density (unobserved) corresponding to same area as 'Din'

     E(C) = E(Dout) - E(Din)

     Model:

     E(Dout) = psi( E(C) )

     Consistency requirement:

     psi(x) >= x

     Implication:

     E(Din) = E(Dout) - E(C) =  E(Dout) - psi.inverse( E(Dout) )
  */
  Type ans = 0;
  for (int i=0; i < Density.size(); i++) {
    Type mu = exp(logmu(Station(i))); // E(Dout)
    psi_t<Type> psi = {a(Gear(i)), b(Gear(i)), c(Gear(i))};
    Type mu_inside = mu - psi.inverse(mu);
    Type mueff = (Inside(i) ? mu_inside : mu);
    Type phi = exp(logphi(AreaFac(i)));
    ans -= dtweedie(Density(i), mueff, phi, power, true);
  }

  if(x.size() == 0) {
    ADREPORT(a);
    ADREPORT(b);
    ADREPORT(c);
  } else {
    psi_t<Type> psi = {a(Gear(0)), b(Gear(0)), c(Gear(0))};
    vector<Type> y(x.size());
    for (int i=0; i<y.size(); i++)
      y(i) = x(i) - psi.inverse(x(i));
    REPORT(x);
    REPORT(y);
    ADREPORT(y);
  }

  return ans;
}
