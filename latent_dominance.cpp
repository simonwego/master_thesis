#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);       // length s
  DATA_IVECTOR(x);      // length K
  DATA_IVECTOR(i);      // length K (0-indexed!)
  DATA_IVECTOR(j);      // length K (0-indexed!)
  DATA_IVECTOR(n);      // length K
  DATA_INTEGER(use_c);

  PARAMETER(log_a);
  PARAMETER(b);
  PARAMETER(c);
  PARAMETER_VECTOR(z);  // length s

  const int K = x.size(); 
  Type cc = use_c ? c : Type(0.0);
  Type nll = Type(0.0);
  Type a = exp(log_a);
  vector<Type> p(K);

  for (int k = 0; k < K; k++) {
    const int ii = i[k];
    const int jj = j[k];

    Type eta =
      a * (z[ii] - z[jj]) +
      b * (y[ii] - y[jj]) +
      cc * (z[ii] * y[jj] - z[jj] * y[ii]);

    p[k] = Type(1) / (Type(1) + exp(-eta));
    nll -= dbinom(Type(x[k]), Type(n[k]), p[k], true);
  }
  nll -= sum(dnorm(z, Type(0.0), Type(1.0), true));
  nll -= dnorm(b,Type(0.0),Type(1.0),true);
  nll -= dnorm(c,Type(0.0),Type(1.0),true);
  nll -= dnorm(log_a,Type(0.0),Type(1.0),true);
  ADREPORT(a);
  ADREPORT(b);
  ADREPORT(c);
  REPORT(p);

  return nll;
}
