data {
  int<lower=1> T; // Length of time series
  int<lower=1> K; // Number of response categories
  int<lower=1, upper=K> Y[T]; // Observations
  int<lower=1> Period; // Known period of seasonal component
}

transformed data {
  real log_unif = -log(T);
}

parameters {
  real A;
  real D;
  real alpha0;
  real alpha1;
  real Delta; // level shift at changepoint
  ordered[K - 2] threshold;
  real<lower=-1,upper=1> rho;
}

transformed parameters {
  vector[K-1] cutpoints = append_row(0, threshold); // for identifiability
  vector[T] lp;
  lp = rep_vector(log_unif, T);
  {
    for (cp in 1:T) {
      vector[T] Z;
      vector[T] mu;
      for (t in 1:T) {
        real s;
        s = A * cos(2 * pi() / Period * t) + D * sin(2 * pi() / Period * t);
        mu[t] = alpha0 + alpha1 * (t * 1.0 / T) + s + (t >= cp ? Delta : 0);
        if (t == 1) {
          Z[t] = mu[t];
        } else {
          Z[t] = mu[t] + rho * (Z[t - 1] - mu[t - 1]);
        }
        lp[cp] = lp[cp] + ordered_probit_lpmf(Y[t] | Z[t], cutpoints);
      }
    }
  }
}

model {
  A ~ normal(0, 3);
  D ~ normal(0, 3);
  alpha0 ~ normal(0, 3);
  alpha1 ~ normal(0, 3);
  Delta ~ normal(0, 3);
  threshold ~ normal(0, 3);
  target += log_sum_exp(lp);
}
