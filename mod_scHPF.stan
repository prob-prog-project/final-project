data {
  int<lower=0> C;
  int<lower=0> G;

  real<lower=0> a_prime;
  real<lower=0> b_prime;
  real<lower=0> a;

  real<lower=0> c_prime;
  real<lower=0> eta_sigma;
  real<lower=0> c;

  int<lower=0> k;

  int<lower=0> X[C,G];

  real<lower=0> mixture_mu[C,G];



}

parameters {
  vector<lower=0> [C] theta;
  matrix<lower=0> [C,G] beta; //for now removed the k
  matrix<lower=0>[ C,G] eta;
  vector<lower=0>[C] xi;

}
model {

for (cc in 1:C){
  xi[cc] ~ gamma(a_prime,b_prime);
#  for (kk in 1:k){
    theta[cc] ~ gamma(a,xi[cc]);
#  }
}

  for (g in 1:G){
    for (cc in 1:C){
      eta[cc,g] ~ normal(mixture_mu[cc,g],eta_sigma);
    #  for (kk in 1:k){
        beta[cc,g]~ gamma(c,eta[cc,g]);
    #  }
    }
  }

  for (cc in 1:C){
    for (g in 1:G){
      X[cc,g] ~ poisson((theta[cc]*beta[cc,g]));
    }
  }
}
