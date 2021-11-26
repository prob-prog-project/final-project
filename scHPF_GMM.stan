data {
int<lower=0> C;
int<lower=0> G;

real<lower=0> a_prime;
real<lower=0> b_prime;
real<lower=0> a;

real<lower=0> c_prime;
real<lower=0> d_prime;
real<lower=0> c;

int<lower=0> k;

int<lower=0> X[C,G];


real infer[C,G];

real prior_mu[3];







}

parameters {
simplex[3] pi; #mixing proportion
vector<lower=0>[3] sigma;
matrix<lower=0> [C,k] theta;
matrix<lower=0> [G,k] beta;
vector<lower=0>[G] eta;
vector<lower=0>[C] xi;



}
model {

  vector[3] contributions;

#  mu ~ normal(prior_mu,.5);
  sigma ~ normal(0,.5);

  pi ~ dirichlet([1.,1.,1.]');

  for (cc in 1:C){
    for (g in 1:G){
      for (kk in 1:3){
        contributions[kk] = log(pi[kk])+normal_lpdf(infer[cc][g] | prior_mu[kk], sigma[kk]) ;
      }
      target += log_sum_exp(contributions) ;
    }
  }

  int mixture[C,G] = rep_array(1, C, G);
  real ps[3];
  for (cc in 1:C){
    for (g in 1:G){
      for (kk in 1:3){
        ps[kk] = log(pi[kk]) + (normal_lpdf(infer[cc][g] | prior_mu[kk], sigma[kk]));
        if(ps[kk] >= ps[mixture[cc][g]]){
          mixture[cc][g] = kk;
        }
      }


    }
  }

  //mixture holds the mean for each value

  for (cc in 1:C){
    xi[cc] ~ gamma(a_prime,b_prime);
    for (kk in 1:k){
      theta[cc,kk] ~ gamma(a,xi[cc]);
    }
  }

  for (g in 1:G){
    eta[g] ~ gamma(c_prime,d_prime);
    for (kk in 1:k){
      beta[g,kk] ~ gamma(c,eta[g]);
    }
  }

  for (cc in 1:C){
    for (g in 1:G){
      X[cc,g] ~ poisson(dot_product(theta[cc],beta[g]));
    }
  }
}
