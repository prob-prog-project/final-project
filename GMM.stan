data{
  int<lower=0> C;
  int<lower=0> G;

  real infer[C,G];

  real prior_mu[3];



}

parameters{
  simplex[3] pi; #mixing proportion
#  vector<lower=0>[3] mu; #lower bound at 0 because all means must be >0
  vector<lower=0>[3] sigma;
}

model{
  vector[3] contributions;

#  mu ~ normal(prior_mu,.5);
  sigma ~ normal(0,.5);

  pi ~ dirichlet([1.,1.,1.]');

  for (c in 1:C){
    for (g in 1:G){
      for (k in 1:3){
        contributions[k] = log(pi[k])+normal_lpdf(infer[c][g] | prior_mu[k], sigma[k]) ;
      }
      target += log_sum_exp(contributions) ;
    }
  }

}
