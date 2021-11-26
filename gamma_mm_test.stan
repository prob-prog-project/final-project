// data{
//
//   int<lower=0> C;
//   int<lower=0> G;
//
//   real<lower=0> eta[C,G];
//
//
//   vector<lower=0>[3] priors;
//
//   //real beta[3]; could be here if we want to set it
//
//
//   vector[3] a_prime;
//   vector[3] b_prime;
//
// #  vector[3] mu_alpha;
// #  vector[3] mu_beta;
//
//   vector[3] a;
//   vector[3] b;
// }
//
// parameters{
//   simplex[3] pi; //mixing proportion
//
//   real<lower=0> beta[3];
//   real<lower=0>  alpha[3];
// }
//
//
//
// model{
//   real ps[3];
//
//   pi ~ dirichlet([1,1,1]');
//
//   for (j in 1:3){
//     alpha[j] ~ gamma(a_prime[j],b_prime[j]);
//     beta[j] ~ gamma(a[j],b[j]);
//   }
//
//   for (c in 1:C){
//     for (g in 1:G){
//       for (k in 1:3)
//         ps[k] = log(pi[k]) + gamma_lpdf(eta[c][g] | 100000*alpha[k], 100000*beta[k]);
//       target += log_sum_exp(ps);
//     }
//   }
// }
//
// generated quantities{
//   int mixture[C,G];
//   real ps[3];
//   for (c in 1:C){
//     for (g in 1:G){
//       for (k in 1:3){
//         ps[k] = log(pi[k]) + log(gamma_cdf(eta[c][g] | alpha[k], beta[k]));
//       }
//       real max = 0;
//       int max_i = 0;
//       for (k in 1:3){
//         if(ps[k] > max) {
//           mixture[c][g] = k;
//           max = ps[k];
//         }
//       }
//     }
//   }
//
// }
