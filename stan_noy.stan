
//data block
//defines the data to be passed into the model
data {
  int<lower=1> N;  // no. observations   
  int<lower=0> p;  // number of x + 1 for intercept (this needs to be a matrix of 1s.)
  vector[N] y;     // log POC
  matrix[N,p] x;   // matrix of x values
}

//parameter block defines what stan will estimate during sampling
//these have no bounds set on them as of now as they should be allowed to be negative
parameters{ 
  vector[p] beta; // the effects on the mu from each parameter
  vector[p] gamma; // the effects on sigma from each paramter
}

//transformed parameters  block defines the transformations of paramters and data
transformed parameters {
  vector<lower=1e-10>[N]     sigma; //the standard deviation
  for(i in 1:N){ 
  sigma[i] = x[i,]*gamma; 
  }
  vector[N]     mu; //the mean of the distribution
  for(i in 1:N){ 
  mu[i] = x[i,]*beta;
  }
}

//model block. Where the data is modelled
model{
  // Priors //I have not set any priors as of yet
  //?????

  // Likelihood, the function of the model paramters given the observed data
  // If these parameters were true, how likely is the data we already observed?
  for(i in 1:N){
  y[i] ~ normal(mu[i], sigma[i]);
  }
}
