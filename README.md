# Dirichlet Laplace priors for optimal shrinkage: an implementation using R

The R code suggests an implementation of the 2014 article by Anirban Bhattacharya, Debdeep Pati, Natesh S. Pillai, David B. Dunson.

It includes both the implementation of the simulation part and the use of the algorithm on Prostate Cancer Data 

Model: y ~ N(theta, I_n)

prior for theta: theta_j ~ N(0, psi_j phi_j^2 tau^2), psi_j i.i.d. exp(1/2),

phi ~ Dir(a,..., a), tau ~ gamma(n*a, 1/2)
