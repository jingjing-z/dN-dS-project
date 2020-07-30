functions {
  // substitution rate matrix
	matrix PDRM(real mu, real kappa, real omega) {
	  matrix[61,61] M;
	  int r1[31];
	  int r2[138];
	  int r3[58];
	  int r4[36];
	  int c1[31];
	  int c2[138];
	  int c3[58];
	  int c4[36];
	  
	  M = rep_matrix(0.,61,61);
	  
	  r1 = {1,3,3,4,5,7,9,11,14,16,18,20,22,24,26,28,30,34,36,38,40,42,44,46,48,50,52,54,56,58,60};
	  r2 = {1,1,1,1,1,1,2,2,2,2,2,2,3,3,4,4,4,5,5,5,5,6,6,6,6,7,7,8,8,8,9,9,10,10,11,11,11,12,12,12,13,13,14,14,14,14,15,15,15,15,16,16,16,16,17,17,17,17,18,18,18,18,19,19,19,19,20,20,20,20,21,21,21,21,22,22,22,22,23,23,23,23,24,24,25,25,26,26,27,27,28,29,30,30,30,31,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,42,42,43,43,46,46,47,47,48,48,49,49,50,50,51,51,52,52,53,53,54,54,55,55};
	  r3 = {1,1,2,2,3,4,5,6,7,8,9,9,10,10,11,12,13,14,15,16,17,22,23,24,25,30,30,31,31,32,32,32,33,33,34,35,36,37,38,38,39,39,40,40,41,41,42,43,44,45,46,47,48,49,54,55,56,57};
	  r4 = {5,5,6,6,14,14,15,15,18,18,19,19,26,26,27,27,28,29,30,31,34,34,35,35,46,46,47,47,50,50,51,51,58,58,59,59};
	  c1 = {2,4,16,17,6,8,10,12,15,17,19,21,23,25,27,29,31,35,37,39,41,43,45,47,49,51,53,55,57,61};
	  c2 = {3,4,9,11,30,46,3,4,10,12,31,47,32,48,13,33,39,9,11,34,50,10,12,35,51,36,52,13,37,53,38,54,39,55,13,42,58,13,43,59,45,61,22,26,30,46,23,27,31,47,24,28,32,48,25,29,33,49,22,26,34,50,23,27,35,51,24,28,36,52,25,29,37,53,24,25,38,54,24,25,39,55,40,56,41,57,42,58,43,59,60,61,33,38,42,33,39,43,40,44,41,45,38,42,39,43,40,44,41,45,40,41,40,41,44,45,44,45,54,58,55,59,56,60,57,61,54,58,55,59,56,60,57,61,56,57,56,57};
	  c3 = {5,14,6,15,7,8,18,19,20,21,11,22,12,23,26,27,29,18,19,20,21,26,27,28,29,34,46,35,47,33,36,48,37,49,50,51,52,53,42,54,43,55,44,56,45,57,58,59,60,61,50,51,52,53,58,59,60,61};
	  c4 = {7,8,7,8,16,17,16,17,20,21,20,21,28,29,28,29,44,45,32,32,36,37,36,37,48,49,48,49,52,53,52,53,60,61,60,61};
	  
	  for (i in 1:31) {
	    M[r1[i],c1[i]] = kappa*mu;
	  }
	  
	  for (i in 1:138){
	    M[r2[i],c2[i]] = omega*mu;
	  }
	  
    for (i in 1:58){
      M[r3[i],c3[i]] = kappa*omega*mu;
    }
    
    for (i in 1:36){
      M[r4[i],c4[i]] = mu;
    }
    
    
    // Fill in the lower triangle
    M = M'+ M;
  
    // Compute the diagonal
    for (i in 1:61){
      M[i, i] = -sum(row(M, i));
    }
    
    return(M);
  }
  
}

data {
  vector[61] x; // number of times allele j be counted
  int<lower=0> n; // sum of total number of times
  vector[61] pi;
  
}


parameters {
  real<lower=0> mu;
  real<lower=0> omega;
  real<lower=0> kappa;
  
}


transformed parameters {
  matrix[61,61] mutmat;
  matrix[61,61] V; //eigenvectors
  matrix[61,61] Vinv;
  vector[61] D; //eigenvalues transformed to -> 1/1-Dkk
  
  mutmat = PDRM(mu, kappa, omega);
  V = eigenvectors_sym(mutmat);
  D = inv(1-eigenvalues_sym(mutmat));
  Vinv = inv(V);
  for (i in 1:61) {
    Vinv[i] = Vinv[i] * D[i];
  }
}



model {
  matrix[61,61] alpha_Ai;
  vector[61] alpha_A;
  vector[61] lik_full;
  real lik;
  real m_Ai;
  real m_AA;
  
  // priors
  target += lognormal_lpdf( kappa | 1, 1.25 );
  target += exponential_lpdf( omega | 1 );
  target += exponential_lpdf( mu | 0.7 );
  
  // likelihood
  alpha_Ai = rep_matrix(0.,61,61);
  alpha_A = rep_vector(0.,61);
  lik_full = rep_vector(0.,61);
  
  for (A in 1:61){
    lik = 0;
  
    m_AA = dot_product(row(V, A), col(Vinv, A));
    if (m_AA < 1e-6) {
      m_AA = 1e-6;
    }
    for (i in 1:61){
      //print("m_Ai[", i, "] = ", m_Ai);
      //print("m_AA[", i, "] = ", m_AA);
      
      if (A == i){
        alpha_Ai[A,i] = 1;
      } else {
        m_Ai = dot_product(row(V, A), col(Vinv, i));
        
        if (m_Ai < 1e-6) {
          m_Ai = 1e-6;
        }

        alpha_Ai[A,i] = m_Ai/m_AA;
      }
      lik += lgamma(x[i]+alpha_Ai[A,i]) - lgamma(alpha_Ai[A,i]) - lgamma(x[i]+1);
    }
  
    alpha_A[A] = sum(alpha_Ai[A,]);
    lik_full[A] = lik + lgamma(alpha_A[A]) - lgamma(n+alpha_A[A]) + lgamma(n+1); 
  }
  target += log_sum_exp(lik_full + log(pi));
}
