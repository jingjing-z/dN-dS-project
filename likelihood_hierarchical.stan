functions {
  // substitution rate matrix
	matrix PDRM(real mu, real kappa, real omega, vector pi) {
	  matrix[61,61] M;
	  row_vector[61] equilibrium;
	  
	  equilibrium = to_row_vector(pi);
	  M = rep_matrix(0.,61,61);
	  
    M[1,2] = kappa*mu;
    M[1,3] = omega*mu;
    M[1,4] = omega*mu;
    M[1,5] = kappa*omega*mu;
    M[1,9] = omega*mu;
    M[1,11] = omega*mu;
    M[1,14] = kappa*omega*mu;
    M[1,30] = omega*mu;
    M[1,46] = omega*mu;
    M[2,3] = omega*mu;
    M[2,4] = omega*mu;
    M[2,6] = kappa*omega*mu;
    M[2,10] = omega*mu;
    M[2,12] = omega*mu;
    M[2,15] = kappa*omega*mu;
    M[2,31] = omega*mu;
    M[2,47] = omega*mu;
    M[3,4] = kappa*mu;
    M[3,7] = kappa*omega*mu;
    M[3,16] = kappa*mu;
    M[3,32] = omega*mu;
    M[3,48] = omega*mu;
    M[4,8] = kappa*omega*mu;
    M[4,13] = omega*mu;
    M[4,17] = kappa*mu;
    M[4,33] = omega*mu;
    M[4,49] = omega*mu;
    M[5,6] = kappa*mu;
    M[5,7] = mu;
    M[5,8] = mu;
    M[5,9] = omega*mu;
    M[5,11] = omega*mu;
    M[5,18] = kappa*omega*mu;
    M[5,34] = omega*mu;
    M[5,50] = omega*mu;
    M[6,7] = mu;
    M[6,8] = mu;
    M[6,10] = omega*mu;
    M[6,12] = omega*mu;
    M[6,19] = kappa*omega*mu;
    M[6,35] = omega*mu;
    M[6,51] = omega*mu;
    M[7,8] = kappa*mu;
    M[7,20] = kappa*omega*mu;
    M[7,36] = omega*mu;
    M[7,52] = omega*mu;
    M[8,13] = omega*mu;
    M[8,21] = kappa*omega*mu;
    M[8,37] = omega*mu;
    M[8,53] = omega*mu;
    M[9,10] = kappa*mu;
    M[9,11] = kappa*omega*mu;
    M[9,22] = kappa*omega*mu;
    M[9,38] = omega*mu;
    M[9,54] = omega*mu;
    M[10,12] = kappa*omega*mu;
    M[10,23] = kappa*omega*mu;
    M[10,39] = omega*mu;
    M[10,55] = omega*mu;
    M[11,12] = kappa*mu;
    M[11,13] = omega*mu;
    M[11,26] = kappa*omega*mu;
    M[11,42] = omega*mu;
    M[11,58] = omega*mu;
    M[12,13] = omega*mu;
    M[12,27] = kappa*omega*mu;
    M[12,43] = omega*mu;
    M[12,59] = omega*mu;
    M[13,29] = kappa*omega*mu;
    M[13,45] = omega*mu;
    M[13,61] = omega*mu;
    M[14,15] = kappa*mu;
    M[14,16] = mu;
    M[14,17] = mu;
    M[14,18] = kappa*omega*mu;
    M[14,22] = omega*mu;
    M[14,26] = omega*mu;
    M[14,30] = omega*mu;
    M[14,46] = omega*mu;
    M[15,16] = mu;
    M[15,17] = mu;
    M[15,19] = kappa*omega*mu;
    M[15,23] = omega*mu;
    M[15,27] = omega*mu;
    M[15,31] = omega*mu;
    M[15,47] = omega*mu;
    M[16,17] = kappa*mu;
    M[16,20] = kappa*omega*mu;
    M[16,24] = omega*mu;
    M[16,28] = omega*mu;
    M[16,32] = omega*mu;
    M[16,48] = omega*mu;
    M[17,21] = kappa*omega*mu;
    M[17,25] = omega*mu;
    M[17,29] = omega*mu;
    M[17,33] = omega*mu;
    M[17,49] = omega*mu;
    M[18,19] = kappa*mu;
    M[18,20] = mu;
    M[18,21] = mu;
    M[18,22] = omega*mu;
    M[18,26] = omega*mu;
    M[18,34] = omega*mu;
    M[18,50] = omega*mu;
    M[19,20] = mu;
    M[19,21] = mu;
    M[19,23] = omega*mu;
    M[19,27] = omega*mu;
    M[19,35] = omega*mu;
    M[19,51] = omega*mu;
    M[20,21] = kappa*mu;
    M[20,24] = omega*mu;
    M[20,28] = omega*mu;
    M[20,36] = omega*mu;
    M[20,52] = omega*mu;
    M[21,25] = omega*mu;
    M[21,29] = omega*mu;
    M[21,37] = omega*mu;
    M[21,53] = omega*mu;
    M[22,23] = kappa*mu;
    M[22,24] = omega*mu;
    M[22,25] = omega*mu;
    M[22,26] = kappa*omega*mu;
    M[22,38] = omega*mu;
    M[22,54] = omega*mu;
    M[23,24] = omega*mu;
    M[23,25] = omega*mu;
    M[23,27] = kappa*omega*mu;
    M[23,39] = omega*mu;
    M[23,55] = omega*mu;
    M[24,25] = kappa*mu;
    M[24,28] = kappa*omega*mu;
    M[24,40] = omega*mu;
    M[24,56] = omega*mu;
    M[25,29] = kappa*omega*mu;
    M[25,41] = omega*mu;
    M[25,57] = omega*mu;
    M[26,27] = kappa*mu;
    M[26,28] = mu;
    M[26,29] = mu;
    M[26,42] = omega*mu;
    M[26,58] = omega*mu;
    M[27,28] = mu;
    M[27,29] = mu;
    M[27,43] = omega*mu;
    M[27,59] = omega*mu;
    M[28,29] = kappa*mu;
    M[28,44] = mu;
    M[28,60] = omega*mu;
    M[29,45] = mu;
    M[29,61] = omega*mu;
    M[30,31] = kappa*mu;
    M[30,32] = mu;
    M[30,33] = omega*mu;
    M[30,34] = kappa*omega*mu;
    M[30,38] = omega*mu;
    M[30,42] = omega*mu;
    M[30,46] = kappa*omega*mu;
    M[31,32] = mu;
    M[31,33] = omega*mu;
    M[31,35] = kappa*omega*mu;
    M[31,39] = omega*mu;
    M[31,43] = omega*mu;
    M[31,47] = kappa*omega*mu;
    M[32,33] = kappa*omega*mu;
    M[32,36] = kappa*omega*mu;
    M[32,40] = omega*mu;
    M[32,44] = omega*mu;
    M[32,48] = kappa*omega*mu;
    M[33,37] = kappa*omega*mu;
    M[33,41] = omega*mu;
    M[33,45] = omega*mu;
    M[33,49] = kappa*omega*mu;
    M[34,35] = kappa*mu;
    M[34,36] = mu;
    M[34,37] = mu;
    M[34,38] = omega*mu;
    M[34,42] = omega*mu;
    M[34,50] = kappa*omega*mu;
    M[35,36] = mu;
    M[35,37] = mu;
    M[35,39] = omega*mu;
    M[35,43] = omega*mu;
    M[35,51] = kappa*omega*mu;
    M[36,37] = kappa*mu;
    M[36,40] = omega*mu;
    M[36,44] = omega*mu;
    M[36,52] = kappa*omega*mu;
    M[37,41] = omega*mu;
    M[37,45] = omega*mu;
    M[37,53] = kappa*omega*mu;
    M[38,39] = kappa*mu;
    M[38,40] = omega*mu;
    M[38,41] = omega*mu;
    M[38,42] = kappa*omega*mu;
    M[38,54] = kappa*omega*mu;
    M[39,40] = omega*mu;
    M[39,41] = omega*mu;
    M[39,43] = kappa*omega*mu;
    M[39,55] = kappa*omega*mu;
    M[40,41] = kappa*mu;
    M[40,44] = kappa*omega*mu;
    M[40,56] = kappa*omega*mu;
    M[41,45] = kappa*omega*mu;
    M[41,57] = kappa*omega*mu;
    M[42,43] = kappa*mu;
    M[42,44] = omega*mu;
    M[42,45] = omega*mu;
    M[42,58] = kappa*omega*mu;
    M[43,44] = omega*mu;
    M[43,45] = omega*mu;
    M[43,59] = kappa*omega*mu;
    M[44,45] = kappa*mu;
    M[44,60] = kappa*omega*mu;
    M[45,61] = kappa*omega*mu;
    M[46,47] = kappa*mu;
    M[46,48] = mu;
    M[46,49] = mu;
    M[46,50] = kappa*omega*mu;
    M[46,54] = omega*mu;
    M[46,58] = omega*mu;
    M[47,48] = mu;
    M[47,49] = mu;
    M[47,51] = kappa*omega*mu;
    M[47,55] = omega*mu;
    M[47,59] = omega*mu;
    M[48,49] = kappa*mu;
    M[48,52] = kappa*omega*mu;
    M[48,56] = omega*mu;
    M[48,60] = omega*mu;
    M[49,53] = kappa*omega*mu;
    M[49,57] = omega*mu;
    M[49,61] = omega*mu;
    M[50,51] = kappa*mu;
    M[50,52] = mu;
    M[50,53] = mu;
    M[50,54] = omega*mu;
    M[50,58] = omega*mu;
    M[51,52] = mu;
    M[51,53] = mu;
    M[51,55] = omega*mu;
    M[51,59] = omega*mu;
    M[52,53] = kappa*mu;
    M[52,56] = omega*mu;
    M[52,60] = omega*mu;
    M[53,57] = omega*mu;
    M[53,61] = omega*mu;
    M[54,55] = kappa*mu;
    M[54,56] = omega*mu;
    M[54,57] = omega*mu;
    M[54,58] = kappa*omega*mu;
    M[55,56] = omega*mu;
    M[55,57] = omega*mu;
    M[55,59] = kappa*omega*mu;
    M[56,57] = kappa*mu;
    M[56,60] = kappa*omega*mu;
    M[57,61] = kappa*omega*mu;
    M[58,59] = kappa*mu;
    M[58,60] = mu;
    M[58,61] = mu;
    M[59,60] = mu;
    M[59,61] = mu;
    M[60,61] = kappa*mu;
    
    // Fill in the lower triangle
    //M = M'+ M;
    for (i in 1:61) {
      for (j in (i+1):61) {
        M[j, i] = M[i, j];
      }
    }
    
    for (i in 1:61) {
      M[i] .*= equilibrium; 
    }
  
    // Compute the diagonal
    for (i in 1:61){
      M[i, i] = -sum(row(M, i));
    }
    
    return(M);
  }
  
}

data {
  int<lower=0> l; // length of the gene
  matrix[l,61] X; // number of times allele j be counted
  vector[l] n; // sum of total number of times
  vector[61] pi;
  
}


parameters {
  vector<lower=0>[l] mu;
  vector<lower=0>[l] omega;
  vector<lower=0>[l] kappa;
  
  vector<lower=0>[l] omega_popsd;
  real<lower=0> log_kappa_popmean;
  real<lower=0> omega_mean;
  vector<lower=0>[l] log_kappa_popsd;
}



model {
  matrix[61,61] mutmat;
  matrix[61,61] V; //eigenvectors
  matrix[61,61] VD;
  row_vector[61] D; //eigenvalues transformed to -> 1/1-Dkk
  
  matrix[61,61] alpha_Ai;
  vector[61] alpha_A;
  vector[61] lik_full;
  
  real lik;
  real m_Ai;
  real m_AA;
  
  
  // make the model fully hierarchical
    for (H in 1:l){
    // priors
    target += normal_lpdf( log_kappa_popmean | 0, 1 );
    target += exponential_lpdf( log_kappa_popsd | 10 );
    target += normal_lpdf( log_kappa_rnde[H] | 0, log_kappa_popsd);
    
    target += normal_lpdf( omega[H] | 0, omega_popsd );
    target += cauchy_lpdf( omega_popsd | 0, 1)
    
    target += exponential_lpdf( mu | 0.7 );
  
    // likelihood
    alpha_Ai = rep_matrix(0.,61,61);
    alpha_A = rep_vector(0.,61);
    lik_full = rep_vector(0.,61);
    
    // transforms
    kappa[H] = exp( log_kappa_popmean + log_kappa_rnde[H])
    
    mutmat = PDRM(mu[H], kappa[H], omega[H], pi);
    V = eigenvectors_sym(mutmat);
    D = to_row_vector(inv(1-eigenvalues_sym(mutmat)));
    VD = V;
    for (i in 1:61) {
    VD[i] .*= D;
    }
    
    for (A in 1:61){
    lik = 0;
  
    m_AA = dot_product(row(V, A), row(VD, A));
    if (m_AA < 1e-6) {
      m_AA = 1e-6;
    }
    //print("m_AA[", i, "] = ", m_AA);
    for (i in 1:61){
      
      if (A == i){
        alpha_Ai[A,i] = 1;
      } else {
        m_Ai = dot_product(row(V, A), row(VD, i));
        
        if (m_Ai < 1e-6) {
          m_Ai = 1e-6;
        }
        //print("m_Ai[", i, "] = ", m_Ai);

        alpha_Ai[A,i] = m_Ai/m_AA;
      }
      lik += lgamma(X[H,i]+alpha_Ai[A,i]) - lgamma(alpha_Ai[A,i]) - lgamma(X[H,i]+1);
    }
  
    alpha_A[A] = sum(alpha_Ai[A,]);
    lik_full[A] = lik + lgamma(alpha_A[A]) - lgamma(n[H]+alpha_A[A]) + lgamma(n[H]+1); 
  }
  target += log_sum_exp(lik_full + log(pi));
 }
    
}
  
