functions {
  // substitution rate matrix
  matrix GTR(real mu, real omega, real AG, real AC, real AT, real GC, real CT, vector pi) {
    matrix[61,61] M;
    row_vector[61] equilibrium;
      
    equilibrium = to_row_vector(pi);
    M = rep_matrix(0.,61,61);
      
    M[1,2] = CT*mu;
    M[1,3] = AT*omega*mu;
    M[1,4] = omega*mu;
    M[1,5] = CT*omega*mu;
    M[1,9] = AT*omega*mu;
    M[1,11] = omega*mu;
    M[1,14] = CT*omega*mu;
    M[1,30] = AT*omega*mu;
    M[1,46] = omega*mu;
    M[2,3] = AC*omega*mu;
    M[2,4] = GC*omega*mu;
    M[2,6] = CT*omega*mu;
    M[2,10] = AT*omega*mu;
    M[2,12] = omega*mu;
    M[2,15] = CT*omega*mu;
    M[2,31] = AT*omega*mu;
    M[2,47] = omega*mu;
    M[3,4] = AG*mu;
    M[3,7] = CT*omega*mu;
    M[3,16] = CT*mu;
    M[3,32] = AT*omega*mu;
    M[3,48] = omega*mu;
    M[4,8] = CT*omega*mu;
    M[4,13] = omega*mu;
    M[4,17] = CT*mu;
    M[4,33] = AT*omega*mu;
    M[4,49] = omega*mu;
    M[5,6] = CT*mu;
    M[5,7] = AT*mu;
    M[5,8] = mu;
    M[5,9] = AC*omega*mu;
    M[5,11] = GC*omega*mu;
    M[5,18] = CT*omega*mu;
    M[5,34] = AT*omega*mu;
    M[5,50] = omega*mu;
    M[6,7] = AC*mu;
    M[6,8] = GC*mu;
    M[6,10] = AC*omega*mu;
    M[6,12] = GC*omega*mu;
    M[6,19] = CT*omega*mu;
    M[6,35] = AT*omega*mu;
    M[6,51] = omega*mu;
    M[7,8] = AG*mu;
    M[7,20] = CT*omega*mu;
    M[7,36] = AT*omega*mu;
    M[7,52] = omega*mu;
    M[8,13] = GC*omega*mu;
    M[8,21] = CT*omega*mu;
    M[8,37] = AT*omega*mu;
    M[8,53] = omega*mu;
    M[9,10] = CT*mu;
    M[9,11] = AG*omega*mu;
    M[9,22] = CT*omega*mu;
    M[9,38] = AT*omega*mu;
    M[9,54] = omega*mu;
    M[10,12] = AG*omega*mu;
    M[10,23] = CT*omega*mu;
    M[10,39] = AT*omega*mu;
    M[10,55] = omega*mu;
    M[11,12] = CT*mu;
    M[11,13] = omega*mu;
    M[11,26] = CT*omega*mu;
    M[11,42] = AT*omega*mu;
    M[11,58] = omega*mu;
    M[12,13] = GC*omega*mu;
    M[12,27] = CT*omega*mu;
    M[12,43] = AT*omega*mu;
    M[12,59] = omega*mu;
    M[13,29] = CT*omega*mu;
    M[13,45] = AT*omega*mu;
    M[13,61] = omega*mu;
    M[14,15] = CT*mu;
    M[14,16] = AT*mu;
    M[14,17] = mu;
    M[14,18] = CT*omega*mu;
    M[14,22] = AT*omega*mu;
    M[14,26] = omega*mu;
    M[14,30] = AC*omega*mu;
    M[14,46] = GC*omega*mu;
    M[15,16] = AC*mu;
    M[15,17] = GC*mu;
    M[15,19] = CT*omega*mu;
    M[15,23] = AT*omega*mu;
    M[15,27] = omega*mu;
    M[15,31] = AC*omega*mu;
    M[15,47] = GC*omega*mu;
    M[16,17] = AG*mu;
    M[16,20] = CT*omega*mu;
    M[16,24] = AT*omega*mu;
    M[16,28] = omega*mu;
    M[16,32] = AC*omega*mu;
    M[16,48] = GC*omega*mu;
    M[17,21] = CT*omega*mu;
    M[17,25] = AT*omega*mu;
    M[17,29] = omega*mu;
    M[17,33] = AC*omega*mu;
    M[17,49] = GC*omega*mu;
    M[18,19] = CT*mu;
    M[18,20] = AT*mu;
    M[18,21] = mu;
    M[18,22] = AC*omega*mu;
    M[18,26] = GC*omega*mu;
    M[18,34] = AC*omega*mu;
    M[18,50] = GC*omega*mu;
    M[19,20] = AC*mu;
    M[19,21] = GC*mu;
    M[19,23] = AC*omega*mu;
    M[19,27] = GC*omega*mu;
    M[19,35] = AC*omega*mu;
    M[19,51] = GC*omega*mu;
    M[20,21] = AG*mu;
    M[20,24] = AC*omega*mu;
    M[20,28] = GC*omega*mu;
    M[20,36] = AC*omega*mu;
    M[20,52] = GC*omega*mu;
    M[21,25] = AC*omega*mu;
    M[21,29] = GC*omega*mu;
    M[21,37] = AC*omega*mu;
    M[21,53] = GC*omega*mu;
    M[22,23] = CT*mu;
    M[22,24] = AT*omega*mu;
    M[22,25] = omega*mu;
    M[22,26] = AG*omega*mu;
    M[22,38] = AC*omega*mu;
    M[22,54] = GC*omega*mu;
    M[23,24] = AC*omega*mu;
    M[23,25] = GC*omega*mu;
    M[23,27] = AG*omega*mu;
    M[23,39] = AC*omega*mu;
    M[23,55] = GC*omega*mu;
    M[24,25] = AG*mu;
    M[24,28] = AG*omega*mu;
    M[24,40] = AC*omega*mu;
    M[24,56] = GC*omega*mu;
    M[25,29] = AG*omega*mu;
    M[25,41] = AC*omega*mu;
    M[25,57] = GC*omega*mu;
    M[26,27] = CT*mu;
    M[26,28] = AT*mu;
    M[26,29] = mu;
    M[26,42] = AC*omega*mu;
    M[26,58] = GC*omega*mu;
    M[27,28] = AC*mu;
    M[27,29] = GC*mu;
    M[27,43] = AC*omega*mu;
    M[27,59] = GC*omega*mu;
    M[28,29] = AG*mu;
    M[28,44] = AC*mu;
    M[28,60] = GC*omega*mu;
    M[29,45] = AC*mu;
    M[29,61] = GC*omega*mu;
    M[30,31] = CT*mu;
    M[30,32] = AT*mu;
    M[30,33] = omega*mu;
    M[30,34] = CT*omega*mu;
    M[30,38] = AT*omega*mu;
    M[30,42] = omega*mu;
    M[30,46] = AG*omega*mu;
    M[31,32] = AC*mu;
    M[31,33] = GC*omega*mu;
    M[31,35] = CT*omega*mu;
    M[31,39] = AT*omega*mu;
    M[31,43] = omega*mu;
    M[31,47] = AG*omega*mu;
    M[32,33] = AG*omega*mu;
    M[32,36] = CT*omega*mu;
    M[32,40] = AT*omega*mu;
    M[32,44] = omega*mu;
    M[32,48] = AG*omega*mu;
    M[33,37] = CT*omega*mu;
    M[33,41] = AT*omega*mu;
    M[33,45] = omega*mu;
    M[33,49] = AG*omega*mu;
    M[34,35] = CT*mu;
    M[34,36] = AT*mu;
    M[34,37] = mu;
    M[34,38] = AC*omega*mu;
    M[34,42] = GC*omega*mu;
    M[34,50] = AG*omega*mu;
    M[35,36] = AC*mu;
    M[35,37] = GC*mu;
    M[35,39] = AC*omega*mu;
    M[35,43] = GC*omega*mu;
    M[35,51] = AG*omega*mu;
    M[36,37] = AG*mu;
    M[36,40] = AC*omega*mu;
    M[36,44] = GC*omega*mu;
    M[36,52] = AG*omega*mu;
    M[37,41] = AC*omega*mu;
    M[37,45] = GC*omega*mu;
    M[37,53] = AG*omega*mu;
    M[38,39] = CT*mu;
    M[38,40] = AT*omega*mu;
    M[38,41] = omega*mu;
    M[38,42] = AG*omega*mu;
    M[38,54] = AG*omega*mu;
    M[39,40] = AC*omega*mu;
    M[39,41] = GC*omega*mu;
    M[39,43] = AG*omega*mu;
    M[39,55] = AG*omega*mu;
    M[40,41] = AG*mu;
    M[40,44] = AG*omega*mu;
    M[40,56] = AG*omega*mu;
    M[41,45] = AG*omega*mu;
    M[41,57] = AG*omega*mu;
    M[42,43] = CT*mu;
    M[42,44] = AT*omega*mu;
    M[42,45] = omega*mu;
    M[42,58] = AG*omega*mu;
    M[43,44] = AC*omega*mu;
    M[43,45] = GC*omega*mu;
    M[43,59] = AG*omega*mu;
    M[44,45] = AG*mu;
    M[44,60] = AG*omega*mu;
    M[45,61] = AG*omega*mu;
    M[46,47] = CT*mu;
    M[46,48] = AT*mu;
    M[46,49] = mu;
    M[46,50] = CT*omega*mu;
    M[46,54] = AT*omega*mu;
    M[46,58] = omega*mu;
    M[47,48] = AC*mu;
    M[47,49] = GC*mu;
    M[47,51] = CT*omega*mu;
    M[47,55] = AT*omega*mu;
    M[47,59] = omega*mu;
    M[48,49] = AG*mu;
    M[48,52] = CT*omega*mu;
    M[48,56] = AT*omega*mu;
    M[48,60] = omega*mu;
    M[49,53] = CT*omega*mu;
    M[49,57] = AT*omega*mu;
    M[49,61] = omega*mu;
    M[50,51] = CT*mu;
    M[50,52] = AT*mu;
    M[50,53] = mu;
    M[50,54] = AC*omega*mu;
    M[50,58] = GC*omega*mu;
    M[51,52] = AC*mu;
    M[51,53] = GC*mu;
    M[51,55] = AC*omega*mu;
    M[51,59] = GC*omega*mu;
    M[52,53] = AG*mu;
    M[52,56] = AC*omega*mu;
    M[52,60] = GC*omega*mu;
    M[53,57] = AC*omega*mu;
    M[53,61] = GC*omega*mu;
    M[54,55] = CT*mu;
    M[54,56] = AT*omega*mu;
    M[54,57] = omega*mu;
    M[54,58] = AG*omega*mu;
    M[55,56] = AC*omega*mu;
    M[55,57] = GC*omega*mu;
    M[55,59] = AG*omega*mu;
    M[56,57] = AG*mu;
    M[56,60] = AG*omega*mu;
    M[57,61] = AG*omega*mu;
    M[58,59] = CT*mu;
    M[58,60] = AT*mu;
    M[58,61] = mu;
    M[59,60] = AC*mu;
    M[59,61] = GC*mu;
    M[60,61] = AG*mu;
    
    for (i in 1:61) {
      M[i] .*= equilibrium; 
    }
    
    // Fill in the lower triangle
    //M = M'+ M;
    for (i in 1:61) {
      for (j in (i+1):61) {
        M[j, i] = M[i, j];
      }
    }
    
    // Compute the diagonal
    for (i in 1:61){
      M[i, i] = -sum(row(M, i));
    }
    
    return(M);
  }
  
  
  // parallelization
  real partial_sum(int[] n_slice,
                   int start, int end,
                   vector AG,
                   vector AC,
                   vector AT,
                   vector GC,
                   vector CT,
                   vector omega,
                   vector mu,
                   vector pi,
                   matrix X) {
                     
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
  
    real partial_ll;
    partial_ll = 0;
    
    for (H in start:end){
      
      // improper
      
      //partial_ll += lognormal_lpdf( beta | 1, 0.5 );
      //partial_ll += lognormal_lpdf( gamma | 1.5, 0.5 );
      //partial_ll += lognormal_lpdf( delta | 2, 0.8 );
      //partial_ll += lognormal_lpdf( epsilon | 4, 2 );
      //partial_ll += lognormal_lpdf( eta | 6, 2.5 );
      
      partial_ll += exponential_lpdf( omega[H] | 1 );
      partial_ll += exponential_lpdf( mu[H] | 1 );  
           
      // transforms 
      print("AG = ", AG);
      print("AC = ", AC);
      print("AT = ", AT);
      print("GC = ", GC);
      print("CT = ", CT);
      
      print("omega = ", omega);
      print("mu = ", mu);
      
      mutmat = GTR(mu[H], omega[H], AG[H], AC[H], AT[H], GC[H], CT[H], pi);
      print("mutmat[", H, "] = ", mutmat);
      
      V = eigenvectors_sym(mutmat); // eigenvectors (eigen decomposition)
      D = to_row_vector(inv(1-eigenvalues_sym(mutmat))); // 1/(1-D) in equation (5)
      
      print("V[", H, "] = ", V);
      print("D[", H, "] = ", D);
        
      VD = V;
      for (i in 1:61) {
        VD[i] .*= D;
        
      }
           
      // likelihood
      alpha_Ai = rep_matrix(0.,61,61);
      alpha_A = rep_vector(0.,61);
      lik_full = rep_vector(0.,61);
    
      for (A in 1:61){
        lik = 0;
        m_AA = dot_product(row(V, A), row(VD, A));
             
        if (m_AA < 1e-6) {
          m_AA = 1e-6;
        }
       
        for (i in 1:61){
          if (A == i){
            alpha_Ai[A,i] = 1;
          } else {
            m_Ai = dot_product(row(V, A), row(VD, i));
            if (m_Ai < 1e-6) {
              m_Ai = 1e-6;
            }
            alpha_Ai[A,i] = m_Ai/m_AA;
          }
                 
          lik += lgamma(X[H,i]+alpha_Ai[A,i]) - lgamma(alpha_Ai[A,i]) - lgamma(X[H,i]+1);
        }
  
        alpha_A[A] = sum(alpha_Ai[A,]);
        lik_full[A] = lik + lgamma(alpha_A[A]) - lgamma(n_slice[H-start]+alpha_A[A]) + lgamma(n_slice[H-start]+1);
           
      }
      
      partial_ll += log_sum_exp(lik_full + log(pi));
    }
    
    return partial_ll;
  }
}


data {
  int<lower=0> l; // length of the gene, say 294
  matrix[l,61] X; // number of times allele j be counted
  int n[l]; // sum of total number of times
  vector[61] pi;
}


parameters {
  vector<lower=0>[l] mu;
  
  vector<lower=0>[l] AG;
  vector<lower=0>[l] AC;
  vector<lower=0>[l] AT;
  vector<lower=0>[l] GC;
  vector<lower=0>[l] CT;
  
  vector<lower=0>[l] omega;
}


model{
  int grainsize = 1;

  target += reduce_sum(partial_sum, n, 
                       grainsize,
                       AG, AC, AT, GC, CT,
                       omega, mu, pi, X);
                       
}