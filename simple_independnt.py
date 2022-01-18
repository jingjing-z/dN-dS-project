import numpy as np
# import pandas as pd
import numpyro
import numpyro.distributions as dist
from scipy.stats import loggamma
from scipy.special import logsumexp
from numpyro.infer import HMC, HMCECS, MCMC, NUTS
from jax import random
import jax.numpy as jnp

def substitution_rate_matrix(mu, kappa, omega, pi):
    M = np.zeros((61, 61))
    M[0, 1] = kappa * mu
    M[0, 2] = omega * mu
    M[0, 3] = omega * mu
    M[0, 4] = kappa * omega * mu
    M[0, 8] = omega * mu
    M[0, 10] = omega * mu
    M[0, 13] = kappa * omega * mu
    M[0, 29] = omega * mu
    M[0, 45] = omega * mu
    M[1, 2] = omega * mu
    M[1, 3] = omega * mu
    M[1, 5] = kappa * omega * mu
    M[1, 9] = omega * mu
    M[1, 11] = omega * mu
    M[1, 14] = kappa * omega * mu
    M[1, 30] = omega * mu
    M[1, 46] = omega * mu
    M[2, 3] = kappa * mu
    M[2, 6] = kappa * omega * mu
    M[2, 15] = kappa * mu
    M[2, 31] = omega * mu
    M[2, 47] = omega * mu
    M[3, 7] = kappa * omega * mu
    M[3, 12] = omega * mu
    M[3, 16] = kappa * mu
    M[3, 32] = omega * mu
    M[3, 48] = omega * mu
    M[4, 5] = kappa * mu
    M[4, 6] = mu
    M[4, 7] = mu
    M[4, 8] = omega * mu
    M[4, 10] = omega * mu
    M[4, 17] = kappa * omega * mu
    M[4, 33] = omega * mu
    M[4, 49] = omega * mu
    M[5, 6] = mu
    M[5, 7] = mu
    M[5, 9] = omega * mu
    M[5, 11] = omega * mu
    M[5, 18] = kappa * omega * mu
    M[5, 34] = omega * mu
    M[5, 50] = omega * mu
    M[6, 7] = kappa * mu
    M[6, 19] = kappa * omega * mu
    M[6, 35] = omega * mu
    M[6, 51] = omega * mu
    M[7, 12] = omega * mu
    M[7, 20] = kappa * omega * mu
    M[7, 36] = omega * mu
    M[7, 52] = omega * mu
    M[8, 9] = kappa * mu
    M[8, 10] = kappa * omega * mu
    M[8, 21] = kappa * omega * mu
    M[8, 37] = omega * mu
    M[8, 53] = omega * mu
    M[9, 11] = kappa * omega * mu
    M[9, 22] = kappa * omega * mu
    M[9, 38] = omega * mu
    M[9, 54] = omega * mu
    M[10, 11] = kappa * mu
    M[10, 12] = omega * mu
    M[10, 25] = kappa * omega * mu
    M[10, 41] = omega * mu
    M[10, 57] = omega * mu
    M[11, 12] = omega * mu
    M[11, 26] = kappa * omega * mu
    M[11, 42] = omega * mu
    M[11, 58] = omega * mu
    M[12, 28] = kappa * omega * mu
    M[12, 44] = omega * mu
    M[12, 60] = omega * mu
    M[13, 14] = kappa * mu
    M[13, 15] = mu
    M[13, 16] = mu
    M[13, 17] = kappa * omega * mu
    M[13, 21] = omega * mu
    M[13, 25] = omega * mu
    M[13, 29] = omega * mu
    M[13, 45] = omega * mu
    M[14, 15] = mu
    M[14, 16] = mu
    M[14, 18] = kappa * omega * mu
    M[14, 22] = omega * mu
    M[14, 26] = omega * mu
    M[14, 30] = omega * mu
    M[14, 46] = omega * mu
    M[15, 16] = kappa * mu
    M[15, 19] = kappa * omega * mu
    M[15, 23] = omega * mu
    M[15, 27] = omega * mu
    M[15, 31] = omega * mu
    M[15, 47] = omega * mu
    M[16, 20] = kappa * omega * mu
    M[16, 24] = omega * mu
    M[16, 28] = omega * mu
    M[16, 32] = omega * mu
    M[16, 48] = omega * mu
    M[17, 18] = kappa * mu
    M[17, 19] = mu
    M[17, 20] = mu
    M[17, 21] = omega * mu
    M[17, 25] = omega * mu
    M[17, 33] = omega * mu
    M[17, 49] = omega * mu
    M[18, 19] = mu
    M[18, 20] = mu
    M[18, 22] = omega * mu
    M[18, 26] = omega * mu
    M[18, 34] = omega * mu
    M[18, 50] = omega * mu
    M[19, 20] = kappa * mu
    M[19, 23] = omega * mu
    M[19, 27] = omega * mu
    M[19, 35] = omega * mu
    M[19, 51] = omega * mu
    M[20, 24] = omega * mu
    M[20, 28] = omega * mu
    M[20, 36] = omega * mu
    M[20, 52] = omega * mu
    M[21, 22] = kappa * mu
    M[21, 23] = omega * mu
    M[21, 24] = omega * mu
    M[21, 25] = kappa * omega * mu
    M[21, 37] = omega * mu
    M[21, 53] = omega * mu
    M[22, 23] = omega * mu
    M[22, 24] = omega * mu
    M[22, 26] = kappa * omega * mu
    M[22, 38] = omega * mu
    M[22, 54] = omega * mu
    M[23, 24] = kappa * mu
    M[23, 27] = kappa * omega * mu
    M[23, 39] = omega * mu
    M[23, 55] = omega * mu
    M[24, 28] = kappa * omega * mu
    M[24, 40] = omega * mu
    M[24, 56] = omega * mu
    M[25, 26] = kappa * mu
    M[25, 27] = mu
    M[25, 28] = mu
    M[25, 41] = omega * mu
    M[25, 57] = omega * mu
    M[26, 27] = mu
    M[26, 28] = mu
    M[26, 42] = omega * mu
    M[26, 58] = omega * mu
    M[27, 28] = kappa * mu
    M[27, 43] = mu
    M[27, 59] = omega * mu
    M[28, 44] = mu
    M[28, 60] = omega * mu
    M[29, 30] = kappa * mu
    M[29, 31] = mu
    M[29, 32] = omega * mu
    M[29, 33] = kappa * omega * mu
    M[29, 37] = omega * mu
    M[29, 41] = omega * mu
    M[29, 45] = kappa * omega * mu
    M[30, 31] = mu
    M[30, 32] = omega * mu
    M[30, 34] = kappa * omega * mu
    M[30, 38] = omega * mu
    M[30, 42] = omega * mu
    M[30, 46] = kappa * omega * mu
    M[31, 32] = kappa * omega * mu
    M[31, 35] = kappa * omega * mu
    M[31, 39] = omega * mu
    M[31, 43] = omega * mu
    M[31, 47] = kappa * omega * mu
    M[32, 36] = kappa * omega * mu
    M[32, 40] = omega * mu
    M[32, 44] = omega * mu
    M[32, 48] = kappa * omega * mu
    M[33, 34] = kappa * mu
    M[33, 35] = mu
    M[33, 36] = mu
    M[33, 37] = omega * mu
    M[33, 41] = omega * mu
    M[33, 49] = kappa * omega * mu
    M[34, 35] = mu
    M[34, 36] = mu
    M[34, 38] = omega * mu
    M[34, 42] = omega * mu
    M[34, 50] = kappa * omega * mu
    M[35, 36] = kappa * mu
    M[35, 39] = omega * mu
    M[35, 43] = omega * mu
    M[35, 51] = kappa * omega * mu
    M[36, 40] = omega * mu
    M[36, 44] = omega * mu
    M[36, 52] = kappa * omega * mu
    M[37, 38] = kappa * mu
    M[37, 39] = omega * mu
    M[37, 40] = omega * mu
    M[37, 41] = kappa * omega * mu
    M[37, 53] = kappa * omega * mu
    M[38, 39] = omega * mu
    M[38, 40] = omega * mu
    M[38, 42] = kappa * omega * mu
    M[38, 54] = kappa * omega * mu
    M[39, 40] = kappa * mu
    M[39, 43] = kappa * omega * mu
    M[39, 55] = kappa * omega * mu
    M[40, 44] = kappa * omega * mu
    M[40, 56] = kappa * omega * mu
    M[41, 42] = kappa * mu
    M[41, 43] = omega * mu
    M[41, 44] = omega * mu
    M[41, 57] = kappa * omega * mu
    M[42, 43] = omega * mu
    M[42, 44] = omega * mu
    M[42, 58] = kappa * omega * mu
    M[43, 44] = kappa * mu
    M[43, 59] = kappa * omega * mu
    M[44, 60] = kappa * omega * mu
    M[45, 46] = kappa * mu
    M[45, 47] = mu
    M[45, 48] = mu
    M[45, 49] = kappa * omega * mu
    M[45, 53] = omega * mu
    M[45, 57] = omega * mu
    M[46, 47] = mu
    M[46, 48] = mu
    M[46, 50] = kappa * omega * mu
    M[46, 54] = omega * mu
    M[46, 58] = omega * mu
    M[47, 48] = kappa * mu
    M[47, 51] = kappa * omega * mu
    M[47, 55] = omega * mu
    M[47, 59] = omega * mu
    M[48, 52] = kappa * omega * mu
    M[48, 56] = omega * mu
    M[48, 60] = omega * mu
    M[49, 50] = kappa * mu
    M[49, 51] = mu
    M[49, 52] = mu
    M[49, 53] = omega * mu
    M[49, 57] = omega * mu
    M[50, 51] = mu
    M[50, 52] = mu
    M[50, 54] = omega * mu
    M[50, 58] = omega * mu
    M[51, 52] = kappa * mu
    M[51, 55] = omega * mu
    M[51, 59] = omega * mu
    M[52, 56] = omega * mu
    M[52, 60] = omega * mu
    M[53, 54] = kappa * mu
    M[53, 55] = omega * mu
    M[53, 56] = omega * mu
    M[53, 57] = kappa * omega * mu
    M[54, 55] = omega * mu
    M[54, 56] = omega * mu
    M[54, 58] = kappa * omega * mu
    M[55, 56] = kappa * mu
    M[55, 59] = kappa * omega * mu
    M[56, 60] = kappa * omega * mu
    M[57, 58] = kappa * mu
    M[57, 59] = mu
    M[57, 60] = mu
    M[58, 59] = mu
    M[58, 60] = mu
    M[59, 60] = kappa * mu

    for i in range(61):
        M[i] *= pi

    # Fill in the lower triangle
    for i in range(61):
        for j in range(61):
            M[j, i] = M[i, j]

    # Compute the diagonal
    rowsum = np.sum(M, axis=1)
    for i in range(61):
        M[i, i] = -rowsum[i]

    return M


def likelihood(X, n, mu, kappa, omega, pi):
    D = [0] * 61
    alpha_Ai = np.zeros((61, 61))
    alpha_A = [0] * 61
    lik_full = [0] * 61

    mutmat = substitution_rate_matrix(mu, kappa, omega, pi)
    a, b = np.linalg.eig(mutmat) # eigenvalues and eigenvectors
    V = b
    D = [1 / (1-A) for A in a.tolist()]
    VD = V
    for i in range(61):
        VD[i] *= D

    for A in range(61):
        lik = 0
        m_AA = np.dot(V[A, :], VD[A, :])

        if m_AA < 1e-6:
            m_AA = 1e-6

        for i in range(61):
            if A==i:
                alpha_Ai[A, i] = 1
            else:
                m_Ai = np.dot(V[A, :], VD[i, :])

                if m_Ai < 1e-6:
                    m_Ai = 1e-6
                alpha_Ai[A, i] = m_Ai / m_AA
            lik += loggamma.pdf(X[i] + alpha_Ai[A, i]) - loggamma.pdf(alpha_Ai[A, i]) -\
                loggamma(X[i] + 1)

        alpha_A[A] = np.sum(alpha_Ai[A, :])
        lik_full[A] = lik + loggamma.pdf(alpha_A[A]) - loggamma.pdf(n + alpha_A[A]) +\
            loggamma(n + 1)

    partial_ll += logsumexp(lik_full + np.log(pi))
    return partial_ll



def model(X, n, pi):
    # parameters
    kappa_prior = dist.LogNormal(1, 1.5)
    omega_prior = dist.Exponential(0.5)
    mu_prior = dist.Exponential(0.5)
    kappa = numpyro.sample("kappa", kappa_prior)
    omega = numpyro.sample("omega", omega_prior)
    mu = numpyro.sample("mu", mu_prior)
    with numpyro.plate('n', n):
        numpyro.sample('obs', likelihood(X, n, kappa, omega, mu, pi))

X = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 166, 4723, 0, 0]
n = 4889
pi = [1 / 61] * 61

nuts_kernel = NUTS(model)
mcmc = MCMC(nuts_kernel, num_warmup=500, num_samples=1000)
rng_key = random.PRNGKey(0)
mcmc.run(rng_key, X, n, pi, extra_fields=('potential_energy',))