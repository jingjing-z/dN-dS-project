import argparse
import os

import jax
import jax.numpy as jnp
import jax.random as random
from jax.scipy.special import logsumexp
from jax.lax import lgamma

import numpyro
import numpyro.distributions as dist
from numpyro.infer import HMC, MCMC, NUTS, SA, Predictive, log_likelihood

def substitution_rate_matrix(mu, kappa, omega, pi):
    M = jnp.zeros((61, 61))
    M = M.at[0, 1].set(kappa * mu)
    M = M.at[0, 2].set(omega * mu)
    M = M.at[0, 3].set(omega * mu)
    M = M.at[0, 4].set(kappa * omega * mu)
    M = M.at[0, 8].set(omega * mu)
    M = M.at[0, 10].set(omega * mu)
    M = M.at[0, 13].set(kappa * omega * mu)
    M = M.at[0, 29].set(omega * mu)
    M = M.at[0, 45].set(omega * mu)
    M = M.at[1, 2].set(omega * mu)
    M = M.at[1, 3].set(omega * mu)
    M = M.at[1, 5].set(kappa * omega * mu)
    M = M.at[1, 9].set(omega * mu)
    M = M.at[1, 11].set(omega * mu)
    M = M.at[1, 14].set(kappa * omega * mu)
    M = M.at[1, 30].set(omega * mu)
    M = M.at[1, 46].set(omega * mu)
    M = M.at[2, 3].set(kappa * mu)
    M = M.at[2, 6].set(kappa * omega * mu)
    M = M.at[2, 15].set(kappa * mu)
    M = M.at[2, 31].set(omega * mu)
    M = M.at[2, 47].set(omega * mu)
    M = M.at[3, 7].set(kappa * omega * mu)
    M = M.at[3, 12].set(omega * mu)
    M = M.at[3, 16].set(kappa * mu)
    M = M.at[3, 32].set(omega * mu)
    M = M.at[3, 48].set(omega * mu)
    M = M.at[4, 5].set(kappa * mu)
    M = M.at[4, 6].set(mu)
    M = M.at[4, 7].set(mu)
    M = M.at[4, 8].set(omega * mu)
    M = M.at[4, 10].set(omega * mu)
    M = M.at[4, 17].set(kappa * omega * mu)
    M = M.at[4, 33].set(omega * mu)
    M = M.at[4, 49].set(omega * mu)
    M = M.at[5, 6].set(mu)
    M = M.at[5, 7].set(mu)
    M = M.at[5, 9].set(omega * mu)
    M = M.at[5, 11].set(omega * mu)
    M = M.at[5, 18].set(kappa * omega * mu)
    M = M.at[5, 34].set(omega * mu)
    M = M.at[5, 50].set(omega * mu)
    M = M.at[6, 7].set(kappa * mu)
    M = M.at[6, 19].set(kappa * omega * mu)
    M = M.at[6, 35].set(omega * mu)
    M = M.at[6, 51].set(omega * mu)
    M = M.at[7, 12].set(omega * mu)
    M = M.at[7, 20].set(kappa * omega * mu)
    M = M.at[7, 36].set(omega * mu)
    M = M.at[7, 52].set(omega * mu)
    M = M.at[8, 9].set(kappa * mu)
    M = M.at[8, 10].set(kappa * omega * mu)
    M = M.at[8, 21].set(kappa * omega * mu)
    M = M.at[8, 37].set(omega * mu)
    M = M.at[8, 53].set(omega * mu)
    M = M.at[9, 11].set(kappa * omega * mu)
    M = M.at[9, 22].set(kappa * omega * mu)
    M = M.at[9, 38].set(omega * mu)
    M = M.at[9, 54].set(omega * mu)
    M = M.at[10, 11].set(kappa * mu)
    M = M.at[10, 12].set(omega * mu)
    M = M.at[10, 25].set(kappa * omega * mu)
    M = M.at[10, 41].set(omega * mu)
    M = M.at[10, 57].set(omega * mu)
    M = M.at[11, 12].set(omega * mu)
    M = M.at[11, 26].set(kappa * omega * mu)
    M = M.at[11, 42].set(omega * mu)
    M = M.at[11, 58].set(omega * mu)
    M = M.at[12, 28].set(kappa * omega * mu)
    M = M.at[12, 44].set(omega * mu)
    M = M.at[12, 60].set(omega * mu)
    M = M.at[13, 14].set(kappa * mu)
    M = M.at[13, 15].set(mu)
    M = M.at[13, 16].set(mu)
    M = M.at[13, 17].set(kappa * omega * mu)
    M = M.at[13, 21].set(omega * mu)
    M = M.at[13, 25].set(omega * mu)
    M = M.at[13, 29].set(omega * mu)
    M = M.at[13, 45].set(omega * mu)
    M = M.at[14, 15].set(mu)
    M = M.at[14, 16].set(mu)
    M = M.at[14, 18].set(kappa * omega * mu)
    M = M.at[14, 22].set(omega * mu)
    M = M.at[14, 26].set(omega * mu)
    M = M.at[14, 30].set(omega * mu)
    M = M.at[14, 46].set(omega * mu)
    M = M.at[15, 16].set(kappa * mu)
    M = M.at[15, 19].set(kappa * omega * mu)
    M = M.at[15, 23].set(omega * mu)
    M = M.at[15, 27].set(omega * mu)
    M = M.at[15, 31].set(omega * mu)
    M = M.at[15, 47].set(omega * mu)
    M = M.at[16, 20].set(kappa * omega * mu)
    M = M.at[16, 24].set(omega * mu)
    M = M.at[16, 28].set(omega * mu)
    M = M.at[16, 32].set(omega * mu)
    M = M.at[16, 48].set(omega * mu)
    M = M.at[17, 18].set(kappa * mu)
    M = M.at[17, 19].set(mu)
    M = M.at[17, 20].set(mu)
    M = M.at[17, 21].set(omega * mu)
    M = M.at[17, 25].set(omega * mu)
    M = M.at[17, 33].set(omega * mu)
    M = M.at[17, 49].set(omega * mu)
    M = M.at[18, 19].set(mu)
    M = M.at[18, 20].set(mu)
    M = M.at[18, 22].set(omega * mu)
    M = M.at[18, 26].set(omega * mu)
    M = M.at[18, 34].set(omega * mu)
    M = M.at[18, 50].set(omega * mu)
    M = M.at[19, 20].set(kappa * mu)
    M = M.at[19, 23].set(omega * mu)
    M = M.at[19, 27].set(omega * mu)
    M = M.at[19, 35].set(omega * mu)
    M = M.at[19, 51].set(omega * mu)
    M = M.at[20, 24].set(omega * mu)
    M = M.at[20, 28].set(omega * mu)
    M = M.at[20, 36].set(omega * mu)
    M = M.at[20, 52].set(omega * mu)
    M = M.at[21, 22].set(kappa * mu)
    M = M.at[21, 23].set(omega * mu)
    M = M.at[21, 24].set(omega * mu)
    M = M.at[21, 25].set(kappa * omega * mu)
    M = M.at[21, 37].set(omega * mu)
    M = M.at[21, 53].set(omega * mu)
    M = M.at[22, 23].set(omega * mu)
    M = M.at[22, 24].set(omega * mu)
    M = M.at[22, 26].set(kappa * omega * mu)
    M = M.at[22, 38].set(omega * mu)
    M = M.at[22, 54].set(omega * mu)
    M = M.at[23, 24].set(kappa * mu)
    M = M.at[23, 27].set(kappa * omega * mu)
    M = M.at[23, 39].set(omega * mu)
    M = M.at[23, 55].set(omega * mu)
    M = M.at[24, 28].set(kappa * omega * mu)
    M = M.at[24, 40].set(omega * mu)
    M = M.at[24, 56].set(omega * mu)
    M = M.at[25, 26].set(kappa * mu)
    M = M.at[25, 27].set(mu)
    M = M.at[25, 28].set(mu)
    M = M.at[25, 41].set(omega * mu)
    M = M.at[25, 57].set(omega * mu)
    M = M.at[26, 27].set(mu)
    M = M.at[26, 28].set(mu)
    M = M.at[26, 42].set(omega * mu)
    M = M.at[26, 58].set(omega * mu)
    M = M.at[27, 28].set(kappa * mu)
    M = M.at[27, 43].set(mu)
    M = M.at[27, 59].set(omega * mu)
    M = M.at[28, 44].set(mu)
    M = M.at[28, 60].set(omega * mu)
    M = M.at[29, 30].set(kappa * mu)
    M = M.at[29, 31].set(mu)
    M = M.at[29, 32].set(omega * mu)
    M = M.at[29, 33].set(kappa * omega * mu)
    M = M.at[29, 37].set(omega * mu)
    M = M.at[29, 41].set(omega * mu)
    M = M.at[29, 45].set(kappa * omega * mu)
    M = M.at[30, 31].set(mu)
    M = M.at[30, 32].set(omega * mu)
    M = M.at[30, 34].set(kappa * omega * mu)
    M = M.at[30, 38].set(omega * mu)
    M = M.at[30, 42].set(omega * mu)
    M = M.at[30, 46].set(kappa * omega * mu)
    M = M.at[31, 32].set(kappa * omega * mu)
    M = M.at[31, 35].set(kappa * omega * mu)
    M = M.at[31, 39].set(omega * mu)
    M = M.at[31, 43].set(omega * mu)
    M = M.at[31, 47].set(kappa * omega * mu)
    M = M.at[32, 36].set(kappa * omega * mu)
    M = M.at[32, 40].set(omega * mu)
    M = M.at[32, 44].set(omega * mu)
    M = M.at[32, 48].set(kappa * omega * mu)
    M = M.at[33, 34].set(kappa * mu)
    M = M.at[33, 35].set(mu)
    M = M.at[33, 36].set(mu)
    M = M.at[33, 37].set(omega * mu)
    M = M.at[33, 41].set(omega * mu)
    M = M.at[33, 49].set(kappa * omega * mu)
    M = M.at[34, 35].set(mu)
    M = M.at[34, 36].set(mu)
    M = M.at[34, 38].set(omega * mu)
    M = M.at[34, 42].set(omega * mu)
    M = M.at[34, 50].set(kappa * omega * mu)
    M = M.at[35, 36].set(kappa * mu)
    M = M.at[35, 39].set(omega * mu)
    M = M.at[35, 43].set(omega * mu)
    M = M.at[35, 51].set(kappa * omega * mu)
    M = M.at[36, 40].set(omega * mu)
    M = M.at[36, 44].set(omega * mu)
    M = M.at[36, 52].set(kappa * omega * mu)
    M = M.at[37, 38].set(kappa * mu)
    M = M.at[37, 39].set(omega * mu)
    M = M.at[37, 40].set(omega * mu)
    M = M.at[37, 41].set(kappa * omega * mu)
    M = M.at[37, 53].set(kappa * omega * mu)
    M = M.at[38, 39].set(omega * mu)
    M = M.at[38, 40].set(omega * mu)
    M = M.at[38, 42].set(kappa * omega * mu)
    M = M.at[38, 54].set(kappa * omega * mu)
    M = M.at[39, 40].set(kappa * mu)
    M = M.at[39, 43].set(kappa * omega * mu)
    M = M.at[39, 55].set(kappa * omega * mu)
    M = M.at[40, 44].set(kappa * omega * mu)
    M = M.at[40, 56].set(kappa * omega * mu)
    M = M.at[41, 42].set(kappa * mu)
    M = M.at[41, 43].set(omega * mu)
    M = M.at[41, 44].set(omega * mu)
    M = M.at[41, 57].set(kappa * omega * mu)
    M = M.at[42, 43].set(omega * mu)
    M = M.at[42, 44].set(omega * mu)
    M = M.at[42, 58].set(kappa * omega * mu)
    M = M.at[43, 44].set(kappa * mu)
    M = M.at[43, 59].set(kappa * omega * mu)
    M = M.at[44, 60].set(kappa * omega * mu)
    M = M.at[45, 46].set(kappa * mu)
    M = M.at[45, 47].set(mu)
    M = M.at[45, 48].set(mu)
    M = M.at[45, 49].set(kappa * omega * mu)
    M = M.at[45, 53].set(omega * mu)
    M = M.at[45, 57].set(omega * mu)
    M = M.at[46, 47].set(mu)
    M = M.at[46, 48].set(mu)
    M = M.at[46, 50].set(kappa * omega * mu)
    M = M.at[46, 54].set(omega * mu)
    M = M.at[46, 58].set(omega * mu)
    M = M.at[47, 48].set(kappa * mu)
    M = M.at[47, 51].set(kappa * omega * mu)
    M = M.at[47, 55].set(omega * mu)
    M = M.at[47, 59].set(omega * mu)
    M = M.at[48, 52].set(kappa * omega * mu)
    M = M.at[48, 56].set(omega * mu)
    M = M.at[48, 60].set(omega * mu)
    M = M.at[49, 50].set(kappa * mu)
    M = M.at[49, 51].set(mu)
    M = M.at[49, 52].set(mu)
    M = M.at[49, 53].set(omega * mu)
    M = M.at[49, 57].set(omega * mu)
    M = M.at[50, 51].set(mu)
    M = M.at[50, 52].set(mu)
    M = M.at[50, 54].set(omega * mu)
    M = M.at[50, 58].set(omega * mu)
    M = M.at[51, 52].set(kappa * mu)
    M = M.at[51, 55].set(omega * mu)
    M = M.at[51, 59].set(omega * mu)
    M = M.at[52, 56].set(omega * mu)
    M = M.at[52, 60].set(omega * mu)
    M = M.at[53, 54].set(kappa * mu)
    M = M.at[53, 55].set(omega * mu)
    M = M.at[53, 56].set(omega * mu)
    M = M.at[53, 57].set(kappa * omega * mu)
    M = M.at[54, 55].set(omega * mu)
    M = M.at[54, 56].set(omega * mu)
    M = M.at[54, 58].set(kappa * omega * mu)
    M = M.at[55, 56].set(kappa * mu)
    M = M.at[55, 59].set(kappa * omega * mu)
    M = M.at[56, 60].set(kappa * omega * mu)
    M = M.at[57, 58].set(kappa * mu)
    M = M.at[57, 59].set(mu)
    M = M.at[57, 60].set(mu)
    M = M.at[58, 59].set(mu)
    M = M.at[58, 60].set(mu)
    M = M.at[59, 60].set(kappa * mu)

    for i in range(61):
        # M[i] *= pi
        M = M.at[i].multiply(pi)

    # Fill in the lower triangle
    for i in range(61):
        for j in range(61):
            # M[j, i] = M[i, j]
            M = M.at[j, i].set(M[i, j])

    # Compute the diagonal
    rowsum = jnp.sum(M, axis=1)
    for i in range(61):
        # M[i, i] = -rowsum[i]
        M = M.at[i, i].set(-rowsum[i])
        
    return M


def likelihood(X, n, mu, kappa, omega, pi):
    alpha_Ai = jnp.zeros((61, 61))
    alpha_A = [0] * 61
    lik_full = [0] * 61

    mutmat = substitution_rate_matrix(mu, kappa, omega, pi)
    a, b = jnp.linalg.eig(mutmat) # eigenvalues and eigenvectors
    V = b
    D = [1 / (1-A) for A in a.tolist()]
    VD = V
    for i in range(61):
        # VD[i] *= D
        VD = VD.at[i].multiply(D)

    for A in range(61):
        lik = 0
        m_AA = jnp.dot(V[A, :], VD[A, :])

        if m_AA < 1e-6:
            m_AA = 1e-6

        for i in range(61):
            if A == i:
                # alpha_Ai[A, i] = 1
                alpha_Ai = alpha_Ai.at[A, i].set(1)
            else:
                m_Ai = jnp.dot(V[A, :], VD[i, :])

                if m_Ai < 1e-6:
                    m_Ai = 1e-6
                # alpha_Ai[A, i] = m_Ai / m_AA
                alpha_Ai = alpha_Ai.at[A, i].set(m_Ai / m_AA)
            lik += lgamma(X[i] + alpha_Ai[A, i]) - lgamma(alpha_Ai[A, i]) -lgamma(X[i] + 1)

        # alpha_A[A] = jnp.sum(alpha_Ai[A, :])
        alpha_A = alpha_A.at[A].set(jnp.sum(alpha_Ai[A, :]))
        # lik_full[A] = lik + lgamma(alpha_A[A]) - lgamma(n + alpha_A[A]) + lgamma(n + 1)
        lik_full = lik_full.at[A].set(lik + lgamma(alpha_A[A]) - lgamma(n + alpha_A[A]) + lgamma(n + 1))
    partial_ll += logsumexp(lik_full + jnp.log(pi))

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


