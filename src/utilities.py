import numpy as np
from scipy.stats import chi2,  norm

def compute_expectation(pg, A, M):
    expec = M / pg * np.sum(1 / np.arange(M - A, M + 1))
    
    return expec 

def compute_variance(pg, A, M):
    probs = pg * np.arange(M - A, M + 1) / M 
    var = np.sum( (1 - probs) / (probs ** 2))

    return var

def compute_n_approx(p, A, M, pg, case):
    if case == "all":
        n_approx = (
            compute_expectation(pg, A, M)
            + M / pg * (np.log(2) - np.euler_gamma)
            - M / pg * np.log(chi2.ppf(1 - p, df=2))
        )
    elif case == "proportion":
        n_approx = compute_expectation(pg, A, M) + np.sqrt(
            compute_variance(pg, A, M)
        ) * norm.ppf(p)
    else:
        raise Exception(f"Invalid value for case parameter: {case}")

    return n_approx
