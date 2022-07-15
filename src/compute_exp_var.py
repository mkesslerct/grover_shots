import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import chi2
import pathlib


FIGURES_DIRECTORY = pathlib.Path("..") / "figures"
RESULTS_DIRECTORY = pathlib.Path("..") / "results"

def compute_expectation(pg, A, M):
    expec = M / pg * np.sum(1 / np.arange(M - A, M + 1))
    
    return expec 

# def approximate_expectation(pg, M):
#      return (M * np.log(M) + M * np.euler_gamma + 0.5)  / pg   

def compute_variance(pg, A, M):
    probs = pg * np.arange(M - A, M + 1) / M 
    var = np.sum( (1 - probs) / (probs ** 2))

    return var

# def approximate_variance(pg, M):
#     expansion = (np.pi ** 2) / ( 6 * pg ** 2) * (M ** 2) #- M * np.log(M ) / pg - np.euler_gamma  * M / pg - 1 / (2 * pg)

#     return expansion

