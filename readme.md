# Introduction
This repository contains the code  that implements the formulae to be found in the manuscript " " (2022) by Mathieu Kessler[^1], Diego Alonso[^2] and Pedro Sánchez[^2].  They allow to compute the required number of shots to ensure, with probability $p$, to observe all, or a fraction of all the solutions to a search problem, using the Grover's algorithm.

# Contents:
Following the index of the manuscript, three levels of results are obtained:

##   Expectation and variance of the required number of shots 


  The Python script utilities.py contains two functions compute_expectation and compute_variance which implement equations (4) and (5) in the manuscript.
  
  They only require numpy and can be used as: 
  
``` python
  print("The average number of shots to observe all 100 solutions to the search problem")
  print("if the probability of success is pg = 0.95 is:")
  print(compute_expectation(0.95, 99, 100)
  print(f"while its variance is {compute_variance(0.95, 99, 100)}")
```

## Approximations to the number of required shots
The function compute_n_approx implements equations (7) and (9) of the manuscript, which provide the number s of shots to ensure,  with  given probability p, to observe all (equation (7)) or a fraction of all (equation (9)) solutions. The function implements an argument "case", which admits the value "all" for the implementation of (7) and "proportion" for the implementation of (9). 

It requires numpy and scipy and can be used as:

``` python
print("The number of shots required to ensure that, with probability 0.9, we observe")
print("all 100 solutions to a search problem, with pg = 0.95 is:")
print("compute_n_approx(0.9, 99, 100, 0.9, "all")
```

``` python
print("The number of shots required to ensure that, with probability 0.9, we observe")
print("at least half of the 100 solutions to a search problem, with pg = 0.95 is:")
print("compute_n_approx(0.9, 49, 100, 0.9, "proportion")
```

##  Exact probability mass function (pmf)
The implementation of formula (10) requires, for large values of $M$ a high precision and range library. This repository contains an implementation in C and an implementation in Python of the exact pmf expression (10). They both use the gmp library. 

### Implementation in C:
To compile in Linux, assuming that you are in the src directory which contain the .c file:

``` shell
gcc -o gmp_FX_A_M gmp_FX_A_M.c -O2 -lgmp 
```

You can then specify A, M, pg and s as command line arguments:
``` shell
./gmp_FX_A_M 99 100 0.7 1500
```

Notes:
1. check that a "results" folder exists at the same level as the containing src folder. The results will be written as a csv file into that folder, each row correspond to an integer $t$
1. You can use the approximate formula, equations (7) or (9), to have an estimate of the value of $s$ you need to reach the probability you are interested in.

### Implementation in Python

You should install gmpy2, for example from conda-forge if you use a conda distribution.

The script gmpy2_FX_A_M.py contains the function F_array which admits as arguments: A, M, pg,  an integer s which provides the upper limit for the computation of the cumulative distribution function. It also admits the parameter F_threshold which allows to stop the computation if a given value is reached for the cumulative distribution function. 
The values of the pmf P(X_{A, M} = t) and cdf P(X_{A, M} <= t) are written to
a csv file in the RESULTS_DIRECTORY.
Example of use: 
  
``` python
M = 100
A = 100 - 1  # number of distinct solutions that have been sampled: A + 1
s = 1350
pg = 0.9
F_array(A, M, s, pg, F_threshold=0.98)
```

### Example of computation of quantile: 
Assume that you have used either the Python script or the C implementation to compute the values of the cumulative distribution F, for A, M, a value pg and s. The results are stored in an csv file called, for example, `pmf_cdf_A_9_M_10_s_200_700.csv` in the `results` folder.  In order to compute the value of $s$ that fulfills $P(X_{A, M} <= s) = p$, the Python function `compute_quantile` from the `utilities.py` file can be used:

``` python
from pathlib import Path
import pandas as pd 
from utilities import compute_quantile

RESULTS_DIR = Path("..") / "results"

# load the cdf dataframe from the csv file
F = pd.read_csv(RESULTS_DIR / "pmf_cdf_A_9_M_10_s_200_700.csv")
# compute the quantile for p = 0.95 for example
p = 0.95
print(f"The required number of shots to observe, with probability {p},")
print(f"all solutions is {compute_quantile(p, F)}.")
```



[^1]: Department of Applied Mathematics and Statistics. Universidad Politécnica de Cartagena, Campus Muralla del Mar, 30202, Cartagena (Murcia), Spain
[^2]: Department of Information Technology and Communications. Universidad Politécnica de Cartagena, Campus Muralla del Mar, 30202, Cartagena (Murcia), Spain

