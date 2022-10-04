import pathlib
import sys

from gmpy2 import mpz, mpfr, bincoef, mul, sub, add

RESULTS_DIRECTORY = pathlib.Path("..") / "results"

def bincoef_array(K):
    """returns a list bincoef[i] = bincoef(K, i + 1),
    0 <= i <= K - 1
    [bincoef(K, 1), bincoef(K, 2), ....., bincoef(K, K)]
    """
    b = [bincoef(K, i + 1) for i in range(K)]

    return b


def factK_Stirling(K, l, bincoef_a):
    """Implements K! Stirling(K, l), K <= l of the second kind
    Subsection.~28.6(i)]{NIST:DLMF}Subsection.~28.6(i)]{NIST:DLMF}http://dlmf.nist.gov/
    bincoef_a is a list, it contains
      [bincoef(K, 1), bincoef(K, 2), ....., bincoef(K, K)]
    """
    f_S = mpz("0")
    for r in range(1, K + 1):
        prod = mul(pow(r, l), bincoef_a[r - 1])
        f_S = add(f_S, pow(-1, (K - r)) * prod)

    return f_S


def factA_Stirling_array(A, m):
    """Returns A! Stirling(A, l), for l in range(A, m + 1) as a list
    A <= m"""
    bincoef_a = bincoef_array(A)
    return [factK_Stirling(A, j, bincoef_a) for j in range(A, m + 1)]


def F_array(A, M, s, pg, F_threshold):
    """Computes [P(X_{A, M} <= t) for t in range(m, s + 1)] where
    X_{A, M} is the drawing on which, for the first time, the number of distinct
    solutions hat have been sampled is A + 1.
    If the cumulative distribution function P(X_{A, M} <= t) exceeds F_threshold
    the computation stops before reaching s.
    The values of the pmf P(X_{A, M} = t) and cdf P(X_{A, M} <= t) are written to
    a csv file in the RESULTS_DIRECTORY.
    It implements the formula (10) of Proposition 2.3 of the paper.
    """
    factA_S_array = factA_Stirling_array(A, s - 1)

    pgg = mpfr(pg)
    prg = mpfr(pg / M)

    sum_proba = mpfr("0")
    file_path = (
        RESULTS_DIRECTORY / f"pmf_cdf_A_{A}_M_{M}_s_{s}_{int(pg * 1000)}.csv"
    )
    results_file = open(file_path, "w")
    results_file.write("t,p,F\n")

    for t in range(A + 1, s + 1):  # t is "s" in the exact computation formula
        if sum_proba > F_threshold:
            break
        proba = mpfr("0")
        for l in range(A, t):
            proba = add(
                proba,
                mul(
                    mul(
                        mul(A + 1, mul(pow(1 - pgg, t - l - 1), pow(prg, l + 1))),
                        bincoef(t - 1, l),
                    ),
                    factA_S_array[l - A],
                ),
            )
            proba = mul(proba, bincoef(M, A + 1))
        sum_proba = add(sum_proba, proba)
        results_file.write(
            f"{t},{'{:.6Nf}'.format(proba)},{'{:.6Nf}'.format(sum_proba)}\n"
        )
        print(f"pg: {pg}, F({t}) computed, value: {'{:.6Nf}'.format(sum_proba)}")

    results_file.close()


if __name__ == "__main__":
    args = sys.argv
    if not (len(args) == 5):
        print(f"Usage: {args[0]} <number A> <number M> <float pg> <number s>")

    A = int(args[1])
    M = int(args[2])
    pg = float(args[3])
    s = int(args[4])

    assert ((A > 0) & (M > A)), "A must be > 0 and < M" 
    assert((pg > 0) & (pg < 1)), "pg must be > 0 and < 1"
    assert(s >= M), "s must be >= M"
                 
    F_array(A, M, s, pg, F_threshold=0.99)
    # Note:  if the cumulative distribution function P(X_{A, M} <= t) exceeds
    # F_threshold the computation stops before reaching s.









