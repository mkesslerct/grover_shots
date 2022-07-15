import gmpy2
from gmpy2 import mpz, mpfr, bincoef, mul, sub, add

import pathlib


FIGURES_DIRECTORY = pathlib.Path("..") / "figures"

RESULTS_DIRECTORY = pathlib.Path("..") / "results"

# gmpy2.get_context().precision = 200000


def bincoef_array(K):
    """returns a list bincoef[i] = bincoef(K, i + 1),
    0 <= i <= K - 1
    [bincoef(K, 1), bincoef(K, 2), ....., bincoef(K, K)]
    """
    b = [bincoef(K, i + 1) for i in range(K)]

    return b


def factK_Stirling(K, l, bincoef_a):
    """Implements K! Stirling(K, l), K <= l of the second kind
    from p 19, H. S. Wilf, generatingfunctionology, 2nd ed., Academic Press, Boston, 1994
    and p 1037, Handbook of Integrals and sums, https://multimediarepository.blob.core.windows.net/imagecontainer/92688b5c8bf4478d80231406c207b221.png
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


def F_array(A, M, n, pg_values, F_threshold):
    """Returns [P(W_M <= t) for t in range(m, n + 1)] where
    W_M is the drawing on which, for the first time, the number of distinct
    red balls that have been sampled is A + 1
    If m is not is A + 1, sum_prob_prev is requested
    Implemented formula:
    https://multimediarepository.blob.core.windows.net/imagecontainer/98085adef8164a7997764dafef9b24dc.jpeg
    """
    factA_S_array = factA_Stirling_array(A, n - 1)

    for pg in pg_values:
        pgg = mpfr(pg)
        prg = mpfr(pg / M)

        sum_proba = mpfr("0")
        file_path = (
            RESULTS_DIRECTORY / f"pmf_cdf_A_{A}_M_{M}_n_{n}_{int(pg * 1000)}_gmpy2.csv"
        )
        results_file = open(file_path, "w")
        results_file.write("t,p,F\n")

        for t in range(A + 1, n + 1):  # t is "n" in the exact computation formula
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
    M = 100
    A = 100 - 1  # number of distinct red balls that have sampled: A + 1
    n = 1350
    pg_values = [0.7, 0.95, 0.999]

    F_array(A, M, n, pg_values, F_threshold=0.98)
