import pathlib

import numpy as np
import pandas as pd

from utilities import compute_expectation, compute_variance, compute_n_approx

FIGURES_DIRECTORY = pathlib.Path("..") / "figures"
RESULTS_DIRECTORY = pathlib.Path("..") / "results"

from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerTuple

def plot_ppf(A, M, n, pg_values, case="exact"):
    prob_filename = {
        pg_values[k]: f"pmf_cdf_A_{A}_M_{M}_n_{n}_{int(pg_values[k] * 1000)}.csv"
        for k in range(3)
    }

    fig, ax = plt.subplots()
    lines = dict()
    lines_approx = dict()

    for i, pg in enumerate(pg_values):
        pr = pg / M

        prob_file_path = RESULTS_DIRECTORY / prob_filename[pg]

        if prob_file_path.is_file():
            cum_prob = pd.read_csv(prob_file_path)
            cum_prob = cum_prob.loc[cum_prob["F"] < 0.99]
            print(f"pg: {pg}, final proba: {cum_prob['F'].iloc[-1]}")
            (lines[pg],) = ax.plot(
                cum_prob["F"], cum_prob["t"], label=f"$p_G$ = {pg}"
            )

            if case in ["all", "proportion"]:
                p = np.linspace(0.05, 0.95, 19)
                n_approx = compute_n_approx(p, A, M, pg, case)
                (lines_approx[pg],) = ax.plot(
                    p,
                    n_approx,
                    linestyle=":",
                    marker="+",
                    color=f"C{i}",
                    label=f"$p_G$ = {pg}, approximation of ppf",
                )

            ax.set_xlim(-0.05, 1.05)
            ax.set_title(f"M: {M}, A: {A}")

        else:
            print(f"{prob_file_path} does not exists")

    l = ax.legend(
        [(l1, l2) for (l1, l2) in zip(lines.values(), lines_approx.values())],
        [f"$p_G$: {pgv}" for pgv in pg_values],
        numpoints=1,
        handler_map={tuple: HandlerTuple(ndivide=None)},
    )

    ax.add_artist(l)
    ax.set_xlabel("$P(X_{A,M} \leq n)$")
    ax.set_ylabel("n")
    fig.tight_layout()

    plt.savefig(FIGURES_DIRECTORY / f"n_vs_F_A_{A}_M_{M}_{case}.pdf", transparent=True)
    plt.show()


if __name__ == "__main__":
    M_values = [10, 10, 100, 100, 1000, 1000]
    A_values = [4, 9, 49, 99, 499, 999]
    n_values = [50, 150, 500, 2000, 5000, 15000]
    pg_values = [0.7, 0.95, 0.999]
    for (A, M, n) in zip(A_values, M_values, n_values):
        plot_ppf(A, M, n, pg_values, case=("all" if (A + 1 == M) else "proportion"))
