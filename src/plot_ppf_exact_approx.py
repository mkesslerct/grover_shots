import pathlib

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

from collections import defaultdict

from compute_exp_var import compute_expectation, compute_variance

FIGURES_DIRECTORY = pathlib.Path("..") / "figures"
RESULTS_DIRECTORY = pathlib.Path("..") / "results"

from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple


def compute_n_approx(p, A, M, pg, case):
    if case == "all":
        n_approx = (
            compute_expectation(pg, A, M)
            + M / pg * (np.log(2) - np.euler_gamma)
            - M / pg * np.log(chi2.ppf(1 - p, df=2))
        )
    elif case == "proportion":
        n_approx = compute_expectation(pg, A, M) + +np.sqrt(
            compute_variance(pg, A, M)
        ) * norm.ppf(p)
    else:
        raise Exception(f"Invalid value for case parameter: {case}")

    return n_approx


def plot_ppf(A, M, n, pg_values, case="exact"):
    cum_prob = dict()
    prob_file_exists = defaultdict(lambda: False)
    num_prob_files = 0
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
            prob_file_exists[pg] = True
            cum_prob[pg] = pd.read_csv(prob_file_path)
            cum_prob[pg] = cum_prob[pg].loc[cum_prob[pg]["F"] < 0.99]
            print(f"pg: {pg}, final proba: {cum_prob[pg]['F'].iloc[-1]}")
            (lines[pg],) = ax.plot(
                cum_prob[pg]["F"], cum_prob[pg]["t"], label=f"$p_G$ = {pg}"
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

            num_prob_files += 1
        else:
            print(f"{prob_file_path} does not exists")

    # if num_prob_files > 0:
    #     for p_threshold in [0.7, 0.85, 0.95]:
    #         ax.axvline(x=p_threshold, ymin=0, ls="--")

    #         ax.annotate(
    #             f"{p_threshold}",
    #             xy=(p_threshold, A + 1),
    #             xytext=(6, 2),  # 4 points vertical offset.
    #             textcoords="offset points",
    #             ha="left",
    #             va="bottom",
    #             fontsize="small",
    #         )

    #         for i, pg in enumerate([p for p in pg_values if prob_file_exists[p]]):
    #             if np.max(cum_prob[pg]["F"]) >= p_threshold:
    #                 idx_threshold = (
    #                     cum_prob[pg]["F"].values >= p_threshold
    #                 ).searchsorted(True, side="left")
    #                 M_threshold = cum_prob[pg]["t"].iloc[idx_threshold]
    #                 xytext = (-16, 4) if i < 2 else (14, -15)
    #                 ax.annotate(
    #                     f"{M_threshold}",
    #                     xy=(p_threshold, M_threshold),
    #                     xytext=xytext,
    #                     textcoords="offset points",
    #                     ha="center",
    #                     va="bottom",
    #                     color=f"C{i}",
    #                 )
    # else:
    #     print(f"No csv file")

    l = ax.legend(
        [(l1, l2) for (l1, l2) in zip(lines.values(), lines_approx.values())],
        [f"$p_G$: {pgv}" for pgv in pg_values],
        numpoints=1,
        handler_map={tuple: HandlerTuple(ndivide=None)},
    )

    ax.add_artist(l)
    # legend_approx = ax.legend(handles=[l for l in lines_approx.values()], loc="lower center")
    # ax.add_artist(legend_approx)
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
