import pathlib

import pickle

import numpy as np
import pandas as pd

from collections import defaultdict

from compute_exp_var import compute_expectation, compute_variance

FIGURES_DIRECTORY = pathlib.Path("..") / "figures"
RESULTS_DIRECTORY = pathlib.Path("..") / "results"

from matplotlib import pyplot as plt

fig, ax = plt.subplots()

cum_prob = dict()

pg_values = [0.7, 0.95, 0.999]

M = 10
A = 9
n = 150
pg_values = [0.7, 0.95, 0.999]

prob_file_exists = defaultdict(lambda: False)
num_prob_files = 0
prob_filename = {
    pg_values[k]: f"pmf_cdf_A_{A}_M_{M}_n_{n}_{int(pg_values[k] * 1000)}.csv"
    for k in range(3)
}

for i, pg in enumerate(pg_values):
    pr = pg / M

    prob_file_path = RESULTS_DIRECTORY / prob_filename[pg]

    if prob_file_path.is_file():
        prob_file_exists[pg] = True
        cum_prob[pg] = pd.read_csv(prob_file_path)
        cum_prob[pg] = cum_prob[pg].loc[cum_prob[pg]["F"] < 0.99]
        print(f"pg: {pg}, final proba: {cum_prob[pg]['F'].iloc[-1]}")

        ax.plot(cum_prob[pg]["F"], cum_prob[pg]["t"], label=f"pg = {pg}")
        ax.set_xlim(-0.05, 1.05)
        ax.set_title(f"M: {M}, A: {A}")

        num_prob_files += 1
    else:
        print(f"{prob_file_path} does not exists")

if num_prob_files > 0:
    for p_threshold in [0.7, 0.85, 0.95]:
        ax.axvline(x=p_threshold, ymin=0, ls="--")

        ax.annotate(
            f"{p_threshold}",
            xy=(p_threshold, A + 1),
            xytext=(6, 2),  # 4 points vertical offset.
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize="small",
        )

        for i, pg in enumerate([p for p in pg_values if prob_file_exists[p]]):
            if np.max(cum_prob[pg]["F"]) >= p_threshold:
                idx_threshold = (cum_prob[pg]["F"].values >= p_threshold).searchsorted(
                    True, side="left"
                )
                M_threshold = cum_prob[pg]["t"].iloc[idx_threshold]
                xytext = (-16, 4) if i < 2 else (14, -15)
                ax.annotate(
                    f"{M_threshold}",
                    xy=(p_threshold, M_threshold),
                    xytext=xytext,
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    color=f"C{i}",
                )

    ax.legend()
    ax.set_xlabel("$P(X_{A,M} \leq n)$")
    ax.set_ylabel("n")
    fig.tight_layout()

    plt.savefig(FIGURES_DIRECTORY / f"n_vs_F_A_{A}_M_{M}_csv.pdf", transparent=True)
    plt.show()

else:
    print(f"No csv file")
