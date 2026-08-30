#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
df = pd.read_csv(ROOT / "results" / "simulation" / "lar_melar_summary.csv")

fig, ax = plt.subplots(figsize=(8.2, 5.2))
for name, g in df.groupby("distortion", sort=False):
    g = g.sort_values("n")
    ax.plot(g["n"], g["lar_bias"], marker="o", label=f"LAR - {name}")
    ax.plot(g["n"], g["melar_bias"], marker="s", linestyle="--", label=f"ME-LAR - {name}")
ax.axhline(0, linewidth=1)
ax.set_xlabel("Sample size n")
ax.set_ylabel("Bias")
ax.set_title("Bias of LAR and ME-LAR under log-scale distortion")
ax.legend(ncol=2, fontsize=8)
fig.tight_layout()
fig.savefig(ROOT / "figures" / "simulation_lar_melar_bias.png", dpi=220)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.2, 5.2))
for name, g in df.groupby("distortion", sort=False):
    g = g.sort_values("n")
    ax.plot(g["n"], g["melar_jel_coverage"], marker="o", label=f"JEL - {name}")
    ax.plot(g["n"], g["melar_ajel_coverage"], marker="s", linestyle="--", label=f"AJEL - {name}")
ax.axhline(0.95, linewidth=1)
ax.set_xlabel("Sample size n")
ax.set_ylabel("Empirical coverage probability")
ax.set_ylim(0.90, 1.00)
ax.set_title("ME-LAR 95% JEL and AJEL coverage")
ax.legend(ncol=2, fontsize=8)
fig.tight_layout()
fig.savefig(ROOT / "figures" / "simulation_lar_melar_coverage.png", dpi=220)
plt.close(fig)
