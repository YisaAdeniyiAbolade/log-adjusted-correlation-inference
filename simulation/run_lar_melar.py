#!/usr/bin/env python3
"""Reproduce the reference-free LAR/ME-LAR Monte Carlo summary."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from lar_melar import lar_correlation, lar_melar_jackknife, jel_ci  # noqa: E402


def generate_sample(n, rho, tau2, u_distribution, distortion, rng):
    x = rng.normal(size=n)
    y = rho * x + math.sqrt(1.0 - rho * rho) * rng.normal(size=n)

    if u_distribution == "Uniform":
        u = rng.uniform(-1.0, 1.0, n)
        eu2 = 1.0 / 3.0
    elif u_distribution == "Normal":
        u = rng.normal(size=n)
        eu2 = 1.0
    elif u_distribution == "Beta":
        u = 2.0 * rng.beta(2.0, 2.0, n) - 1.0
        eu2 = 1.0 / 5.0
    else:
        raise ValueError(f"unknown U distribution: {u_distribution}")

    if distortion == "Linear":
        q = u
        gx, gy = 0.70 * q, -0.55 * q
    elif distortion == "Quadratic":
        q = u * u - eu2
        if u_distribution == "Normal":
            q = q / math.sqrt(2.0)
        gx, gy = 0.55 * q, -0.45 * q
    elif distortion == "Periodic":
        q = np.sin(np.pi * u)
        gx, gy = 0.75 * q, -0.60 * q
    else:
        raise ValueError(f"unknown distortion: {distortion}")

    ex = rng.normal(scale=math.sqrt(tau2), size=n)
    ey = rng.normal(scale=math.sqrt(tau2), size=n)
    w = x + gx + ex
    z = y + gy + ey
    return w, z, u


def summarize_scenario(row, reps):
    sid = row.scenario_id
    n = int(row.n)
    rho = float(row.rho_true)
    tau2 = float(row.tau2)
    E = np.diag([tau2, tau2])
    records = []
    scenario_index = int(sid[1:]) - 1

    for b in range(reps):
        rng = np.random.default_rng(20260830 + scenario_index * 100000 + b)
        w, z, u = generate_sample(n, rho, tau2, row.u_distribution, row.distortion, rng)
        naive = float(np.corrcoef(w, z)[0, 1])
        lar = lar_correlation(w, z, u)[0]
        melar, pseudo = lar_melar_jackknife(w, z, u, E, method="melar")
        jci = jel_ci(pseudo, adjusted=False)
        aci = jel_ci(pseudo, adjusted=True)
        finite_j = bool(np.all(np.isfinite(jci)))
        finite_a = bool(np.all(np.isfinite(aci)))
        records.append({
            "naive": naive,
            "lar": lar,
            "melar": melar,
            "jel_cover": float(finite_j and jci[0] <= rho <= jci[1]),
            "ajel_cover": float(finite_a and aci[0] <= rho <= aci[1]),
            "jel_length": float(jci[1] - jci[0]) if finite_j else np.nan,
            "ajel_length": float(aci[1] - aci[0]) if finite_a else np.nan,
            "finite_jel": float(finite_j),
            "finite_ajel": float(finite_a),
        })

    d = pd.DataFrame(records)
    return {
        "scenario_id": sid,
        "u_distribution": row.u_distribution,
        "distortion": row.distortion,
        "n": n,
        "rho_true": rho,
        "tau2": tau2,
        "naive_bias": float((d.naive - rho).mean()),
        "naive_rmse": float(np.sqrt(((d.naive - rho) ** 2).mean())),
        "lar_bias": float((d.lar - rho).mean()),
        "lar_rmse": float(np.sqrt(((d.lar - rho) ** 2).mean())),
        "melar_bias": float((d.melar - rho).mean()),
        "melar_rmse": float(np.sqrt(((d.melar - rho) ** 2).mean())),
        "melar_jel_coverage": float(d.jel_cover.mean()),
        "melar_ajel_coverage": float(d.ajel_cover.mean()),
        "melar_jel_length": float(d.jel_length.mean()),
        "melar_ajel_length": float(d.ajel_length.mean()),
        "finite_point": float(np.isfinite(d.melar).mean()),
        "finite_jel": float(d.finite_jel.mean()),
        "finite_ajel": float(d.finite_ajel.mean()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=1000)
    parser.add_argument("--scenarios", default=str(ROOT / "simulation" / "lar_melar_scenarios.csv"))
    parser.add_argument("--output", default=str(ROOT / "results" / "simulation" / "lar_melar_summary.csv"))
    args = parser.parse_args()
    scenarios = pd.read_csv(args.scenarios)
    rows = []
    for row in scenarios.itertuples(index=False):
        result = summarize_scenario(row, args.reps)
        rows.append(result)
        print(row.scenario_id, f"ME-LAR bias={result['melar_bias']:.4f}", f"AJEL={result['melar_ajel_coverage']:.3f}")
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
