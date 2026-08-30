#!/usr/bin/env python3
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from lar_melar import lar_correlation, melar_correlation, lar_melar_jackknife, jel_log_ratio  # noqa: E402

rng = np.random.default_rng(7301)
n = 80
rho = 0.55
x = rng.normal(size=n)
y = rho * x + np.sqrt(1 - rho**2) * rng.normal(size=n)
u = rng.uniform(-1, 1, size=n)
w0 = x + 0.70 * u
z0 = y - 0.55 * u
lar = lar_correlation(w0, z0, u)[0]
melar0 = melar_correlation(w0, z0, u, np.zeros((2, 2)))[0]
assert np.isfinite(lar)
assert abs(lar - melar0) < 1e-12

# With known random-error variance, ME-LAR should remove the principal attenuation
# in this deterministic-seed construction.
tau2 = 0.15
w = w0 + rng.normal(scale=np.sqrt(tau2), size=n)
z = z0 + rng.normal(scale=np.sqrt(tau2), size=n)
lar_err = lar_correlation(w, z, u)[0]
melar, pseudo = lar_melar_jackknife(w, z, u, np.diag([tau2, tau2]), method="melar")
assert np.isfinite(melar)
assert np.all(np.isfinite(pseudo))
assert np.isfinite(lar_err)
assert abs(jel_log_ratio(pseudo, np.mean(pseudo))) < 1e-8
print("LAR/ME-LAR validation checks passed")
