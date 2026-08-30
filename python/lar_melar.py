"""Reference-free LAR/ME-LAR estimators and one-sample JEL/AJEL utilities."""

from __future__ import annotations

import math
import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2


def epanechnikov_kernel(t: np.ndarray) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    out = 0.75 * (1.0 - t * t)
    out[np.abs(t) > 1.0] = 0.0
    return out


def corrected_correlation(x, y, error_cov=None, clip=True):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 3:
        raise ValueError("x and y must be one-dimensional arrays of equal length >= 3")
    S = np.cov(np.column_stack([x, y]), rowvar=False, ddof=1)
    if error_cov is None:
        error_cov = np.zeros((2, 2))
    E = np.asarray(error_cov, dtype=float)
    if E.shape != (2, 2):
        raise ValueError("error_cov must be 2 x 2")
    L = (S - E + (S - E).T) / 2.0
    if not np.all(np.isfinite(L)) or L[0, 0] <= 1e-12 or L[1, 1] <= 1e-12:
        return np.nan
    rho = L[0, 1] / math.sqrt(L[0, 0] * L[1, 1])
    if clip:
        rho = float(np.clip(rho, -0.999999, 0.999999))
    return float(rho)


def _log_mean_exp(x):
    x = np.asarray(x, dtype=float)
    a = float(np.max(x))
    return a + math.log(float(np.mean(np.exp(x - a))))


def lar_adjust(values, u, bandwidth=None, return_kernel=False):
    values = np.asarray(values, dtype=float)
    u = np.asarray(u, dtype=float)
    n = len(values)
    if values.ndim != 1 or u.ndim != 1 or len(u) != n or n < 3:
        raise ValueError("values and u must have equal one-dimensional length >= 3")
    if bandwidth is None:
        bandwidth = max(float(np.std(u, ddof=1)) * n ** (-1.0 / 3.0), 1e-8)
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError("bandwidth must be positive")

    D = (u[:, None] - u[None, :]) / bandwidth
    K = epanechnikov_kernel(D)
    den = K.sum(axis=1)
    if np.any(den <= 0):
        raise FloatingPointError("kernel denominator is zero")

    vmax = float(np.max(values))
    e = np.exp(values - vmax)
    conditional_log_mean = vmax + np.log((K @ e) / den)
    global_log_mean = _log_mean_exp(values)
    adjusted = values - conditional_log_mean + global_log_mean
    if return_kernel:
        return adjusted, float(bandwidth), K
    return adjusted


def lar_correlation(w, z, u, bandwidth=None):
    x_adj = lar_adjust(w, u, bandwidth)
    y_adj = lar_adjust(z, u, bandwidth)
    return corrected_correlation(x_adj, y_adj, np.zeros((2, 2))), x_adj, y_adj


def melar_correlation(w, z, u, error_cov, bandwidth=None):
    x_adj = lar_adjust(w, u, bandwidth)
    y_adj = lar_adjust(z, u, bandwidth)
    return corrected_correlation(x_adj, y_adj, error_cov), x_adj, y_adj


def lar_melar_jackknife(w, z, u, error_cov=None, method="melar", bandwidth=None):
    """Delete-one pseudo-values with the full-sample bandwidth held fixed.

    The kernel smoother and adjusted second moments are recomputed for every deletion.
    The implementation evaluates all delete-one smoothers by removing one row/column
    contribution from the full kernel matrix, which is algebraically equivalent to
    refitting with the common bandwidth and is substantially faster than an explicit loop.
    """
    w = np.asarray(w, dtype=float)
    z = np.asarray(z, dtype=float)
    u = np.asarray(u, dtype=float)
    n = len(w)
    if len(z) != n or len(u) != n or n < 4:
        raise ValueError("w, z and u must have equal length >= 4")
    if method not in {"lar", "melar"}:
        raise ValueError("method must be 'lar' or 'melar'")
    if error_cov is None:
        error_cov = np.zeros((2, 2))
    E = np.asarray(error_cov, dtype=float)
    if E.shape != (2, 2):
        raise ValueError("error_cov must be 2 x 2")
    if method == "lar":
        E = np.zeros((2, 2))
    if bandwidth is None:
        bandwidth = max(float(np.std(u, ddof=1)) * n ** (-1.0 / 3.0), 1e-8)

    D = (u[:, None] - u[None, :]) / bandwidth
    K = epanechnikov_kernel(D)
    den = K.sum(axis=1)

    def full_adjust(values):
        vmax = float(np.max(values))
        e = np.exp(values - vmax)
        num = K @ e
        conditional = vmax + np.log(num / den)
        global_mean = vmax + math.log(float(e.mean()))
        return values - conditional + global_mean

    x_full = full_adjust(w)
    y_full = full_adjust(z)
    theta = corrected_correlation(x_full, y_full, E)

    def leave_adjust_matrix(values):
        vmax = float(np.max(values))
        e = np.exp(values - vmax)
        num = K @ e
        den_minus = den[:, None] - K
        num_minus = num[:, None] - K * e[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            conditional = vmax + np.log(num_minus / den_minus)
        global_mean = vmax + np.log((e.sum() - e) / (n - 1))
        adjusted = values[:, None] - conditional + global_mean[None, :]
        np.fill_diagonal(adjusted, 0.0)
        return adjusted

    X = leave_adjust_matrix(w)
    Y = leave_adjust_matrix(z)
    count = n - 1
    sx, sy = X.sum(axis=0), Y.sum(axis=0)
    sxx, syy, sxy = (X * X).sum(axis=0), (Y * Y).sum(axis=0), (X * Y).sum(axis=0)
    vx = (sxx - sx * sx / count) / (count - 1)
    vy = (syy - sy * sy / count) / (count - 1)
    cv = (sxy - sx * sy / count) / (count - 1)
    ax = vx - E[0, 0]
    ay = vy - E[1, 1]
    ac = cv - E[0, 1]
    leave = np.full(n, np.nan)
    ok = (ax > 1e-12) & (ay > 1e-12) & np.isfinite(ax) & np.isfinite(ay) & np.isfinite(ac)
    leave[ok] = ac[ok] / np.sqrt(ax[ok] * ay[ok])
    leave[ok] = np.clip(leave[ok], -0.999999, 0.999999)
    if not np.isfinite(theta) or not np.all(np.isfinite(leave)):
        return float(theta), np.full(n, np.nan)
    return float(theta), n * theta - (n - 1) * leave

def _el_lambda(g):
    g = np.asarray(g, dtype=float)
    if not np.all(np.isfinite(g)):
        return None
    if np.all(np.abs(g) < 1e-14):
        return 0.0
    gmin, gmax = float(g.min()), float(g.max())
    if gmin >= 0 or gmax <= 0:
        return None
    lo = -1.0 / gmax + 1e-12
    hi = -1.0 / gmin - 1e-12

    def score(lam):
        return float(np.sum(g / (1.0 + lam * g)))

    try:
        return float(brentq(score, lo, hi, maxiter=200))
    except ValueError:
        return None


def jel_log_ratio(pseudo, theta, adjusted=False):
    g = np.asarray(pseudo, dtype=float) - float(theta)
    if adjusted:
        a = math.log(max(len(g), 2)) / 2.0
        g = np.r_[g, -a * g.mean()]
    lam = _el_lambda(g)
    if lam is None:
        return np.inf
    den = 1.0 + lam * g
    if np.any(den <= 0):
        return np.inf
    return float(2.0 * np.log(den).sum())


def jel_ci(pseudo, level=0.95, adjusted=False, domain=(-0.999, 0.999)):
    pseudo = np.asarray(pseudo, dtype=float)
    if not np.all(np.isfinite(pseudo)):
        return np.array([np.nan, np.nan])
    crit = float(chi2.ppf(level, 1))
    lo, hi = map(float, domain)
    center = float(np.clip(np.mean(pseudo), lo + 1e-10, hi - 1e-10))

    def root_value(theta):
        lr = jel_log_ratio(pseudo, theta, adjusted)
        return (1e12 if not np.isfinite(lr) else lr) - crit

    if root_value(center) > 1e-6:
        return np.array([np.nan, np.nan])
    try:
        left = lo if root_value(lo) <= 0 else brentq(root_value, lo, center, maxiter=100)
        right = hi if root_value(hi) <= 0 else brentq(root_value, center, hi, maxiter=100)
        return np.array([left, right], dtype=float)
    except ValueError:
        grid = np.linspace(lo, hi, 201)
        vals = np.array([root_value(t) for t in grid])
        ok = vals <= 0
        if not np.any(ok):
            return np.array([np.nan, np.nan])
        ids = np.flatnonzero(ok)
        return np.array([grid[ids[0]], grid[ids[-1]]], dtype=float)
