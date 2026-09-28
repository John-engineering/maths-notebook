"""Kesten tail exponent of GARCH(1,1): solve E[(alpha z^2 + beta)^kappa] = 1."""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm, t as student_t
from scipy.integrate import quad


def moment_A(kappa, alpha, beta, nu=None):
    """E[(alpha z^2 + beta)^kappa] for standard normal or unit-variance t_nu z."""
    if nu is None:
        f = lambda z: np.exp(kappa * np.log(alpha * z * z + beta) + norm.logpdf(z))
    else:
        s = np.sqrt((nu - 2) / nu)
        f = lambda z: np.exp(kappa * np.log(alpha * z * z + beta)
                             + student_t.logpdf(z / s, nu) - np.log(s))
    return 2 * quad(f, 0, np.inf, limit=200)[0]


def kappa(alpha, beta, nu=None):
    if alpha + beta >= 1:
        return 1.0 if abs(alpha + beta - 1) < 1e-12 else np.nan
    hi = 30.0 if nu is None else nu / 2 - 1e-6
    g = lambda k: moment_A(k, alpha, beta, nu) - 1
    if g(hi) < 0:
        return np.inf
    return brentq(g, 1e-6, hi)
