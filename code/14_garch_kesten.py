"""Note 14: GARCH predicts its own tail exponent.

sigma^2_{t+1} = omega + (alpha z_t^2 + beta) sigma^2_t is a Kesten recursion.
Kesten-Goldie: P(sigma^2 > x) ~ x^-kappa with E[(alpha z^2 + beta)^kappa] = 1,
so returns r = sigma z have tail exponent 2 kappa (Breiman, if z is lighter).
Small-alpha approximation with persistence p = alpha + beta:
    2 kappa ~ 2 + 4 (1 - p) (-log p)/(1-p) p^2 / (alpha^2 Var z^2)
"""
import numpy as np
from numba import njit
from scipy.optimize import brentq
from scipy.stats import norm, t as student_t
from scipy.integrate import quad
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg)
    out.append(msg)


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


def approx_2kappa(alpha, p, var_z2):
    return 2 + 4 * p**2 * (-np.log(p)) / (alpha**2 * var_z2)


# ---------------------------------------------------------------- table
log("alpha  beta   persistence  2kappa (Gaussian z)  approx   2kappa (t6 z)  kurtosis finite?")
for alpha, beta in [(0.05, 0.94), (0.08, 0.90), (0.09, 0.90), (0.10, 0.89), (0.09, 0.905),
                    (0.12, 0.87), (0.07, 0.925), (0.15, 0.80)]:
    p = alpha + beta
    k_n = kappa(alpha, beta)
    k_t = kappa(alpha, beta, nu=6)
    kurt = (3 * alpha**2 + 2 * alpha * beta + beta**2) < 1
    log(f"{alpha:5.2f}  {beta:5.3f}  {p:8.3f}     {2*k_n:10.2f}        {approx_2kappa(alpha, p, 2.0):6.2f}   "
        f"{min(2*k_t, 6):8.2f}        {kurt}")

# ---------------------------------------------------------------- contour
alphas = np.linspace(0.02, 0.20, 37)
ps = np.linspace(0.90, 0.999, 45)
Z = np.full((len(ps), len(alphas)), np.nan)
for i, p in enumerate(ps):
    for j, a in enumerate(alphas):
        if p - a > 0:
            Z[i, j] = 2 * kappa(a, p - a)
fig, ax = plt.subplots(figsize=(7.8, 5.0))
levels = [2.5, 3, 3.5, 4, 5, 6, 8, 12]
cs = ax.contour(alphas, ps, np.clip(Z, 0, 30), levels=levels,
                colors=[style.BLUE_RAMP[min(6, k)] for k in range(len(levels))], linewidths=1.6)
ax.clabel(cs, fmt="%.1f", fontsize=8)
ax.fill_between([0.06, 0.11], 0.975, 0.995, color=style.SERIES[1], alpha=0.15, lw=0)
ax.annotate("typical daily\nequity-index fits", (0.062, 0.9765), fontsize=8, color=style.INK)
ax.set_xlabel("α (reaction to squared shocks)")
ax.set_ylabel("persistence α + β")
ax.set_title("Return tail exponent 2κ implied by Gaussian GARCH(1,1)")
style.save(fig, "14_contours.png")

# ---------------------------------------------------------------- simulation
@njit(cache=True)
def sim_garch(omega, alpha, beta, n, seed):
    np.random.seed(seed)
    s2 = omega / (1 - alpha - beta)
    r = np.empty(n)
    for t in range(n):
        z = np.random.randn()
        r[t] = np.sqrt(s2) * z
        s2 = omega + alpha * r[t] ** 2 + beta * s2
    return r


def hill(x, k):
    xs = np.sort(np.abs(x))[::-1]
    return 1 / np.mean(np.log(xs[:k] / xs[k]))


fig, ax = plt.subplots(figsize=(7.4, 4.4))
for i, (a, b) in enumerate([(0.09, 0.90), (0.09, 0.905), (0.05, 0.94)]):
    r = sim_garch(1e-6, a, b, 20_000_000, 3 + i)
    th = 2 * kappa(a, b)
    xs = np.sort(np.abs(r))[::-1]
    n = len(xs)
    surv = np.arange(1, n + 1) / n
    step = max(1, n // 20000)
    ax.loglog(xs[:n // 10:step] / np.std(r), surv[:n // 10:step], color=style.SERIES[i],
              label=f"α={a}, β={b}: 2κ = {th:.2f}")
    # reference slope
    x0 = xs[n // 1000] / np.std(r)
    xx = np.geomspace(x0, xs[0] / np.std(r), 20)
    ax.loglog(xx, 1e-3 * (xx / x0) ** (-th), color=style.SERIES[i], ls="--", lw=1)
    hills = {k: hill(r, k) for k in (2000, 20000, 200000)}
    log(f"sim alpha={a} beta={b}: theory 2kappa {th:.2f}; Hill k=2e3 {hills[2000]:.2f}, "
        f"k=2e4 {hills[20000]:.2f}, k=2e5 {hills[200000]:.2f}")
ax.set_xlabel("|return| / sd (log)")
ax.set_ylabel("P(|r| > x) (log)")
ax.set_title("Simulated GARCH tails (dashed: Kesten slope)")
ax.legend(fontsize=8)
style.plain_log(ax)
style.plain_log(ax, "x")
style.save(fig, "14_tails.png")

# ---------------------------------------------------------------- sensitivity
# How precisely does a GARCH fit pin down the tail? delta method on persistence.
a = 0.09
for p in [0.985, 0.99, 0.995]:
    k1, k2 = 2 * kappa(a, p - 0.0025 - a), 2 * kappa(a, p + 0.0025 - a)
    log(f"alpha={a}, persistence {p}±0.0025: 2kappa ranges {k2:.2f} .. {k1:.2f}")

open(style.os.path.join(style.FIG_DIR, "14_output.txt"), "w").write("\n".join(out) + "\n")
