"""Note 15: does fitting GARCH(1,1) to long-memory volatility invent power-law tails?

Simulate 30-year daily samples from
  (a) a true GARCH(1,1) (alpha=.09, beta=.90): implied 2kappa should be ~4.4;
  (b) lognormal SV with three AR(1) components (long-memory-like clustering,
      no power-law tail: all moments finite);
  (c) rough lognormal SV (log vol ~ fBM-like, H = 0.1).
Fit GARCH(1,1) by Gaussian MLE, compute the Kesten tail exponent of each fit,
and compare with the true tail behaviour.
"""
import numpy as np
from numba import njit
from scipy.optimize import minimize
from scipy.signal import lfilter
import matplotlib.pyplot as plt
import style

style.apply()
import kesten as kest
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


N = 7560          # 30 years of daily data
REPS = 60


@njit(cache=True)
def negloglik(params, r):
    omega, alpha, beta = params
    s2 = np.var(r)
    ll = 0.0
    for t in range(len(r)):
        ll += np.log(s2) + r[t] ** 2 / s2
        s2 = omega + alpha * r[t] ** 2 + beta * s2
    return 0.5 * ll


def fit_garch(r):
    v = np.var(r)

    def f(x):
        # x -> (omega, alpha, beta) with alpha, beta >= 0, alpha + beta < 1
        omega = v * np.exp(x[0])
        a = 1 / (1 + np.exp(-x[1]))
        pers = 1 / (1 + np.exp(-x[2]))
        alpha, beta = a * pers, (1 - a) * pers
        return negloglik(np.array([omega, alpha, beta]), r)

    best = None
    for x0 in ([-4.0, -2.3, 4.0], [-6.0, -2.0, 6.0], [-3.0, -1.5, 3.0]):
        res = minimize(f, np.array(x0), method="Nelder-Mead",
                       options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
        if best is None or res.fun < best.fun:
            best = res
    x = best.x
    a = 1 / (1 + np.exp(-x[1]))
    pers = 1 / (1 + np.exp(-x[2]))
    return v * np.exp(x[0]), a * pers, (1 - a) * pers


@njit(cache=True)
def sim_garch(omega, alpha, beta, n, seed):
    np.random.seed(seed)
    s2 = omega / (1 - alpha - beta)
    r = np.empty(n + 2000)
    for t in range(n + 2000):
        z = np.random.randn()
        r[t] = np.sqrt(s2) * z
        s2 = omega + alpha * r[t] ** 2 + beta * s2
    return r[2000:]


def sim_multi_sv(n, r_, comps=((0.5, 0.25), (0.97, 0.25), (0.998, 0.25))):
    """log sigma = sum of AR(1) components (persistence, stationary sd)."""
    burn = 5000
    ls = np.zeros(n + burn)
    for phi, sd in comps:
        e = r_.standard_normal(n + burn) * sd * np.sqrt(1 - phi**2)
        ls += lfilter([1], [1, -phi], e)
    ls = ls[burn:]
    sig = 0.01 * np.exp(ls - ls.var())                 # daily vol ~1%
    return sig * r_.standard_normal(n), ls


def fgn(n, H, r_):
    k = np.arange(n + 1)
    gam = 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))
    row = np.concatenate([gam, gam[-2:0:-1]])
    lam = np.clip(np.fft.fft(row).real, 0, None)
    Z = r_.standard_normal(len(row)) + 1j * r_.standard_normal(len(row))
    return np.fft.fft(np.sqrt(lam / len(row)) * Z).real[:n]


def sim_rough_sv(n, r_, H=0.1, nu=0.35, lam=1 / 2000):
    """log vol = mean-reverting integral of fGn increments (rough, H=0.1)."""
    burn = 20000
    inc = fgn(n + burn, H, r_) * nu * (1 / 1.0) ** H
    x = lfilter([1], [1, -(1 - lam)], inc)          # slow mean reversion
    x = x[burn:]
    sig = 0.01 * np.exp(x - x.var())
    return sig * r_.standard_normal(n), x


rows = {}
for name in ["true GARCH(1,1)", "3-component lognormal SV", "rough lognormal SV (H=0.1)"]:
    fits = []
    for k in range(REPS):
        r_ = np.random.default_rng(1000 + k)
        if name.startswith("true"):
            r = sim_garch(1e-6 * 0.01, 0.09, 0.90, N, 1000 + k)
        elif name.startswith("3-"):
            r, _ = sim_multi_sv(N, r_)
        else:
            r, _ = sim_rough_sv(N, r_)
        om, a, b = fit_garch(r)
        z2k = 2 * kest.kappa(a, b) if a + b < 0.99999 else 2.0
        fits.append((a, b, a + b, z2k))
    fits = np.array(fits)
    rows[name] = fits
    q = lambda c: np.percentile(fits[:, c], [10, 50, 90])
    log(f"{name:28s} alpha {q(0)[1]:.3f} [{q(0)[0]:.3f},{q(0)[2]:.3f}]  persistence {q(2)[1]:.4f} "
        f"[{q(2)[0]:.4f},{q(2)[2]:.4f}]  implied 2kappa {q(3)[1]:.2f} [{q(3)[0]:.2f},{q(3)[2]:.2f}]")

# true tails of the SV processes: local Hill slopes on a very long sample
def hill_curve(x, ks):
    xs = np.sort(np.abs(x))[::-1]
    return np.array([1 / np.mean(np.log(xs[:k] / xs[k])) for k in ks])


ks = np.unique(np.geomspace(100, 2_000_000, 30).astype(int))
r_big_m, _ = sim_multi_sv(20_000_000, np.random.default_rng(7))
r_big_r, _ = sim_rough_sv(20_000_000 // 4, np.random.default_rng(8))
r_big_g = sim_garch(1e-6 * 0.01, 0.09, 0.90, 20_000_000, 9)
hm = hill_curve(r_big_m, ks)
hr = hill_curve(r_big_r, ks[ks < len(r_big_r) // 5])
hg = hill_curve(r_big_g, ks)
for lab, rr, hh, kk in [("true GARCH", r_big_g, hg, ks), ("3-comp SV", r_big_m, hm, ks),
                        ("rough SV", r_big_r, hr, ks[ks < len(r_big_r) // 5])]:
    fr = kk / len(rr)
    vals = [hh[np.argmin(np.abs(fr - f))] for f in (0.01, 0.02, 0.05)]
    log(f"Hill at top 1% / 2% / 5% of the sample ({lab}): " + " / ".join(f"{v:.2f}" for v in vals))
log("Hill estimates vs number of tail observations k (true tail behaviour):")
for k_, a_, b_ in zip(ks[::5], hm[::5], hg[::5]):
    log(f"   k={k_:8d}: 3-comp SV {a_:5.2f}   true GARCH {b_:5.2f}")

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
names = list(rows.keys())
parts = axes[0].violinplot([rows[n][:, 3] for n in names], showmedians=True)
for i, pc in enumerate(parts["bodies"]):
    pc.set_facecolor(style.SERIES[i])
    pc.set_alpha(0.6)
for key in ("cmedians", "cbars", "cmins", "cmaxes"):
    parts[key].set_color(style.INK_2)
axes[0].axhline(2 * kest.kappa(0.09, 0.90), color=style.SERIES[0], ls=":", lw=1)
axes[0].set_xticks([1, 2, 3])
axes[0].set_xticklabels(["true GARCH", "3-component\nlognormal SV", "rough\nlognormal SV"], fontsize=9)
axes[0].set_ylabel("tail exponent implied by the GARCH fit")
axes[0].set_title("30-year samples: what the fitted GARCH claims")
axes[0].set_ylim(1.5, 10)
axes[1].semilogx(ks, hg, color=style.SERIES[0], marker="o", ms=3, label="true GARCH (Pareto, 2κ = 4.41)")
axes[1].semilogx(ks, hm, color=style.SERIES[1], marker="o", ms=3, label="3-component lognormal SV")
axes[1].semilogx(ks[ks < len(r_big_r) // 5], hr, color=style.SERIES[2], marker="o", ms=3,
                 label="rough lognormal SV")
axes[1].axhline(2 * kest.kappa(0.09, 0.90), color=style.SERIES[0], ls=":", lw=1)
axes[1].set_xlabel("tail sample size k (log)")
axes[1].set_ylabel("Hill estimate")
axes[1].set_title("What the tails actually do (20M observations)")
axes[1].legend(fontsize=8)
style.plain_log(axes[1], "x")
style.save(fig, "15_invented_tails.png")

open(style.os.path.join(style.FIG_DIR, "15_output.txt"), "w").write("\n".join(out) + "\n")
