"""Note 05: Kelly with an estimated edge.

1. One asset, mean estimated from T years: E[growth] * T = k t^2 - k^2 (t^2+1)/2,
   t = SR sqrt(T). Full Kelly needs t > 1; optimal fraction t^2/(1+t^2).
2. N assets: plug-in Kelly growth = (SR^2 - N/T)/2; James-Stein shrinkage.
3. Selection: pick the best of M backtests. Tweedie's formula, fitted on the
   cross-section of backtests, gives the posterior-mean edge.
"""
import numpy as np
from scipy.stats import norm, gaussian_kde
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(5)
out = []


def log(msg):
    print(msg)
    out.append(msg)


# ------------------------------------------------------------------ 1. one asset
# Units: time in multiples of the sample length T, so growth*T depends only on t.
def growth_T(k, t, that):
    """Realised growth*T of betting f = k * mu_hat / sigma^2 when truth is t."""
    return k * that * t - 0.5 * (k * that) ** 2


ts = np.linspace(0, 4, 81)
Z = rng.standard_normal(200_000)
curves = {"full Kelly (plug-in)": [], "half Kelly": [], "shrink 1 − 1/t̂² (plug-in)": [],
          "oracle t²/(1+t²)": []}
for t in ts:
    th = t + Z
    curves["full Kelly (plug-in)"].append(growth_T(1.0, t, th).mean())
    curves["half Kelly"].append(growth_T(0.5, t, th).mean())
    kjs = np.clip(1 - 1 / th**2, 0, None)
    curves["shrink 1 − 1/t̂² (plug-in)"].append(growth_T(kjs, t, th).mean())
    curves["oracle t²/(1+t²)"].append(growth_T(t**2 / (1 + t**2), t, th).mean())
best = 0.5 * ts**2                                      # known-mu Kelly
fig, ax = plt.subplots(figsize=(7.5, 4.4))
ax.plot(ts, best, color=style.NEUTRAL, ls="--", lw=1.2, label="known edge: t²/2")
for i, (k_, v) in enumerate(curves.items()):
    ax.plot(ts, v, color=style.SERIES[i], label=k_)
ax.axhline(0, color=style.NEUTRAL, lw=1)
ax.set_ylim(-1.2, 4.5)
ax.set_xlabel("true t-statistic of the edge, t = SR·√T")
ax.set_ylabel("expected growth × T (above cash)")
ax.set_title("Betting an estimated edge: full Kelly needs t > 1")
ax.legend(fontsize=8, loc="upper left")
style.save(fig, "05_one_asset.png")
for t in [0.5, 1, 1.5, 2, 3]:
    i = np.argmin(np.abs(ts - t))
    log(f"1  t={t:3.1f}: full {curves['full Kelly (plug-in)'][i]:+.3f} (theory {0.5*(t*t-1):+.3f}), "
        f"half {curves['half Kelly'][i]:+.3f} (theory {(3*t*t-1)/8:+.3f}), "
        f"oracle {curves['oracle t²/(1+t²)'][i]:+.3f} (theory {0.5*t**4/(1+t*t):+.3f}), "
        f"JS {curves['shrink 1 − 1/t̂² (plug-in)'][i]:+.3f}")

# ------------------------------------------------------------------ 2. N assets
# One-factor covariance, alphas drawn so the tangency Sharpe is SR_TAN.
def universe(N, SR_TAN, r_):
    beta = r_.uniform(0.6, 1.4, N)
    sig_m, sig_i = 0.16, 0.25
    Sig = sig_m**2 * np.outer(beta, beta) + sig_i**2 * np.eye(N)
    mu = beta * 0.04 + r_.normal(0, 0.02, N)
    sr2 = mu @ np.linalg.solve(Sig, mu)
    mu *= SR_TAN / np.sqrt(sr2)
    return mu, Sig


T_YEARS = 20
SR_TAN = 0.8
Ns = [2, 5, 10, 20, 50, 100, 200]
rows = []
for N in Ns:
    mu, Sig = universe(N, SR_TAN, np.random.default_rng(N))
    Sinv = np.linalg.inv(Sig)
    L = np.linalg.cholesky(Sig / T_YEARS)
    g_full, g_js, g_orc, g_ew = [], [], [], []
    k_orc = SR_TAN**2 * T_YEARS / (SR_TAN**2 * T_YEARS + N)
    ew = np.ones(N) / N
    for _ in range(4000):
        mh = mu + L @ rng.standard_normal(N)
        f = Sinv @ mh
        G = lambda f_: f_ @ mu - 0.5 * f_ @ Sig @ f_
        g_full.append(G(f))
        q = T_YEARS * mh @ Sinv @ mh                  # = sum of squared t-stats
        kjs = max(0.0, 1 - (N - 2) / q) if N > 2 else max(0.0, 1 - 1 / q)
        g_js.append(G(kjs * f))
        g_orc.append(G(k_orc * f))
        # 1/N portfolio, levered by its own estimated Kelly (1-D problem)
        m_ew = ew @ mh
        v_ew = ew @ Sig @ ew
        g_ew.append(G(ew * (m_ew / v_ew)))
    rows.append((N, np.mean(g_full), np.mean(g_js), np.mean(g_orc), np.mean(g_ew)))
    log(f"2  N={N:3d}: plug-in {np.mean(g_full):+.4f} (theory {0.5*(SR_TAN**2 - N/T_YEARS):+.4f}), "
        f"James-Stein {np.mean(g_js):+.4f}, oracle-scaled {np.mean(g_orc):+.4f} "
        f"(k={k_orc:.2f}), levered 1/N {np.mean(g_ew):+.4f}; known-mu max {0.5*SR_TAN**2:.3f}")
rows = np.array(rows)
fig, ax = plt.subplots(figsize=(7.5, 4.4))
labels = ["plug-in Kelly", "James–Stein shrunk", "oracle-scaled", "1/N, levered"]
for j in range(4):
    ax.plot(rows[:, 0], 100 * rows[:, j + 1], marker="o", ms=4, color=style.SERIES[j],
            label=labels[j])
ax.axhline(100 * 0.5 * SR_TAN**2, color=style.NEUTRAL, ls="--", lw=1.2)
ax.annotate("known means", (2, 100 * 0.5 * SR_TAN**2), xytext=(0, 4),
            textcoords="offset points", fontsize=8, color=style.INK_2)
ax.axhline(0, color=style.NEUTRAL, lw=1)
ax.set_xscale("log")
ax.set_yscale("symlog", linthresh=10)
ax.set_xlabel("number of assets N (log)")
ax.set_ylabel("expected growth above cash (%/yr)")
ax.set_title(f"20 years of data, tangency Sharpe {SR_TAN}: optimisation vs 1/N")
ax.legend(fontsize=8, loc="lower left")
style.plain_log(ax, "x")
style.save(fig, "05_n_assets.png")

# ------------------------------------------------------------------ 3. selection
# A research shop tests M strategies and bets Kelly on the best backtest.
# How much should it shrink? Three worlds, differing only in how common and how
# strong real edges are.
from scipy.optimize import minimize

M = 2000
T_BT = 10
REGIMES = {
    "no real edges": (0.00, 0.0),
    "2% weak edges": (0.02, 0.15),
    "5% strong edges": (0.05, 0.40),
}


def research_round(r_, frac, mean_sr):
    real = r_.random(M) < frac
    sr = np.where(real, r_.exponential(mean_sr if mean_sr > 0 else 1, M), 0.0)
    t_true = sr * np.sqrt(T_BT)
    t_hat = t_true + r_.standard_normal(M)
    return t_true, t_hat


def tweedie_lindsey(t_hat, x_eval, degree=5, bins=80):
    """Efron's Tweedie correction E[t|t_hat] = t_hat + (log f)'(t_hat), with the
    marginal density f fitted by Lindsey's method (Poisson regression of
    histogram counts on a polynomial)."""
    counts, edges = np.histogram(t_hat, bins=bins)
    xc = 0.5 * (edges[1:] + edges[:-1])
    loc, sc = xc.mean(), xc.std()
    X = np.vander((xc - loc) / sc, degree + 1, increasing=True)

    def nll(b):
        eta = X @ b
        return np.sum(np.exp(eta) - counts * eta)

    def grad(b):
        return X.T @ (np.exp(X @ b) - counts)

    b0 = np.zeros(degree + 1)
    b0[0] = np.log(counts.mean() + 1)
    b = minimize(nll, b0, jac=grad, method="BFGS").x
    u = (np.asarray(x_eval) - loc) / sc
    dpoly = sum(k * b[k] * u ** (k - 1) for k in range(1, degree + 1)) / sc
    return np.asarray(x_eval) + dpoly


def npmle_posterior_mean(t_hat, x_eval, grid=np.linspace(-2, 12, 141), iters=600):
    """Kiefer-Wolfowitz nonparametric MLE of the prior on true t (a discrete
    distribution on a grid), fitted by EM on binned data; returns E[t | t_hat]."""
    counts, edges = np.histogram(t_hat, bins=np.arange(t_hat.min() - 0.05,
                                                         t_hat.max() + 0.1, 0.05))
    xc = 0.5 * (edges[1:] + edges[:-1])
    keep = counts > 0
    xc, counts = xc[keep], counts[keep]
    Lik = norm.pdf(xc[:, None] - grid[None, :])
    w = np.full(len(grid), 1 / len(grid))
    for _ in range(iters):
        P = Lik * w
        P /= P.sum(1, keepdims=True)
        w = (counts[:, None] * P).sum(0) / counts.sum()
    x_eval = np.asarray(x_eval)
    Le = norm.pdf(x_eval[:, None] - grid[None, :]) * w
    return (Le @ grid) / Le.sum(1)


def tweedie_kde(t_hat, x_eval):
    kde = gaussian_kde(t_hat, bw_method=0.25)
    h = 1e-3
    x_eval = np.asarray(x_eval)
    dlog = (np.log(kde(x_eval + h)) - np.log(kde(x_eval - h))) / (2 * h)
    return x_eval + dlog


# illustrative round for the figure: the weak-edge world
t_true, t_hat = research_round(np.random.default_rng(99), *REGIMES["2% weak edges"])
grid = np.linspace(t_hat.min(), t_hat.max(), 300)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
fig.subplots_adjust(wspace=0.25)
axes[0].hist(t_hat, bins=60, color=style.SERIES[0], alpha=0.85)
axes[0].set_xlabel("backtest t-statistic")
axes[0].set_ylabel("number of strategies")
axes[0].set_title(f"{M} backtests, 2% with a (weak) real edge")
axes[1].plot(grid, grid, color=style.NEUTRAL, ls="--", lw=1.2, label="naive: believe the backtest")
axes[1].plot(grid, tweedie_lindsey(t_hat, grid), color=style.SERIES[1],
             label="Tweedie (f-modelling)")
axes[1].plot(grid, npmle_posterior_mean(t_hat, grid), color=style.SERIES[2],
             label="NPMLE prior (g-modelling)")
axes[1].scatter(t_hat, t_true, s=5, alpha=0.3, color=style.SERIES[0], edgecolors="none",
                label="(t̂, true t) for each strategy")
axes[1].set_xlabel("backtest t-statistic t̂")
axes[1].set_ylabel("true t")
axes[1].set_title("Selection bias, corrected without knowing the truth")
axes[1].legend(fontsize=8, loc="upper left")
style.save(fig, "05_tweedie.png")

# many rounds per regime: realised growth of betting on the top strategy
for regime, params in REGIMES.items():
    res = {"naive": [], "half Kelly": [], "single JS": [], "Tweedie (KDE)": [],
           "Tweedie (Lindsey)": [], "NPMLE": [], "oracle": []}
    tops = []
    for rep in range(300):
        r_ = np.random.default_rng(1000 + rep)
        t_true, t_hat = research_round(r_, *params)
        i = np.argmax(t_hat)
        th, tt = t_hat[i], t_true[i]
        tops.append((th, tt))
        post_k = tweedie_kde(t_hat, [th])[0]
        post_l = tweedie_lindsey(t_hat, [th])[0]
        post_n = npmle_posterior_mean(t_hat, [th])[0]
        for name, k in [("naive", 1.0), ("half Kelly", 0.5),
                        ("single JS", max(0, 1 - 1 / th**2)),
                        ("Tweedie (KDE)", max(0.0, post_k) / th),
                        ("Tweedie (Lindsey)", max(0.0, post_l) / th),
                        ("NPMLE", max(0.0, post_n) / th),
                        ("oracle", tt / th)]:
            res[name].append(growth_T(k, tt, th))
    tops = np.array(tops)
    log(f"3  [{regime}] best backtest t̂ averages {tops[:,0].mean():.2f}, its true t {tops[:,1].mean():.2f}")
    for name, v in res.items():
        v = np.array(v)
        log(f"3  [{regime}] {name:17s} mean growth*T {v.mean():+7.3f}   P(growth<0) {np.mean(v < 0):.2f}")

open(style.os.path.join(style.FIG_DIR, "05_output.txt"), "w").write("\n".join(out) + "\n")
