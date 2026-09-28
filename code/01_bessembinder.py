"""Note 01: is Bessembinder's "4% of stocks create all the wealth" a theorem?

Simulates a CAPM world with *zero* alpha anywhere, lognormal idiosyncratic noise,
random listing dates and random lifetimes, and measures the same statistics
Bessembinder (2018) measured on CRSP.
"""
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
import style

style.apply()
rng = np.random.default_rng(1926)

# ---------------------------------------------------------------- analytic part
def share_needed(s, drift_excess=0.0):
    """Lognormal gross excess return X = exp(Y), Y~N(a, s^2) with a chosen so
    E[X] = exp(drift_excess). Benchmark gross = 1 (T-bills, in excess units).
    Returns the fraction of stocks (from the top) whose net creation (X-1)
    equals the aggregate net creation, i.e. the bottom rest nets to zero."""
    a = drift_excess - s * s / 2
    EX = np.exp(drift_excess)
    if EX <= 1:
        return np.nan

    def bottom_net(y0):  # E[(X-1) 1{Y<y0}]
        return EX * norm.cdf((y0 - a - s * s) / s) - norm.cdf((y0 - a) / s)

    # bottom_net -> 0 as y0 -> -inf, dips negative, then rises to E[X]-1 > 0;
    # we want the nontrivial root to the right of the minimum.
    ys = np.linspace(a - 10 * s, a + 10 * s + 5, 20001)
    vals = bottom_net(ys)
    i_min = np.argmin(vals)
    j = i_min + np.argmax(vals[i_min:] >= 0)
    y0 = brentq(bottom_net, ys[j - 1], ys[j])
    return 1 - norm.cdf((y0 - a) / s)


def frac_beat_bills(sigma, T, erp):
    """P(buy-and-hold beats T-bills) under GBM with arithmetic excess drift erp."""
    return norm.cdf((erp - sigma**2 / 2) * np.sqrt(T) / sigma)


# ---------------------------------------------------------------- simulation
MONTHS = 90 * 12            # 1926-2016
N = 25_000
RF_ANNUAL = 0.033
ERP_ANNUAL = 0.065          # arithmetic equity premium on the market
SIG_M = 0.19

rf = RF_ANNUAL / 12
# one market path (simple returns), repeated for a few seeds later
def simulate(seed, idio_scale=1.0, n=N):
    r = np.random.default_rng(seed)
    N = n
    mkt_ex = ERP_ANNUAL / 12 + SIG_M / np.sqrt(12) * r.standard_normal(MONTHS)
    birth = r.integers(0, MONTHS - 1, N)
    life = np.ceil(r.exponential(10.8 * 12, N)).astype(int)  # median ~7.5y
    death = np.minimum(birth + life, MONTHS)
    beta = np.clip(r.normal(1.0, 0.35, N), 0.2, 2.5)
    cap0 = np.exp(r.normal(0, 1.6, N))                         # starting size
    # annual idiosyncratic vol: smaller firms are noisier (~25% for giants,
    # ~80% for micro-caps), plus firm-level scatter
    sig_id = 0.45 * cap0 ** -0.12 * np.exp(0.25 * r.standard_normal(N))
    sig_id = idio_scale * sig_id

    t = np.arange(MONTHS)
    alive = (t[None, :] >= birth[:, None]) & (t[None, :] < death[:, None])
    s = sig_id[:, None] / np.sqrt(12)
    eps = r.standard_normal((N, MONTHS)) * s - s * s / 2
    gross = (1 + rf + beta[:, None] * mkt_ex[None, :]) * np.exp(eps)
    gross = np.where(alive, gross, 1.0)
    month_beat = (gross - 1 - rf > 0)[alive].mean()

    # buy-and-hold vs bills over each stock's life
    log_bh = np.log(gross).sum(1)
    n_months = death - birth
    log_bill = n_months * np.log1p(rf)
    beat = log_bh > log_bill

    # dollar wealth creation: sum_t cap_{t-1} (R_t - rf), capitalised to the end
    # of sample at the bill rate (as in Bessembinder).
    cap_path = cap0[:, None] * np.cumprod(gross, axis=1) / gross  # cap at t-1
    wc_t = np.where(alive, cap_path * (gross - 1 - rf), 0.0)
    to_end = (1 + rf) ** (MONTHS - 1 - t)
    wc = (wc_t * to_end[None, :]).sum(1)

    def top_share(w):
        order = np.argsort(-w)
        cum = np.cumsum(w[order])
        total = w.sum()
        k_all = int(np.argmax(cum >= total)) + 1 if total > 0 else N
        k_half = int(np.argmax(cum >= total / 2)) + 1
        return k_all / N, k_half, order

    frac_all, k_half, order = top_share(wc)
    # same statistic if every stock had started at $1 (removes size dispersion)
    wc_unit = np.exp(log_bh) - np.exp(log_bill)
    frac_all_unit, _, _ = top_share(wc_unit)
    return dict(
        frac_beat=beat.mean(),
        frac_month_beat=month_beat,
        frac_all=frac_all,
        frac_all_unit=frac_all_unit,
        k_half=k_half,
        median_life=np.median(n_months) / 12,
        wc_sorted=wc[order],
        log_rel=log_bh - log_bill,
        n_months=n_months,
    )


results = [simulate(s) for s in range(30)]
keys = ["frac_beat", "frac_month_beat", "frac_all", "frac_all_unit", "k_half",
        "median_life"]
print("baseline over 30 market histories: median [10th, 90th pct]")
for k in keys:
    v = np.array([res[k] for res in results], float)
    print(f"  {k:16s} {np.median(v):9.4f}  [{np.percentile(v, 10):.4f}, "
          f"{np.percentile(v, 90):.4f}]")

# the causal knob: scale idiosyncratic volatility and watch concentration move
scales = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25]
knob = {}
for sc in scales:
    runs = [simulate(100 + s, idio_scale=sc, n=10_000) for s in range(12)]
    knob[sc] = {k: np.median([r_[k] for r_ in runs]) for k in keys}
    print(f"idio x{sc:4.2f}: beat bills {knob[sc]['frac_beat']:.3f}  "
          f"all-wealth share {knob[sc]['frac_all']:.4f}  "
          f"(unit-start {knob[sc]['frac_all_unit']:.4f})")

# ---------------------------------------------------------------- figures
# Fig 1: analytic -- fraction of stocks accounting for all net wealth creation
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
fig.subplots_adjust(wspace=0.3)
ax = axes[0]
Ts = np.linspace(1, 40, 200)
for i, sig in enumerate([0.25, 0.40, 0.55]):
    ax.plot(Ts, 100 * frac_beat_bills(sig, Ts, 0.065), color=style.SERIES[i],
            label=f"σ = {sig:.0%}")
ax.axhline(50, color=style.NEUTRAL, lw=1, ls="--")
ax.set_xlabel("holding period T (years)")
ax.set_ylabel("% of stocks beating T-bills")
ax.set_title("Median stock loses to bills if σ²/2 > premium")
ax.legend(loc="upper right")

ax = axes[1]
ss = np.linspace(0.15, 3.0, 120)
# expected gross excess return held at exp(6.5% x 8y); only the spread varies
shares = [share_needed(s, drift_excess=0.065 * 8) for s in ss]
ax.plot(ss, 100 * np.array(shares), color=style.SERIES[0])
ax.set_yscale("log")
ax.set_xlabel("dispersion of log lifetime return, s = σ√T")
ax.set_ylabel("% of stocks that create all net wealth")
ax.set_title("Concentration depends on σ√T")
for sv in [0.4 * np.sqrt(7.5), 0.55 * np.sqrt(15)]:
    ax.axvline(sv, color=style.NEUTRAL, lw=1, ls=":")
ax.annotate("σ=40%, T=7.5y", (0.4 * np.sqrt(7.5), 60), xytext=(5, 0),
            textcoords="offset points", fontsize=8, color=style.INK_2)
ax.annotate("σ=55%, T=15y", (0.55 * np.sqrt(15), 30), xytext=(5, 0),
            textcoords="offset points", fontsize=8, color=style.INK_2)
style.plain_log(ax)
style.save(fig, "01_analytic.png")

# Fig 2: simulated Lorenz-style curve of wealth creation
fa = np.array([r_["frac_all"] for r_ in results])
res = results[int(np.argsort(fa)[len(fa) // 2])]       # the median history
wc = res["wc_sorted"]
fig, ax = plt.subplots(figsize=(6.5, 4.2))
x = 100 * np.arange(1, len(wc) + 1) / len(wc)
ax.plot(x, np.cumsum(wc) / wc.sum() * 100, color=style.SERIES[0])
ax.axhline(100, color=style.NEUTRAL, lw=1, ls="--")
ax.axvline(100 * res["frac_all"], color=style.SERIES[1], lw=1.2)
ax.annotate(f"top {100*res['frac_all']:.1f}% of stocks\n= 100% of net wealth",
            (100 * res["frac_all"], 20), xytext=(8, 0), textcoords="offset points",
            fontsize=9, color=style.INK)
ax.set_xscale("log")
style.plain_log(ax, "x")
ax.set_xlabel("top x% of stocks, ranked by dollar wealth creation (log)")
ax.set_ylabel("cumulative % of aggregate net wealth")
ax.set_title("Zero-alpha CAPM world, 25,000 simulated listings")
style.save(fig, "01_simulated_lorenz.png")

# Fig 2b: the knob
fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
axes[0].plot(scales, [100 * knob[sc]["frac_all"] for sc in scales], marker="o",
             ms=5, color=style.SERIES[0], label="dollar-weighted")
axes[0].plot(scales, [100 * knob[sc]["frac_all_unit"] for sc in scales],
             marker="o", ms=5, color=style.SERIES[1], label="every stock starts at $1")
axes[0].set_yscale("log")
axes[0].set_xlabel("idiosyncratic volatility, × baseline")
axes[0].set_ylabel("% of stocks creating all net wealth (log)")
axes[0].set_title("Turn up idiosyncratic noise → concentration")
style.plain_log(axes[0])
axes[0].legend()
axes[1].plot(scales, [100 * knob[sc]["frac_beat"] for sc in scales], marker="o",
             ms=5, color=style.SERIES[0])
axes[1].axhline(50, color=style.NEUTRAL, lw=1, ls="--")
axes[1].set_xlabel("idiosyncratic volatility, × baseline")
axes[1].set_ylabel("% of stocks beating T-bills over life")
axes[1].set_title("…and the typical stock loses to bills")
style.save(fig, "01_idio_knob.png")

# Fig 3: the "concentration tax" -- median growth of k-stock buy-and-hold
# portfolios versus k, for 20-year horizons, idiosyncratic noise only
T = 20
sig = 0.40
ks = [1, 2, 3, 5, 10, 20, 30, 50, 100, 300]
trials = 20_000
med = []
p_under = []
for k in ks:
    Y = sig * np.sqrt(T) * rng.standard_normal((trials, k)) - sig**2 * T / 2
    port = np.exp(Y).mean(1)                         # buy and hold, equal start
    med.append(np.median(np.log(port)) / T)
    p_under.append((port < 1).mean())
fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
axes[0].plot(ks, 100 * np.array(med), marker="o", ms=5, color=style.SERIES[0])
axes[0].set_xscale("log")
axes[0].set_xlabel("stocks held, k (log)")
axes[0].set_ylabel("median annual growth relative to index (%)")
axes[0].set_title("The typical k-stock investor lags the index")
axes[1].plot(ks, 100 * np.array(p_under), marker="o", ms=5, color=style.SERIES[1])
axes[1].axhline(50, color=style.NEUTRAL, lw=1, ls="--")
axes[1].set_xscale("log")
axes[1].set_xlabel("stocks held, k (log)")
axes[1].set_ylabel("% of portfolios behind the index after 20y")
axes[1].set_title("…and usually ends behind it")
style.save(fig, "01_concentration_tax.png")
for k, m, p in zip(ks, med, p_under):
    print(f"k={k:4d}  median growth vs index {100*m:6.2f}%/yr   P(behind) {p:.3f}")
