"""Note 28: the variance risk premium and volatility targeting on real data.

1. VRP: VIX^2 (annualised implied variance for the next 30 calendar days) vs the
   realised variance of the S&P 500 over the next 21 trading days (1999-2018).
   Short-variance P&L, its skew, and the Kelly fraction from the exact
   empirical distribution vs the Gaussian (mean-variance) Kelly fraction.
2. Volatility targeting (note 09): estimate the risk-return exponent p from
   E[r | sigma] ~ sigma^p, then compare Sharpe ratios of exposure ~ sigma^-m
   for m = 0, 1, 2 with the lognormal prediction exp(-s^2 (m + p - 2)^2 / 2).
   Daily S&P 1999-2018 (vol = lagged VIX or EWMA) and monthly Fama-French
   market 1926-2018 (vol = lagged 6-month realised vol).
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import data
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


spx = data.sp500_daily()["Close"]
r = np.log(spx).diff().dropna()
vix = data.vix_daily()["close"].reindex(r.index).ffill()

# ------------------------------------------------------------ 1. VRP
H = 21
rv_fwd = (r**2)[::-1].rolling(H).sum()[::-1].shift(-1) * (252 / H)   # next 21 days, annualised
iv = (vix / 100) ** 2
d = pd.DataFrame({"iv": iv, "rv": rv_fwd}).dropna()
vrp = d["iv"] - d["rv"]
log(f"VRP sample {d.index[0]:%Y-%m} to {d.index[-1]:%Y-%m}: mean implied var {d['iv'].mean():.4f} "
    f"(vol {np.sqrt(d['iv'].mean()):.3f}), mean realised var {d['rv'].mean():.4f} "
    f"(vol {np.sqrt(d['rv'].mean()):.3f}); implied > realised on {100*(vrp>0).mean():.0f}% of days")
# non-overlapping monthly short-variance returns, per unit of implied variance
mon = d.iloc[::H]
ret_sv = (mon["iv"] - mon["rv"]) / mon["iv"]            # payoff of a var swap scaled by its strike
log(f"short variance swap (strike-scaled, {len(ret_sv)} non-overlapping months): mean {ret_sv.mean():+.3f}, "
    f"sd {ret_sv.std():.3f}, Sharpe (annualised) {ret_sv.mean()/ret_sv.std()*np.sqrt(12):.2f}, "
    f"skew {ret_sv.skew():.2f}, worst {ret_sv.min():+.2f} ({ret_sv.idxmin():%Y-%m})")
worst = ret_sv.sort_values().head(5)
log("worst months: " + ", ".join(f"{k:%Y-%m} {v:+.2f}" for k, v in worst.items()))


def kelly_exact(x):
    """argmax_f E log(1 + f x) over the empirical distribution."""
    lo_bound = -1 / x.max() + 1e-6 if x.max() > 0 else -10
    hi_bound = -1 / x.min() - 1e-6 if x.min() < 0 else 10
    res = minimize_scalar(lambda f: -np.mean(np.log1p(f * x)), bounds=(max(lo_bound, -10), min(hi_bound, 10)),
                          method="bounded")
    return res.x, -res.fun


x = ret_sv.values
f_gauss = x.mean() / x.var()
f_exact, g_exact = kelly_exact(x)
g_at_gauss = np.mean(np.log1p(f_gauss * x)) if np.all(1 + f_gauss * x > 0) else -np.inf
log(f"Kelly leverage on strike-scaled short variance: mean-variance {f_gauss:.2f} vs exact empirical "
    f"{f_exact:.2f} (ruin bound 1/max loss = {1/(-x.min()):.2f}); growth at exact Kelly {12*g_exact:+.3f}/yr, "
    f"at mean-variance Kelly {'ruin' if not np.isfinite(g_at_gauss) else f'{12*g_at_gauss:+.3f}/yr'}")

# ------------------------------------------------------------ 2. vol targeting, daily S&P
def sharpe(x):
    return x.mean() / x.std() * np.sqrt(252)


ew = np.sqrt((r**2).ewm(halflife=21, min_periods=60).mean() * 252).shift(1)
vol_vix = (vix / 100).shift(1)
dd = pd.DataFrame({"r": r, "ew": ew, "vix": vol_vix}).dropna()
log(f"daily S&P {dd.index[0]:%Y}-{dd.index[-1]:%Y}:")
for vname in ["ew", "vix"]:
    sig = dd[vname]
    s_log = np.log(sig).std()
    # risk-return exponent: regress next-day return on log sigma: E[r] ~ sigma^p means
    # E[r / sigma^p] constant; estimate p by maximising fit of E[r|bucket] ~ c sigma^p
    q = pd.qcut(sig, 5, labels=False)
    mu_b = dd.groupby(q)["r"].mean() * 252
    sig_b = sig.groupby(q).mean()
    line = ", ".join(f"σ≈{s_:.2f}: μ={m_:+.2f}" for s_, m_ in zip(sig_b, mu_b))
    base = sharpe(dd["r"])
    res_m = {m: sharpe(dd["r"] * sig ** (-m)) for m in (0, 1, 2)}
    log(f"  vol proxy {vname}: sd(log sigma) = {s_log:.2f}; annualised mean return by vol quintile: {line}")
    log(f"  Sharpe: constant {res_m[0]:.3f}, vol-targeted {res_m[1]:.3f}, variance-targeted {res_m[2]:.3f}")

# ------------------------------------------------------------ monthly Fama-French market, 1926-2018
ff = data.french_monthly()
mk = ff["Mkt-RF"]
rv6 = np.sqrt((mk**2).rolling(6).mean() * 12).shift(1)
dm = pd.DataFrame({"r": mk, "s": rv6}).dropna()
s_log = np.log(dm["s"]).std()
q = pd.qcut(dm["s"], 5, labels=False)
mu_b = dm.groupby(q)["r"].mean() * 12
sig_b = dm.groupby(q)["s"].mean()
# p from log-log slope of mean excess return on vol across quintiles (means are positive)
p_hat = np.polyfit(np.log(sig_b), np.log(np.maximum(mu_b, 1e-4)), 1)[0]
sh = {m: (dm["r"] * dm["s"] ** (-m)).mean() / (dm["r"] * dm["s"] ** (-m)).std() * np.sqrt(12) for m in (0, 1, 2)}
pred = {m: np.exp(-s_log**2 * (m + p_hat - 2) ** 2 / 2) for m in (0, 1, 2)}
log(f"monthly FF market {dm.index[0]:%Y}-{dm.index[-1]:%Y}: sd(log sigma) = {s_log:.2f}; mean excess return by "
    f"vol quintile: " + ", ".join(f"σ≈{a:.2f}: {b:+.3f}" for a, b in zip(sig_b, mu_b)))
log(f"  estimated risk-return exponent p ≈ {p_hat:.2f}  (note 09: vol targeting helps iff p < 1.5)")
log(f"  Sharpe: constant {sh[0]:.3f}, vol-targeted {sh[1]:.3f}, variance-targeted {sh[2]:.3f}; "
    f"ratio vol/const {sh[1]/sh[0]:.2f} (lognormal theory {pred[1]/pred[0]:.2f}), "
    f"var/const {sh[2]/sh[0]:.2f} (theory {pred[2]/pred[0]:.2f})")

# ------------------------------------------------------------ figures
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3))
fig.subplots_adjust(wspace=0.25)
axes[0].plot(d.index, 100 * np.sqrt(d["iv"]), color=style.SERIES[0], lw=0.8, label="VIX (implied vol)")
axes[0].plot(d.index, 100 * np.sqrt(d["rv"]), color=style.SERIES[1], lw=0.8, label="realised vol, next 21 days")
axes[0].set_yscale("log")
axes[0].set_ylabel("annualised vol (%, log)")
axes[0].set_title(f"Implied above realised on {100*(vrp>0).mean():.0f}% of days")
axes[0].legend(fontsize=8)
axes[0].set_yticks([5, 10, 20, 40, 80])
style.plain_log(axes[0])
axes[1].hist(ret_sv, bins=40, color=style.SERIES[0], alpha=0.85)
axes[1].axvline(0, color=style.INK_2, lw=0.8)
axes[1].set_xlabel("monthly short-variance return (per unit of strike)")
axes[1].set_ylabel("months")
axes[1].set_title(f"Picking up nickels: skew {ret_sv.skew():.1f}")
style.save(fig, "28_vrp.png")

open(style.os.path.join(style.FIG_DIR, "28_output.txt"), "w").write("\n".join(out) + "\n")
