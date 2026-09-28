"""Note 24: 150 years of the S&P 500. What trends are real?

Uses Shiller's monthly S&P composite (1871-2023 with dividends/CPI; price to 2026).
Caveat built into every test: prices are MONTHLY AVERAGES of daily closes. For
a random walk this gives first differences with lag-1 autocorrelation 1/4 and
k-month change variance sigma^2 (k - 1/3), so the null variance ratio is not 1.
All nulls here are simulated with the same averaging.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import data
import style

style.apply()
rng = np.random.default_rng(24)
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


df = data.sp500_monthly()
tr = df["ret_real"].dropna()                        # real total return, 1871-2023
lr = np.log1p(tr)
yrs = len(tr) / 12
log(f"real total returns: {tr.index[0]:%Y-%m} to {tr.index[-1]:%Y-%m} ({yrs:.0f} years)")

# ------------------------------------------------------------ 1. growth rates
ar = 12 * tr.mean()
vol = np.sqrt(12) * tr.std()
geo = np.exp(12 * lr.mean()) - 1
log(f"arithmetic mean {100*ar:.2f}%/yr, vol {100*vol:.1f}%, geometric {100*geo:.2f}%/yr; "
    f"drag {100*(ar-np.log1p(geo)):.2f}% vs sigma^2/2 = {100*vol**2/2:.2f}%")
# note: averaging understates monthly vol by sqrt(2/3); annual non-overlapping
ann = np.exp(lr.groupby(lr.index.year).sum()) - 1
ann = ann[ann.index < 2024]
log(f"calendar-year real returns: mean {100*ann.mean():.2f}%, sd {100*ann.std():.1f}%, "
    f"geometric {100*(np.exp(np.log1p(ann).mean())-1):.2f}%, share of negative years {100*(ann<0).mean():.0f}%")

# by era
eras = [(1871, 1913), (1914, 1945), (1946, 1981), (1982, 2000), (2001, 2023)]
for a, b in eras:
    x = lr[(lr.index.year >= a) & (lr.index.year <= b)]
    log(f"   {a}-{b}: real geometric {100*(np.exp(12*x.mean())-1):5.2f}%/yr, "
        f"vol (x sqrt(3/2) averaging fix) {100*np.sqrt(12*1.5)*x.std():5.1f}%")

# ------------------------------------------------------------ 2. variance ratios
def vr_curve(logret, ks):
    s = np.concatenate([[0], np.cumsum(logret)])
    v1 = np.var(np.diff(s))
    return np.array([np.var(s[k:] - s[:-k]) / (k * v1) for k in ks])


ks = np.array([1, 2, 3, 6, 12, 24, 36, 60, 84, 120, 180, 240])
vr = vr_curve(lr.values, ks)


def averaged_rw(n_months, sigma, mu, r_, days=21):
    """Monthly averages of a daily random walk (in logs)."""
    daily = mu / days + sigma / np.sqrt(days) * r_.standard_normal(n_months * days)
    lp = np.cumsum(daily).reshape(n_months, days).mean(axis=1)
    return np.diff(lp)


sig_true = lr.std() * np.sqrt(1.5)
sims = np.array([vr_curve(averaged_rw(len(lr) + 1, sig_true, lr.mean(), rng), ks) for _ in range(1000)])
lo, med, hi = np.percentile(sims, [2.5, 50, 97.5], axis=0)
log("variance ratio VR(k) = Var(k-month change) / (k Var(1-month change)):")
for k, v, a, m_, b in zip(ks, vr, lo, med, hi):
    log(f"   k={k:4d} months: data {v:5.2f}   averaged random walk: median {m_:4.2f}, 95% band [{a:4.2f}, {b:4.2f}]")

# ------------------------------------------------------------ 3. autocorrelation after averaging
def acf(x, lags):
    x = x - x.mean()
    return np.array([np.dot(x[l:], x[:-l]) / np.dot(x, x) for l in lags])


lags = np.arange(1, 25)
ac = acf(lr.values, lags)
ac_sims = np.array([acf(averaged_rw(len(lr) + 1, sig_true, 0, rng), lags) for _ in range(500)])
ac_lo, ac_hi = np.percentile(ac_sims, [2.5, 97.5], axis=0)
log(f"lag-1 autocorrelation of monthly (averaged) returns: {ac[0]:.3f}; averaged-RW null band "
    f"[{ac_lo[0]:.3f}, {ac_hi[0]:.3f}] (theory 0.25 for continuous averaging)")
outside = [int(l) for l, a, b, c in zip(lags, ac, ac_lo, ac_hi) if a < b or a > c]
log(f"lags (1-24) outside the null band: {outside}; values: "
    + ", ".join(f"{l}:{a:+.3f}" for l, a in zip(lags, ac) if l in outside))
# joint test: sum of squared autocorrelations at lags 2-24 vs the null distribution
q = np.sum(ac[1:] ** 2)
q_null = np.sum(ac_sims[:, 1:] ** 2, axis=1)
log(f"portmanteau (lags 2-24): data {q:.4f}, null 95th pct {np.percentile(q_null, 95):.4f}, "
    f"p = {np.mean(q_null >= q):.3f}")
# heteroskedasticity-robust null: averaged random walk whose daily vol follows
# the realised monthly vol path (centred 13-month rolling sd, x sqrt(3/2))
vol_path = (lr.rolling(13, center=True, min_periods=6).std() * np.sqrt(1.5)).bfill().ffill().values


def averaged_rw_hetero(r_, days=21):
    n = len(vol_path) + 1
    vp = np.concatenate([[vol_path[0]], vol_path])
    daily = (vp[:, None] / np.sqrt(days)) * r_.standard_normal((n, days))
    lp = np.cumsum(daily.ravel()).reshape(n, days).mean(axis=1)
    return np.diff(lp)


ac_h = np.array([acf(averaged_rw_hetero(rng), lags) for _ in range(1000)])
hlo, hhi = np.percentile(ac_h, [2.5, 97.5], axis=0)
q_h = np.sum(ac_h[:, 1:] ** 2, axis=1)
out_h = [int(l) for l, a, b, c in zip(lags, ac, hlo, hhi) if a < b or a > c]
log(f"heteroskedastic null: lags outside 95% band {out_h}; portmanteau p = {np.mean(q_h >= q):.3f} "
    f"(95th pct {np.percentile(q_h, 95):.4f})")
vr_h = np.array([vr_curve(averaged_rw_hetero(rng), ks) for _ in range(1000)])
vhlo, vhmed, vhhi = np.percentile(vr_h, [2.5, 50, 97.5], axis=0)
log("full-sample VR vs heteroskedastic null: " + ", ".join(
    f"{k}m {v:.2f} [{a:.2f},{b:.2f}]" for k, v, a, b in zip(ks, vr, vhlo, vhhi) if k >= 12))
# post-1946 subsample autocorrelations
ac_post = acf(lr[lr.index.year >= 1946].values, lags)
log("post-1946 autocorrelations at lags 5,6,13,14,15,20,21: " + ", ".join(
    f"{l}:{ac_post[l-1]:+.3f}" for l in (5, 6, 13, 14, 15, 20, 21)))
ac_pre = acf(lr[lr.index.year < 1946].values, lags)
log("pre-1946 autocorrelations at the same lags: " + ", ".join(
    f"{l}:{ac_pre[l-1]:+.3f}" for l in (5, 6, 13, 14, 15, 20, 21)))

# post-war robustness of the variance ratio
post = lr[lr.index.year >= 1946].values
vr_post = vr_curve(post, ks)
sims_post = np.array([vr_curve(averaged_rw(len(post) + 1, post.std() * np.sqrt(1.5), post.mean(), rng), ks)
                      for _ in range(1000)])
plo, pmed, phi = np.percentile(sims_post, [2.5, 50, 97.5], axis=0)
log("post-1946 variance ratios: " + ", ".join(f"{k}m {v:.2f} [{a:.2f},{b:.2f}]"
                                                for k, v, a, b in zip(ks, vr_post, plo, phi) if k >= 12))

# ------------------------------------------------------------ 4. tails
def hill(x, frac):
    xs = np.sort(np.abs(x))[::-1]
    k = int(frac * len(xs))
    return 1 / np.mean(np.log(xs[:k] / xs[k]))


z = (lr - lr.rolling(1, min_periods=1).mean() * 0).values
lz = (lr.values - lr.mean())
log(f"monthly real log-returns: kurtosis {pd.Series(lz).kurt()+3:.1f}; Hill tail exponent (both tails) at top 2.5% / 5%: "
    f"{hill(lz, 0.025):.2f} / {hill(lz, 0.05):.2f}")
worst = lr.sort_values().head(5)
log("worst months: " + ", ".join(f"{d:%Y-%m} {100*(np.exp(v)-1):.0f}%" for d, v in worst.items()))

# ------------------------------------------------------------ 5. drawdowns
wealth = np.exp(lr.cumsum())
peak = wealth.cummax()
dd = wealth / peak - 1
log(f"max real drawdown {100*dd.min():.0f}% (trough {dd.idxmin():%Y-%m})")
# underwater spells
under = dd < -1e-9
spells, start = [], None
for d, u in under.items():
    if u and start is None:
        start = d
    if not u and start is not None:
        spells.append((start, d, (d - start).days / 365.25, dd[start:d].min()))
        start = None
spells = sorted(spells, key=lambda s: -s[2])[:5]
for s0, s1, length, depth in spells:
    log(f"   underwater {s0:%Y-%m} to {s1:%Y-%m}: {length:4.1f} years, depth {100*depth:.0f}%")

# ------------------------------------------------------------ figures
fig, ax = plt.subplots(figsize=(7.5, 4.0))
ax.fill_between(lags, hlo, hhi, color=style.NEUTRAL, alpha=0.35, lw=0,
                label="95% band: averaged random walk with the actual volatility history")
ax.bar(lags, ac, width=0.6, color=style.SERIES[0], label="S&P 500 monthly real returns")
ax.axhline(0, color=style.INK_2, lw=0.8)
ax.set_xlabel("lag (months)")
ax.set_ylabel("autocorrelation")
ax.set_title("The 0.26 at lag 1 is an averaging artefact")
ax.legend(fontsize=8)
style.save(fig, "24_acf.png")

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
axes[0].fill_between(ks / 12, vhlo, vhhi, color=style.NEUTRAL, alpha=0.35, lw=0,
                     label="95% band: averaged random walk, actual vol history")
axes[0].plot(ks / 12, vhmed, color=style.INK_2, lw=1, ls="--")
axes[0].plot(ks / 12, vr_post, color=style.SERIES[2], marker="s", ms=4, label="post-1946 only")
axes[0].plot(ks / 12, vr, color=style.SERIES[0], marker="o", ms=4, label="S&P 500 real total return, 1871–2023")
axes[0].axhline(1.0, color=style.NEUTRAL, lw=0.8, ls=":")
axes[0].set_xscale("log")
axes[0].set_xlabel("horizon (years, log)")
axes[0].set_ylabel("variance ratio VR(k)")
axes[0].set_title("Variance ratios vs the right null")
axes[0].legend(fontsize=8, loc="upper left")
style.plain_log(axes[0], "x")
axes[1].fill_between(dd.index, 100 * dd.values, 0, color=style.SERIES[1], alpha=0.5, lw=0)
axes[1].set_ylabel("real drawdown from previous peak (%)")
axes[1].set_title("Real drawdowns, 1871–2023")
style.save(fig, "24_variance_ratio_drawdowns.png")

open(style.os.path.join(style.FIG_DIR, "24_output.txt"), "w").write("\n".join(out) + "\n")
