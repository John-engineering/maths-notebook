"""Note 25: does valuation (CAPE) predict 10-year returns, and how much of that is real?

1. In-sample: regress subsequent 10-year annualised real total return on
   log(CAPE). Report slope and R^2.
2. Null: random-walk prices (bootstrapped real returns), CAPE REBUILT from the
   simulated prices and the actual earnings history, same regression. This
   keeps the mechanical link (price in CAPE's numerator) that biases such
   regressions (Stambaugh), so it is the fair comparison.
3. Out-of-sample (Goyal-Welch): expanding-window forecasts vs the historical mean.
4. The level of CAPE has drifted up: era means.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import data
import style

style.apply()
rng = np.random.default_rng(25)
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


df = data.sp500_monthly()
df = df[df.index <= "2023-06-01"].copy()
H = 120                                              # 10 years
lr = np.log1p(df["ret_real"])                        # real total log return (month t)
fut = lr[::-1].rolling(H).sum()[::-1].shift(-1) / (H / 12)   # next 10y, annualised log
cape = df["cape"]
d = pd.DataFrame({"fut": fut, "lcape": np.log(cape)}).dropna()
log(f"sample: {d.index[0]:%Y-%m} to {d.index[-1]:%Y-%m} start dates ({len(d)} months, "
    f"~{len(d)/H:.0f} non-overlapping decades)")


def ols(x, y):
    X = np.column_stack([np.ones(len(x)), x])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ b
    return b, 1 - r.var() / y.var()


b, r2 = ols(d["lcape"].values, d["fut"].values)
log(f"in-sample: fut10 = {b[0]:.3f} + {b[1]:.3f} log(CAPE); R^2 = {r2:.3f}")
# non-overlapping decades
nono = d.iloc[::H]
bn, r2n = ols(nono["lcape"].values, nono["fut"].values)
log(f"non-overlapping decades ({len(nono)} points): slope {bn[1]:.3f}, R^2 {r2n:.3f}")

# ------------------------------------------------------------ null: rebuilt CAPE
# real price index and 10y average real earnings, as Shiller constructs CAPE
real_price = df["real_price"]
real_earn = df["real_earn"] if "real_earn" in df else df["Real Earnings"]
e10 = real_earn.rolling(H).mean()
ret_pool = lr.dropna().values
ret_idx = lr.dropna().index


def null_r2(n_sims=1000):
    r2s, slopes = [], []
    base = real_price.loc[ret_idx[0]]
    for _ in range(n_sims):
        # stationary block bootstrap of real total returns (12-month blocks)
        n = len(ret_pool)
        idx = []
        while len(idx) < n:
            s = rng.integers(0, n - 12)
            idx.extend(range(s, s + 12))
        sim = ret_pool[np.array(idx[:n])]
        sim = sim - sim.mean() + ret_pool.mean()
        # price-only part: remove the average dividend yield so CAPE levels are comparable
        price_lr = sim - np.log1p(df["div"] / df["price"] / 12).reindex(ret_idx).fillna(0).values
        p = base * np.exp(np.cumsum(price_lr))
        p = pd.Series(p, index=ret_idx)
        cape_sim = p / e10.reindex(ret_idx)
        fut_sim = pd.Series(sim, index=ret_idx)[::-1].rolling(H).sum()[::-1].shift(-1) / (H / 12)
        dd = pd.DataFrame({"fut": fut_sim, "lcape": np.log(cape_sim)}).dropna()
        bb, rr = ols(dd["lcape"].values, dd["fut"].values)
        r2s.append(rr)
        slopes.append(bb[1])
    return np.array(r2s), np.array(slopes)


r2_null, sl_null = null_r2()
log(f"null (random-walk prices, CAPE rebuilt from actual earnings): R^2 median {np.median(r2_null):.3f}, "
    f"95th pct {np.percentile(r2_null, 95):.3f}; P(R^2 >= data) = {np.mean(r2_null >= r2):.3f}; "
    f"P(slope <= data) = {np.mean(sl_null <= b[1]):.3f}")

log(f"null slope distribution: median {np.median(sl_null):.3f}, 5th pct {np.percentile(sl_null, 5):.3f} "
    f"(data {b[1]:.3f}); the null median is not zero because of Stambaugh bias")

# ------------------------------------------------------------ out-of-sample
start = d.index.get_loc(d.index[d.index >= "1920-01-01"][0])
preds, bench, actual, dates = [], [], [], []
for i in range(start, len(d)):
    # only use pairs whose 10y outcome was known at time t: start dates <= t - H
    train = d.iloc[: max(0, i - H + 1)]
    if len(train) < 240:
        continue
    bb, _ = ols(train["lcape"].values, train["fut"].values)
    preds.append(bb[0] + bb[1] * d["lcape"].iloc[i])
    bench.append(train["fut"].mean())
    actual.append(d["fut"].iloc[i])
    dates.append(d.index[i])
preds, bench, actual = map(np.array, (preds, bench, actual))
oos = 1 - np.sum((actual - preds) ** 2) / np.sum((actual - bench) ** 2)
log(f"out-of-sample R^2 vs expanding historical mean ({dates[0]:%Y}-{dates[-1]:%Y} forecasts): {oos:.3f}")
for a, bnd in [("1920", "1959"), ("1960", "1989"), ("1990", "2013")]:
    m = (pd.DatetimeIndex(dates) >= a) & (pd.DatetimeIndex(dates) <= bnd + "-12-31")
    o = 1 - np.sum((actual[m] - preds[m]) ** 2) / np.sum((actual[m] - bench[m]) ** 2)
    bias = np.mean(actual[m] - preds[m])
    log(f"   forecasts made {a}-{bnd}: OOS R^2 {o:+.3f}, mean error (actual - forecast) {100*bias:+.2f}%/yr")

# ------------------------------------------------------------ drift in the level of CAPE
for a, bnd in [(1881, 1919), (1920, 1959), (1960, 1989), (1990, 2023)]:
    c = cape[(cape.index.year >= a) & (cape.index.year <= bnd)]
    log(f"CAPE {a}-{bnd}: mean {c.mean():.1f}, median {c.median():.1f}")
log(f"last CAPE in data ({cape.dropna().index[-1]:%Y-%m}): {cape.dropna().iloc[-1]:.1f}; "
    f"model-implied next-10y real return {100*(np.exp(b[0] + b[1]*np.log(cape.dropna().iloc[-1]))-1):.1f}%/yr")

full = data.sp500_monthly()
p0, p1 = full.loc["2023-06-01", "price"], full["price"].iloc[-1]
yrs_after = (full.index[-1] - pd.Timestamp("2023-06-01")).days / 365.25
log(f"what happened next (price only, nominal): 2023-06 to {full.index[-1]:%Y-%m}: "
    f"{100*((p1/p0)**(1/yrs_after)-1):.1f}%/yr over {yrs_after:.1f} years")

# ------------------------------------------------------------ figure
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
x = d["lcape"].values
yv = 100 * (np.exp(d["fut"].values) - 1)
yr = d.index.year
sc = axes[0].scatter(np.exp(x), yv, c=yr, cmap="Blues", vmin=1850, s=6, alpha=0.7, edgecolors="none")
xx = np.linspace(x.min(), x.max(), 50)
axes[0].plot(np.exp(xx), 100 * (np.exp(b[0] + b[1] * xx) - 1), color=style.INK, lw=1.5)
axes[0].set_xscale("log")
axes[0].set_xlabel("CAPE at start (log scale)")
axes[0].set_ylabel("next 10-year real return (%/yr)")
axes[0].set_title(f"In sample: R² = {r2:.2f}")
cb = fig.colorbar(sc, ax=axes[0])
cb.set_label("start year")
style.plain_log(axes[0], "x")
axes[1].hist(r2_null, bins=40, color=style.NEUTRAL, alpha=0.8, label="random-walk prices, rebuilt CAPE")
axes[1].axvline(r2, color=style.SERIES[1], lw=2, label=f"actual R² = {r2:.2f}")
axes[1].set_xlabel("R² of the 10-year return regression")
axes[1].set_ylabel("simulations")
axes[1].set_title(f"How unusual is that under a random walk? p = {np.mean(r2_null >= r2):.2f}")
axes[1].legend(fontsize=8)
style.save(fig, "25_cape.png")

fig, ax = plt.subplots(figsize=(9, 3.8))
ax.plot(pd.DatetimeIndex(dates), 100 * (np.exp(actual) - 1), color=style.SERIES[0], lw=1.2, label="realised next-10y real return")
ax.plot(pd.DatetimeIndex(dates), 100 * (np.exp(preds) - 1), color=style.SERIES[1], lw=1.2, label="CAPE forecast (real time)")
ax.plot(pd.DatetimeIndex(dates), 100 * (np.exp(bench) - 1), color=style.NEUTRAL, lw=1.2, ls="--", label="historical mean (real time)")
ax.axhline(0, color=style.INK_2, lw=0.8)
ax.set_ylabel("%/yr")
ax.set_title("Real-time forecasts of the next decade's real return")
ax.legend(fontsize=8, ncol=3, loc="lower left")
style.save(fig, "25_oos.png")

open(style.os.path.join(style.FIG_DIR, "25_output.txt"), "w").write("\n".join(out) + "\n")
