"""Note 35: what follows stress? Drawdowns and credit spreads, 1926-2018.

1. 'Buy the dip': forward 12/36/60-month market excess returns after month-ends
   at which the market is >= 20% (or 30%) below its prior peak, vs all months.
   Null: stationary block bootstrap of monthly returns (blocks of 12), same
   conditioning rule applied to each bootstrap path.
2. Credit spread (Moody's BAA - AAA): does it predict next-12-month excess
   returns (with a persistent-regressor null: AR(1) spread independent of
   returns) and next-12-month volatility?
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from arch.data import default as moody
import data
import style

style.apply()
rng = np.random.default_rng(35)
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


ff = data.french_monthly()
ex = ff["Mkt-RF"]
tot = ex + ff["RF"]
lx = np.log1p(ex)


def fwd_sum(x, h):
    return x[::-1].rolling(h).sum()[::-1].shift(-1)


def drawdown(tot_):
    w = (1 + tot_).cumprod()
    return w / w.cummax() - 1


dd = drawdown(tot)
res = {}
for thr in (0.2, 0.3):
    cond = dd <= -thr
    for h in (12, 36, 60):
        f = fwd_sum(lx, h) * 12 / h
        a = f[cond].dropna().mean()
        b = f.dropna().mean()
        res[(thr, h)] = (a, b, cond.sum())
        log(f"drawdown >= {int(100*thr)}% ({cond.mean()*100:.0f}% of months): next {h}m excess return "
            f"{100*a:+.1f}%/yr vs unconditional {100*b:+.1f}%/yr")

# block bootstrap null for the difference (thr 20%, h = 36)
def block_boot(x, block=12):
    n = len(x)
    idx = []
    while len(idx) < n:
        s = rng.integers(0, n - block)
        idx.extend(range(s, s + block))
    return x[np.array(idx[:n])]


diffs = {h: [] for h in (12, 36, 60)}
for _ in range(2000):
    e = pd.Series(block_boot(ex.values), index=ex.index)
    t_ = e + ff["RF"].values
    d_ = drawdown(t_)
    l_ = np.log1p(e)
    c_ = d_ <= -0.2
    for h in (12, 36, 60):
        f_ = fwd_sum(l_, h) * 12 / h
        if c_.sum() > 5:
            diffs[h].append(f_[c_].dropna().mean() - f_.dropna().mean())
for h in (12, 36, 60):
    a, b, n = res[(0.2, h)]
    dv = np.array(diffs[h])
    log(f"   null (block bootstrap) for 20% drawdown, {h}m: 95% band of difference "
        f"[{100*np.percentile(dv,2.5):+.1f}, {100*np.percentile(dv,97.5):+.1f}]%/yr; observed {100*(a-b):+.1f}; "
        f"p(one-sided) = {np.mean(dv >= a - b):.3f}")
# episodes: distinct drawdown spells
spells = (dd <= -0.2).astype(int).diff().fillna(0)
log(f"distinct 20%+ drawdown episodes since 1926: {int((spells == 1).sum())}")

# ------------------------------------------------------------ credit spreads
m = moody.load()
spread = (m["BAA"] - m["AAA"]).rename("spread")
spread.index = spread.index.to_period("M").to_timestamp()
d = pd.DataFrame({"spread": spread, "ex": ex, "lx": lx}).dropna()
d["fwd12"] = fwd_sum(d["lx"], 12)
d["vol12"] = np.sqrt((d["ex"] ** 2)[::-1].rolling(12).mean()[::-1].shift(-1) * 12)
z = d.dropna()
X = np.column_stack([np.ones(len(z)), z["spread"]])
b_ret, *_ = np.linalg.lstsq(X, z["fwd12"].values, rcond=None)
r2_ret = 1 - np.var(z["fwd12"] - X @ b_ret) / np.var(z["fwd12"])
b_vol, *_ = np.linalg.lstsq(X, np.log(z["vol12"]).values, rcond=None)
r2_vol = 1 - np.var(np.log(z["vol12"]) - X @ b_vol) / np.var(np.log(z["vol12"]))
log(f"credit spread {z.index[0]:%Y}-{z.index[-1]:%Y}: mean {z['spread'].mean():.2f}pp, AR(1) "
    f"{np.corrcoef(z['spread'][1:], z['spread'][:-1])[0,1]:.3f}")
log(f"next-12m excess return on spread: slope {100*b_ret[1]:+.1f}%/yr per pp, R^2 {r2_ret:.3f}")
log(f"next-12m log volatility on spread: slope {b_vol[1]:+.3f} per pp, R^2 {r2_vol:.3f}")
# persistent-regressor null: AR(1) spreads independent of returns
phi = np.corrcoef(z["spread"][1:], z["spread"][:-1])[0, 1]
sd = z["spread"].std() * np.sqrt(1 - phi**2)
r2_null = []
for _ in range(2000):
    s = np.empty(len(z))
    s[0] = z["spread"].mean()
    e = rng.standard_normal(len(z)) * sd
    for i in range(1, len(z)):
        s[i] = z["spread"].mean() * (1 - phi) + phi * s[i - 1] + e[i]
    Xs = np.column_stack([np.ones(len(z)), s])
    bs, *_ = np.linalg.lstsq(Xs, z["fwd12"].values, rcond=None)
    r2_null.append(1 - np.var(z["fwd12"] - Xs @ bs) / np.var(z["fwd12"]))
r2_null = np.array(r2_null)
log(f"   return R^2 null (independent AR(1) spread): median {np.median(r2_null):.3f}, 95th pct "
    f"{np.percentile(r2_null,95):.3f}; p = {np.mean(r2_null >= r2_ret):.3f}")
r2v_null = []
for _ in range(1000):
    s = np.empty(len(z))
    s[0] = z["spread"].mean()
    e = rng.standard_normal(len(z)) * sd
    for i in range(1, len(z)):
        s[i] = z["spread"].mean() * (1 - phi) + phi * s[i - 1] + e[i]
    Xs = np.column_stack([np.ones(len(z)), s])
    y = np.log(z["vol12"]).values
    bs, *_ = np.linalg.lstsq(Xs, y, rcond=None)
    r2v_null.append(1 - np.var(y - Xs @ bs) / np.var(y))
r2v_null = np.array(r2v_null)
log(f"   volatility R^2 null: median {np.median(r2v_null):.3f}, 95th pct {np.percentile(r2v_null,95):.3f}; "
    f"p = {np.mean(r2v_null >= r2_vol):.3f}")

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2))
fig.subplots_adjust(wspace=0.25)
hs = [12, 36, 60]
x = np.arange(3)
axes[0].bar(x - 0.2, [100 * res[(0.2, h)][1] for h in hs], width=0.2, color=style.NEUTRAL, label="all months")
axes[0].bar(x, [100 * res[(0.2, h)][0] for h in hs], width=0.2, color=style.SERIES[0], label="after ≥ 20% drawdown")
axes[0].bar(x + 0.2, [100 * res[(0.3, h)][0] for h in hs], width=0.2, color=style.SERIES[1], label="after ≥ 30% drawdown")
axes[0].set_xticks(x)
axes[0].set_xticklabels([f"next {h} months" for h in hs])
axes[0].set_ylabel("annualised excess log return (%)")
axes[0].set_title("Buying the dip, 1926–2018")
axes[0].legend(fontsize=8)
axes[1].scatter(z["spread"], 100 * z["vol12"], s=5, alpha=0.4, color=style.SERIES[0], edgecolors="none")
ss = np.linspace(z["spread"].min(), z["spread"].max(), 50)
axes[1].plot(ss, 100 * np.exp(b_vol[0] + b_vol[1] * ss), color=style.SERIES[1], lw=1.5)
axes[1].set_yscale("log")
axes[1].set_xlabel("BAA − AAA credit spread (pp)")
axes[1].set_ylabel("next-12-month market volatility (%, log)")
axes[1].set_title(f"Credit stress predicts volatility (R² = {r2_vol:.2f})")
axes[1].set_yticks([5, 10, 20, 40, 80])
style.plain_log(axes[1])
style.save(fig, "35_after_stress.png")

open(style.os.path.join(style.FIG_DIR, "35_output.txt"), "w").write("\n".join(out) + "\n")
