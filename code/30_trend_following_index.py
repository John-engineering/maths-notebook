"""Note 30: does trend following work on the US stock market, 1926-2018?

Month-end Fama-French market excess returns (Mkt-RF) and T-bill (RF).
Rules (signal known at the end of month t, applied to month t+1):
  TSMOM-12:  long if trailing 12-month excess return > 0, else cash (or short)
  SMA-10:    long if the total-return index is above its 10-month average
Nulls:
  (a) 'random timing': a long/cash rule that is long the same fraction of
      months but at random times (does timing add anything beyond exposure?);
  (b) iid bootstrap of returns (kills autocorrelation, keeps the drift).
Plus note 04's prediction: long/short trend P&L is convex in the market's move.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import data
import style

style.apply()
rng = np.random.default_rng(30)
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


ff = data.french_monthly()
ex = ff["Mkt-RF"]
rf = ff["RF"]
tot = ex + rf
idx_level = (1 + tot).cumprod()


def signals(ex_, tot_):
    lvl = (1 + tot_).cumprod()
    ts = (np.log1p(ex_).rolling(12).sum() > 0).astype(float)
    sma = (lvl > lvl.rolling(10).mean()).astype(float)
    return ts.shift(1), sma.shift(1)


def stats(x):
    x = x.dropna()
    sr = x.mean() / x.std() * np.sqrt(12)
    w = (1 + x).cumprod()
    mdd = (w / w.cummax() - 1).min()
    return sr, 12 * x.mean(), np.sqrt(12) * x.std(), mdd, x.skew()


ts, sma = signals(ex, tot)
strategies = {
    "buy and hold": ex,
    "TSMOM-12 long/cash": ts * ex,
    "SMA-10 long/cash": sma * ex,
    "TSMOM-12 long/short": (2 * ts - 1) * ex,
}
log(f"sample {ex.index[12]:%Y-%m} to {ex.index[-1]:%Y-%m}; excess returns over T-bills")
log("strategy                Sharpe   mean %/yr   vol %/yr   max DD   skew   months long")
for name, x in strategies.items():
    x = x.iloc[13:]
    sr, m, v, mdd, sk = stats(x)
    frac = {"buy and hold": 1.0, "TSMOM-12 long/cash": ts.iloc[13:].mean(), "SMA-10 long/cash": sma.iloc[13:].mean(),
            "TSMOM-12 long/short": np.nan}[name]
    log(f"{name:22s} {sr:6.2f}   {100*m:8.2f}   {100*v:8.1f}   {100*mdd:6.0f}%  {sk:+5.2f}   {frac:.2f}")

# eras
for a, b in [("1927", "1962"), ("1963", "1992"), ("1993", "2018")]:
    line = []
    for name in ["buy and hold", "TSMOM-12 long/cash", "SMA-10 long/cash"]:
        sr = stats(strategies[name].loc[a:b])[0]
        line.append(f"{name} {sr:.2f}")
    log(f"   {a}-{b}: " + ", ".join(line))

# ------------------------------------------------------------ null (a): random timing with same exposure
x_bh = ex.iloc[13:].values
for name, sig in [("TSMOM-12 long/cash", ts), ("SMA-10 long/cash", sma)]:
    s = sig.iloc[13:].values
    sr_real = stats(pd.Series(s * x_bh))[0]
    frac = s.mean()
    sims = []
    for _ in range(5000):
        # random regime timing with the same number of months long and similar
        # persistence: circularly shift the real signal (keeps its run lengths)
        k = rng.integers(12, len(s) - 12)
        sr_ = pd.Series(np.roll(s, k) * x_bh)
        sims.append(sr_.mean() / sr_.std() * np.sqrt(12))
    sims = np.array(sims)
    log(f"{name}: Sharpe {sr_real:.2f}; same signal shifted to random dates: median {np.median(sims):.2f}, "
        f"95th pct {np.percentile(sims, 95):.2f}; p = {np.mean(sims >= sr_real):.3f}")

# null (b): iid bootstrap of returns (no autocorrelation at all)
sims_b = []
for _ in range(2000):
    xb = pd.Series(rng.choice(ex.values, len(ex)), index=ex.index)
    tb = xb + rf.values
    tsb, _ = signals(xb, tb)
    y = (tsb * xb).iloc[13:]
    sims_b.append(y.mean() / y.std() * np.sqrt(12) - xb.iloc[13:].mean() / xb.iloc[13:].std() * np.sqrt(12))
sims_b = np.array(sims_b)
d_real = stats(strategies["TSMOM-12 long/cash"].iloc[13:])[0] - stats(ex.iloc[13:])[0]
log(f"Sharpe improvement of TSMOM-12 over buy-and-hold: {d_real:+.2f}; iid-bootstrap null: median "
    f"{np.median(sims_b):+.2f}, 95th pct {np.percentile(sims_b, 95):+.2f}; p = {np.mean(sims_b >= d_real):.3f}")

# ------------------------------------------------------------ where does it come from?
x_ = ex.iloc[13:]
s_ = ts.iloc[13:]
vol_tr = np.sqrt((ex**2).rolling(12).mean() * 12).shift(1).iloc[13:]
for lab, m in [("signal long", s_ == 1), ("signal out", s_ == 0)]:
    log(f"{lab:12s}: {m.mean()*100:4.0f}% of months, mean excess {1200*x_[m].mean():+6.2f}%/yr, "
        f"vol {100*np.sqrt(12)*x_[m].std():5.1f}%, Sharpe {x_[m].mean()/x_[m].std()*np.sqrt(12):+.2f}, "
        f"trailing vol {100*vol_tr[m].mean():.1f}%")
# a pure volatility-timing rule with the same exposure: out when trailing vol is in its top 30%
thr = vol_tr.expanding(60).quantile(0.70)
vt = (vol_tr < thr).astype(float)
vt[thr.isna()] = 1.0
xv = vt * x_
log(f"pure vol-timing rule (out when trailing vol in top 30%): long {vt.mean():.2f} of months, "
    f"Sharpe {xv.mean()/xv.std()*np.sqrt(12):.2f}")
both = s_ * vt
xb_ = both * x_
log(f"TSMOM and vol-timing agree on 'out' in {100*((s_==0)&(vt==0)).sum()/max(1,(s_==0).sum()):.0f}% of TSMOM's out-months; "
    f"corr of the two signals {np.corrcoef(s_, vt)[0,1]:.2f}")
# predictability test: regress next-month excess return on the signal, controlling for trailing vol
X = np.column_stack([np.ones(len(x_)), s_.values, vol_tr.values])
beta, *_ = np.linalg.lstsq(X, x_.values, rcond=None)
resid = x_.values - X @ beta
cov = np.linalg.inv(X.T @ X) * resid.var()
log(f"regression next-month excess return = a + b*signal + c*trailing vol: b = {1200*beta[1]:+.2f}%/yr "
    f"(t {beta[1]/np.sqrt(cov[1,1]):+.2f}), c = {beta[2]:+.3f} (t {beta[2]/np.sqrt(cov[2,2]):+.2f})")

# ------------------------------------------------------------ convexity (note 04)
ann_mkt = np.log1p(ex).groupby(ex.index.year).sum()
ann_ls = np.log1p(strategies["TSMOM-12 long/short"].fillna(0)).groupby(ex.index.year).sum()
ann = pd.DataFrame({"mkt": ann_mkt, "ls": ann_ls}).iloc[2:]
quad = np.polyfit(ann["mkt"], ann["ls"], 2)
log(f"annual long/short TSMOM vs market: quadratic fit ls = {quad[0]:+.2f} m^2 {quad[1]:+.2f} m {quad[2]:+.3f}; "
    f"corr(ls, m^2) = {np.corrcoef(ann['ls'], ann['mkt']**2)[0,1]:.2f}")
for yr in [1931, 1937, 1974, 2002, 2008]:
    if yr in ann.index:
        log(f"   {yr}: market {100*(np.exp(ann.loc[yr,'mkt'])-1):+.0f}%, TSMOM long/short {100*(np.exp(ann.loc[yr,'ls'])-1):+.0f}%")

# ------------------------------------------------------------ figures
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
for i, name in enumerate(["buy and hold", "TSMOM-12 long/cash", "SMA-10 long/cash"]):
    w = (1 + strategies[name].iloc[13:].fillna(0)).cumprod()
    axes[0].plot(w.index, w, color=style.SERIES[i], lw=1.3, label=name)
axes[0].set_yscale("log")
axes[0].set_ylabel("growth of $1 (excess of T-bills, log)")
axes[0].set_title("US market, 1927–2018")
axes[0].legend(fontsize=8)
style.plain_log(axes[0])
axes[1].scatter(100 * (np.exp(ann["mkt"]) - 1), 100 * (np.exp(ann["ls"]) - 1), s=14, color=style.SERIES[0], alpha=0.8)
mm = np.linspace(ann["mkt"].min(), ann["mkt"].max(), 100)
axes[1].plot(100 * (np.exp(mm) - 1), 100 * (np.exp(np.polyval(quad, mm)) - 1), color=style.SERIES[1], lw=1.5)
axes[1].axhline(0, color=style.INK_2, lw=0.8)
axes[1].axvline(0, color=style.INK_2, lw=0.8)
axes[1].set_xlabel("market excess return in the year (%)")
axes[1].set_ylabel("long/short TSMOM return (%)")
axes[1].set_title("Trend following is a straddle on the market (note 04)")
style.save(fig, "30_trend_following.png")

open(style.os.path.join(style.FIG_DIR, "30_output.txt"), "w").write("\n".join(out) + "\n")
