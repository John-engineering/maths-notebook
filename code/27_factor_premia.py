"""Note 27: Fama-French factor premia 1926-2018 through the notebook's own tools.

1. Sharpe ratios and t-statistics by era (pre-1963 / 1963-1992 in-sample of
   Fama-French 1993 / post-publication 1993-2018).
2. Always-valid evidence (notes 07, 19): normal-mixture Kelly wealth for each
   factor, started in 1926 and again in 1993; date it first reaches 20x.
3. Quickest detection of death (note 18): CUSUM tuned on the pre-1993 Sharpe,
   run from 1993 with one false alarm per 20 years; date it fires.
4. Kelly shrinkage t^2/(1+t^2) (note 05) at the end of each era.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import data
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


ff = data.french_monthly()
factors = ["Mkt-RF", "SMB", "HML"]
eras = [("1926-07", "1962-12"), ("1963-01", "1992-12"), ("1993-01", "2018-11")]
log(f"data {ff.index[0]:%Y-%m} to {ff.index[-1]:%Y-%m}")
log("factor   era          mean %/yr  vol %/yr  Sharpe   t-stat  Kelly shrink t²/(1+t²)")
for f in factors:
    for a, b in eras + [("1926-07", "2018-11")]:
        x = ff.loc[a:b, f]
        m, s = 12 * x.mean(), np.sqrt(12) * x.std()
        sr = m / s
        t = sr * np.sqrt(len(x) / 12)
        log(f"{f:7s}  {a[:4]}-{b[:4]}   {100*m:8.2f}  {100*s:8.1f}  {sr:6.2f}  {t:6.2f}   {t*t/(1+t*t):.2f}")


# ------------------------------------------------------------ always-valid mixture wealth
def mixture_wealth(full, start, rho=12.0):
    """Normal-mixture Kelly wealth (note 19) on monthly returns from `start`,
    each standardised by a lagged EWMA volatility (half-life 12 months, seeded
    with the first 24 months; no look-ahead), as a vol-aware Kelly bettor
    would. Returns index, wealth, cumulative sum."""
    ew = (full**2).ewm(halflife=12, min_periods=24).mean()
    vol = np.sqrt(ew).shift(1)
    z = (full / vol).loc[start:].dropna()
    S = np.cumsum(z.values)
    t = np.arange(1, len(z) + 1)
    M = np.sqrt(rho / (t + rho)) * np.exp(S**2 / (2 * (t + rho)))
    return z.index, M, S


fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
for i, f in enumerate(factors):
    for start, ls in [("1928-07", "-"), ("1993-01", "--")]:
        idx, M, S = mixture_wealth(ff[f], start)
        hit = np.where((M >= 20) & (S > 0))[0]
        when = f"{idx[hit[0]]:%Y-%m}" if len(hit) else "never"
        log(f"mixture-Kelly wealth for {f} monitored from {start[:4]}: first reaches 20x (positive premium) "
            f"at {when}; final log-wealth {np.log(M[-1]):.2f} (sign of cumulative evidence {'+' if S[-1] > 0 else '-'})")
        axes[0].plot(idx, np.log(M), color=style.SERIES[i], ls=ls,
                     lw=1.4, label=f"{f} from {start[:4]}")
axes[0].axhline(np.log(20), color=style.INK_2, lw=0.8, ls=":")
axes[0].set_ylabel("log mixture-Kelly wealth (evidence)")
axes[0].set_title("Always-valid evidence: reject 'no premium' at 20×")
axes[0].legend(fontsize=7.5, ncol=2, loc="upper left")


# ------------------------------------------------------------ CUSUM for death, monitored from 1993
def cusum_alarm(x, sr_alive, h):
    """x: monthly returns standardised by the pre-1993 vol; alive mean d per month."""
    d = sr_alive / np.sqrt(12)
    W, path = 0.0, []
    for v in x:
        W = max(0.0, W - d * (v - d / 2))
        path.append(W)
    path = np.array(path)
    hit = np.where(path >= h)[0]
    return path, (hit[0] if len(hit) else None)


def h_for_arl(sr, arl_years):
    """Continuous-time CUSUM (note 18): ARL = tau (e^h - h - 1), tau = 2/SR^2 years."""
    tau = 2 / sr**2
    hs = np.linspace(0.01, 15, 20000)
    return hs[np.argmin(np.abs(tau * (np.exp(hs) - hs - 1) - arl_years))]


for i, f in enumerate(factors):
    pre = ff.loc[:"1992-12", f]
    sr = np.sqrt(12) * pre.mean() / pre.std()
    post = ff.loc["1993-01":, f] / pre.std()
    h = h_for_arl(sr, 20)
    path, k = cusum_alarm(post.values, sr, h)
    when = f"{post.index[k]:%Y-%m}" if k is not None else "no alarm by 2018-11"
    # probability of at least one false alarm in this monitoring window if the premium were alive
    rr = np.random.default_rng(i)
    d = sr / np.sqrt(12)
    sims = d + rr.standard_normal((4000, len(post)))
    Wm = np.zeros(4000)
    fired = np.zeros(4000, bool)
    for kk in range(len(post)):
        Wm = np.maximum(0.0, Wm - d * (sims[:, kk] - d / 2))
        fired |= Wm >= h
    after = ""
    if k is not None:
        xa = ff.loc[post.index[k]:, f].iloc[1:]
        after = f"; Sharpe AFTER the alarm to 2018: {np.sqrt(12)*xa.mean()/xa.std():+.2f} over {len(xa)/12:.0f}y"
    log(f"CUSUM 'premium died' for {f}: alive Sharpe (1926-92) {sr:.2f}, h={h:.2f} "
        f"(1 false alarm / 20y): alarm {when}; P(alarm by 2018 | still alive) = {fired.mean():.2f}{after}")
    axes[1].plot(post.index, path / h, color=style.SERIES[i], lw=1.4, label=f"{f} (pre-1993 Sharpe {sr:.2f})")
axes[1].axhline(1.0, color=style.INK_2, lw=0.8, ls=":")
axes[1].set_ylabel("CUSUM / alarm threshold")
axes[1].set_title("From 1993: has the premium died? (alarm at 1)")
axes[1].legend(fontsize=8, loc="upper left")
style.save(fig, "27_factor_evidence.png")

# rolling 10-year Sharpe
fig, ax = plt.subplots(figsize=(10, 3.8))
for i, f in enumerate(factors):
    rs = ff[f].rolling(120).mean() / ff[f].rolling(120).std() * np.sqrt(12)
    ax.plot(rs.index, rs, color=style.SERIES[i], lw=1.3, label=f)
ax.axhline(0, color=style.INK_2, lw=0.8)
ax.axvline(pd.Timestamp("1993-01-01"), color=style.NEUTRAL, ls=":", lw=1)
ax.annotate("Fama & French (1993)", (pd.Timestamp("1993-06-01"), 1.1), fontsize=8, color=style.INK_2)
ax.set_ylabel("trailing 10-year Sharpe")
ax.set_title("Factor premia over time")
ax.legend(fontsize=8)
style.save(fig, "27_rolling_sharpe.png")

open(style.os.path.join(style.FIG_DIR, "27_output.txt"), "w").write("\n".join(out) + "\n")
