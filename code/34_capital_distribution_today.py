"""Note 34: the shape of today's market. Zipf, concentration and valuation.

Data: a recent snapshot of S&P 500 constituents (market cap, P/E, sector),
github.com/datasets/s-and-p-500-companies-financials (prices imply ~2026).
1. Rank-size (capital distribution) curve: Pareto tail exponent by log-log
   regression and Hill; compare with Zipf (1).
2. Concentration: top-10 share, effective number of stocks 1/sum(w^2).
3. Does a zero-alpha compounding market (note 01's simulation) produce the same
   shape? Compare the simulated end-of-sample capital distribution.
4. Valuation dispersion: log P/E vs size.
"""
import os
import urllib.request
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


path = os.path.join(data.RAW, "sp500_financials.csv")
if not os.path.exists(path):
    urllib.request.urlretrieve("https://raw.githubusercontent.com/datasets/s-and-p-500-companies-financials/"
                               "main/data/constituents-financials.csv", path)
df = pd.read_csv(path)
cap = df["Market Cap"].dropna()
cap = cap[cap > 1e9].sort_values(ascending=False).reset_index(drop=True)   # drop obvious glitches
n = len(cap)
w = cap / cap.sum()
rank = np.arange(1, n + 1)
log(f"{n} constituents with market cap; total ${cap.sum()/1e12:.1f}T; largest ${cap.iloc[0]/1e12:.2f}T, "
    f"median ${cap.median()/1e9:.0f}B")
log(f"top-1 share {100*w.iloc[0]:.1f}%, top-10 {100*w.iloc[:10].sum():.1f}%, top-50 {100*w.iloc[:50].sum():.1f}%; "
    f"effective number of stocks 1/sum(w^2) = {1/np.sum(w**2):.0f}")


def tail_fit(c, k):
    lr, lc = np.log(np.arange(1, k + 1)), np.log(c[:k].values)
    slope = np.polyfit(lc, lr, 1)[0]          # log rank = a - zeta log cap
    hill = 1 / np.mean(np.log(c[:k].values / c.iloc[k]))
    return -slope, hill


for k in (20, 50, 100, 250):
    z_ols, z_hill = tail_fit(cap, k)
    log(f"Pareto tail exponent, top {k}: rank-size regression {z_ols:.2f}, Hill {z_hill:.2f} (Zipf = 1)")

# Gabaix-Ibragimov rank-1/2 regression for the top 100 (less biased)
k = 100
lr = np.log(np.arange(1, k + 1) - 0.5)
slope = np.polyfit(np.log(cap[:k].values), lr, 1)[0]
log(f"Gabaix-Ibragimov (rank - 1/2) estimate, top 100: {-slope:.2f} ± {abs(slope)*np.sqrt(2/k):.2f}")

# ------------------------------------------------------------ zero-alpha compounding comparison
import importlib.util
spec = importlib.util.spec_from_file_location("b01", os.path.join(style.HERE, "01_bessembinder.py"))
src = open(os.path.join(style.HERE, "01_bessembinder.py")).read().split("results = [simulate(s)")[0]
ns = {}
exec(compile(src, "01_bessembinder.py", "exec"), ns)


def end_caps(seed, n_firms=25_000):
    """Re-run note 01's zero-alpha market and return the caps of firms alive at the end."""
    r = np.random.default_rng(seed)
    MONTHS = ns["MONTHS"]
    mkt_ex = ns["ERP_ANNUAL"] / 12 + ns["SIG_M"] / np.sqrt(12) * r.standard_normal(MONTHS)
    birth = r.integers(0, MONTHS - 1, n_firms)
    life = np.ceil(r.exponential(10.8 * 12, n_firms)).astype(int)
    death = birth + life
    alive_end = death >= MONTHS
    beta = np.clip(r.normal(1.0, 0.35, n_firms), 0.2, 2.5)
    cap0 = np.exp(r.normal(0, 1.6, n_firms))
    sig_id = 0.45 * cap0 ** -0.12 * np.exp(0.25 * r.standard_normal(n_firms))
    rf = ns["rf"]
    capT = np.empty(n_firms)
    cum_m = np.concatenate([[0], np.cumsum(np.log1p(rf + mkt_ex))])
    for i in np.where(alive_end)[0]:
        T = MONTHS - birth[i]
        s = sig_id[i] / np.sqrt(12)
        idio = r.standard_normal(T) * s - s * s / 2
        mk = np.log1p(rf + beta[i] * mkt_ex[birth[i]:])
        capT[i] = cap0[i] * np.exp(mk.sum() + idio.sum())
    return np.sort(capT[alive_end])[::-1]


sim_curves = []
for seed in range(5):
    c = end_caps(seed)
    top = c[:n]
    ws = top / top.sum()
    z_s = -np.polyfit(np.log(top[:100]), np.log(np.arange(1, 101) - 0.5), 1)[0]
    sim_curves.append(ws)
    log(f"zero-alpha market (note 01), seed {seed}: largest-{n} firms: top-10 share {100*ws[:10].sum():.1f}%, "
        f"effective N {1/np.sum(ws**2):.0f}, tail exponent (top 100) {z_s:.2f}")

# ------------------------------------------------------------ valuation vs size
pe = df.loc[df["Market Cap"] > 1e9, ["Market Cap", "Price/Earnings", "Sector"]].dropna()
pe = pe[(pe["Price/Earnings"] > 0) & (pe["Price/Earnings"] < 500)]
b = np.polyfit(np.log(pe["Market Cap"]), np.log(pe["Price/Earnings"]), 1)
top10 = pe.nlargest(10, "Market Cap")
log(f"P/E: median {pe['Price/Earnings'].median():.1f}; cap-weighted harmonic P/E "
    f"{pe['Market Cap'].sum()/(pe['Market Cap']/pe['Price/Earnings']).sum():.1f}; "
    f"elasticity of P/E to size {b[0]:+.3f} (a 10x larger firm trades at {np.exp(b[0]*np.log(10)):.2f}x the P/E); "
    f"top-10 median P/E {top10['Price/Earnings'].median():.1f}")

# ------------------------------------------------------------ figure
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
axes[0].loglog(rank, w, color=style.SERIES[0], lw=2, label="S&P 500 today")
for j, ws in enumerate(sim_curves):
    axes[0].loglog(np.arange(1, len(ws) + 1), ws, color=style.NEUTRAL, lw=0.8,
                   label="zero-alpha compounding (note 01)" if j == 0 else None)
axes[0].loglog(rank, w.iloc[0] * rank ** -1.0, color=style.SERIES[1], ls="--", lw=1, label="Zipf: weight ∝ 1/rank")
axes[0].set_xlabel("rank")
axes[0].set_ylabel("index weight")
axes[0].set_title("Capital distribution curve")
axes[0].legend(fontsize=8)
style.plain_log(axes[0])
style.plain_log(axes[0], "x")
axes[1].scatter(pe["Market Cap"] / 1e9, pe["Price/Earnings"], s=8, alpha=0.5, color=style.SERIES[0], edgecolors="none")
xx = np.geomspace(pe["Market Cap"].min(), pe["Market Cap"].max(), 50)
axes[1].plot(xx / 1e9, np.exp(np.polyval(b, np.log(xx))), color=style.SERIES[1], lw=1.5)
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlabel("market cap ($bn, log)")
axes[1].set_ylabel("P/E (log)")
axes[1].set_title("Valuation vs size")
style.plain_log(axes[1])
style.plain_log(axes[1], "x")
style.save(fig, "34_capital_distribution.png")

open(style.os.path.join(style.FIG_DIR, "34_output.txt"), "w").write("\n".join(out) + "\n")
