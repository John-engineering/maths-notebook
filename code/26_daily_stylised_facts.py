"""Note 26: is real volatility rough? Do GARCH fits imply the right tails?

Data: daily S&P 500 and NASDAQ Composite OHLC 1999-2018 (arch package), daily
VIX 1990-2026. (The S&P 'Open' field is just the previous close for much of
1999-2005, so nothing here uses opens.)

1. Roughness of log volatility. Daily variance from the Parkinson range
   estimator (ln H/L)^2 / (4 ln 2). Estimator noise adds a constant to the
   second-moment structure function, so fit m2(D) = c D^{2H} + 2 eta^2.
2. Stylised facts: tails (Hill), slow decay of |r| autocorrelation, leverage.
3. GARCH(1,1) and FIGARCH fits (arch): Kesten-implied tail exponent (note 14)
   vs Hill, and GARCH persistence vs FIGARCH (note 15).
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from arch import arch_model
import data
import kesten
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


spx = data.sp500_daily()
ndx = data.nasdaq_daily()
vix = data.vix_daily()


def parkinson_logvol(df):
    v = np.log(df["High"] / df["Low"]) ** 2 / (4 * np.log(2))
    v = v[v > 0]
    return 0.5 * np.log(v)


def structure(lv, lags, q=2):
    x = lv.values
    return np.array([np.mean(np.abs(x[k:] - x[:-k]) ** q) for k in lags])


lags = np.unique(np.geomspace(1, 250, 40).astype(int))


def fit_H(m2, lags, max_lag=100):
    sel = lags <= max_lag
    f = lambda L, c, H, n2: c * L ** (2 * H) + n2
    p, cov = curve_fit(f, lags[sel], m2[sel], p0=[0.05, 0.15, 0.1], bounds=([0, 0.001, 0], [10, 1, 10]))
    return p, np.sqrt(np.diag(cov))


series = {
    "S&P 500 range vol": parkinson_logvol(spx),
    "NASDAQ range vol": parkinson_logvol(ndx),
    "VIX (implied vol)": np.log(vix["close"]),
}
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
for i, (name, lv) in enumerate(series.items()):
    m2 = structure(lv, lags)
    (c, H, n2), se = fit_H(m2, lags)
    naive = np.polyfit(np.log(lags[(lags >= 1) & (lags <= 100)]),
                       np.log(m2[(lags >= 1) & (lags <= 100)]), 1)[0] / 2
    log(f"{name:20s} {lv.index[0]:%Y}-{lv.index[-1]:%Y}: noise-corrected H = {H:.3f} ± {se[1]:.3f} "
        f"(noise var 2eta^2 = {n2:.3f}); naive log-log slope/2 = {naive:.3f}")
    axes[0].loglog(lags, m2, marker="o", ms=3, lw=0, color=style.SERIES[i], label=f"{name}: H ≈ {H:.2f}")
    LL = np.geomspace(1, 250, 100)
    axes[0].loglog(LL, c * LL ** (2 * H) + n2, color=style.SERIES[i], lw=1.2)
    axes[1].loglog(lags, np.maximum(m2 - n2, 1e-6), marker="o", ms=3, color=style.SERIES[i],
                   label=f"{name}, noise removed")
LL = np.geomspace(1, 250, 10)
axes[1].loglog(LL, 0.02 * LL ** 0.2, color=style.NEUTRAL, ls="--", lw=1, label="slope 2H = 0.2 (H = 0.1)")
axes[1].loglog(LL, 0.004 * LL ** 1.0, color=style.NEUTRAL, ls=":", lw=1, label="slope 2H = 1 (H = ½)")
axes[0].set_xlabel("lag (trading days)")
axes[0].set_ylabel("E[(Δ log σ)²]")
axes[0].set_title("Structure function of log volatility")
axes[0].legend(fontsize=8)
axes[1].set_xlabel("lag (trading days)")
axes[1].set_ylabel("E[(Δ log σ)²] − noise")
axes[1].set_title("After removing estimator noise")
axes[1].legend(fontsize=7.5)
for ax in axes:
    style.plain_log(ax, "x")
style.save(fig, "26_roughness.png")

for name, df_ in [("S&P 500", spx), ("NASDAQ", ndx)]:
    lv = parkinson_logvol(df_)
    for a_, b_ in [("1999", "2008"), ("2009", "2018")]:
        sub = lv[a_:b_]
        (c, H, n2), se = fit_H(structure(sub, lags), lags)
        log(f"   {name} {a_}-{b_}: noise-corrected H = {H:.3f} ± {se[1]:.3f}")

# ------------------------------------------------------------ stylised facts
r = np.log(spx["Close"]).diff().dropna()
rn = np.log(ndx["Close"]).diff().dropna()


def hill(x, frac):
    xs = np.sort(np.abs(x))[::-1]
    k = int(frac * len(xs))
    return 1 / np.mean(np.log(xs[:k] / xs[k]))


for name, x in [("S&P 500", r), ("NASDAQ", rn)]:
    z = x.values - x.mean()
    log(f"{name} daily returns: kurtosis {pd.Series(z).kurt()+3:.1f}; Hill at top 1/2/5%: "
        f"{hill(z,0.01):.2f} / {hill(z,0.02):.2f} / {hill(z,0.05):.2f}")
# |r| autocorrelation decay
lag_a = np.unique(np.geomspace(1, 500, 30).astype(int))
a = np.abs(r.values) - np.abs(r.values).mean()
acf_abs = np.array([np.dot(a[k:], a[:-k]) / np.dot(a, a) for k in lag_a])
sel = (lag_a >= 5) & (lag_a <= 250)
decay = -np.polyfit(np.log(lag_a[sel]), np.log(acf_abs[sel]), 1)[0]
log(f"S&P |r| autocorrelation: lag 1 {acf_abs[0]:.2f}, lag 20 {acf_abs[np.argmin(np.abs(lag_a-20))]:.2f}, "
    f"lag 250 {acf_abs[np.argmin(np.abs(lag_a-250))]:.2f}; power-law decay exponent (5-250d) {decay:.2f}")
# leverage: corr(r_t, |r|_{t+k}) vs corr(|r|_t, r_{t+k})
lev_f = [np.corrcoef(r.values[:-k], np.abs(r.values[k:]))[0, 1] for k in range(1, 21)]
lev_b = [np.corrcoef(np.abs(r.values[:-k]), r.values[k:])[0, 1] for k in range(1, 21)]
log(f"leverage: corr(r_t, |r|_(t+k)) for k=1..5: " + ", ".join(f"{v:+.3f}" for v in lev_f[:5])
    + f"; reverse corr(|r|_t, r_(t+k)): " + ", ".join(f"{v:+.3f}" for v in lev_b[:5]))

# ------------------------------------------------------------ GARCH / FIGARCH
res = {}
for name, x in [("S&P 500", r), ("NASDAQ", rn)]:
    y = 100 * x
    g = arch_model(y, vol="GARCH", p=1, q=1, dist="normal").fit(disp="off")
    gt = arch_model(y, vol="GARCH", p=1, q=1, dist="t").fit(disp="off")
    fg = arch_model(y, vol="FIGARCH", p=1, q=1, dist="normal").fit(disp="off")
    a_, b_ = g.params["alpha[1]"], g.params["beta[1]"]
    at, bt, nu = gt.params["alpha[1]"], gt.params["beta[1]"], gt.params["nu"]
    k_n = 2 * kesten.kappa(a_, b_)
    k_t = 2 * kesten.kappa(at, bt, nu=nu)
    log(f"{name} GARCH(1,1) normal: alpha {a_:.3f} beta {b_:.3f} persistence {a_+b_:.4f} -> implied tail 2kappa {k_n:.2f}")
    log(f"{name} GARCH(1,1) t: alpha {at:.3f} beta {bt:.3f} persistence {at+bt:.4f}, nu {nu:.1f} -> implied "
        f"tail min(2kappa, nu) = {min(k_t, nu):.2f}")
    log(f"{name} FIGARCH: d = {fg.params['d']:.3f} (long-memory parameter); loglik FIGARCH {fg.loglikelihood:.1f} "
        f"vs GARCH {g.loglikelihood:.1f}")

open(style.os.path.join(style.FIG_DIR, "26_output.txt"), "w").write("\n".join(out) + "\n")
