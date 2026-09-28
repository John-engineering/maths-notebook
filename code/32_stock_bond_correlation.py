"""Note 32: the stock-bond correlation, 1871-2023. Does it have regimes, and
what drives them?

Bond returns are built from Shiller's long-rate series (10-year Treasury yield
from 1953; long government yields before): each month, hold a par bond with
coupon equal to last month's yield and ~10 years to maturity, reprice it at
this month's yield. (Both series are monthly averages, so both are smoothed in
the same way.)
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


df = data.sp500_monthly()
df = df[df.index <= "2023-09-01"]
y = df["rate"] / 100


def par_bond_return(y_prev, y_now, years=10.0):
    """Semiannual par bond (coupon = y_prev), repriced at y_now after one month."""
    n = int(round(2 * years))
    c = y_prev / 2
    t = np.arange(1, n + 1) / 2 - 1 / 12
    t = t[t > 0]
    price = np.sum(c * (1 + y_now / 2) ** (-2 * t)) + (1 + y_now / 2) ** (-2 * t[-1])
    return price - 1          # dirty price: the month's coupon accrual is included


bond = pd.Series([np.nan] + [par_bond_return(a, b) for a, b in zip(y.values[:-1], y.values[1:])], index=y.index)
stock = df["ret_total"]
infl_yoy = df["cpi"].pct_change(12)
d = pd.DataFrame({"s": stock, "b": bond, "infl": infl_yoy, "y": y}).dropna()
log(f"sample {d.index[0]:%Y-%m} to {d.index[-1]:%Y-%m}; bond total return {100*12*d['b'].mean():.2f}%/yr "
    f"(vol {100*np.sqrt(12)*d['b'].std():.1f}%)")

roll = d["s"].rolling(60).corr(d["b"])
infl_roll = d["infl"].rolling(60).mean()
infl_vol = d["infl"].rolling(60).std()
for a, b in [("1871", "1913"), ("1914", "1945"), ("1946", "1965"), ("1966", "1999"), ("2000", "2021"), ("2022", "2023")]:
    x = d.loc[a:b]
    c = np.corrcoef(x["s"], x["b"])[0, 1] if len(x) > 12 else np.nan
    log(f"   {a}-{b}: stock-bond correlation {c:+.2f}, mean inflation {100*x['infl'].mean():.1f}%, "
        f"inflation sd {100*x['infl'].std():.1f}%")

# what explains the rolling correlation? regress on trailing inflation level and variability
z = pd.DataFrame({"corr": roll, "infl": infl_roll, "ivol": infl_vol, "y": d["y"]}).dropna()
zz = z.iloc[::60]                                        # non-overlapping 5-year windows
X = np.column_stack([np.ones(len(zz)), zz["infl"], zz["ivol"]])
beta, *_ = np.linalg.lstsq(X, zz["corr"].values, rcond=None)
res = zz["corr"].values - X @ beta
r2 = 1 - res.var() / zz["corr"].var()
se = np.sqrt(np.diag(np.linalg.inv(X.T @ X) * res.var() * len(zz) / (len(zz) - 3)))
log(f"non-overlapping 5-year windows ({len(zz)}): corr = {beta[0]:+.2f} {beta[1]:+.2f}*inflation "
    f"{beta[2]:+.2f}*inflation-vol; t = {beta[1]/se[1]:+.2f}, {beta[2]/se[2]:+.2f}; R^2 {r2:.2f}")
# sign agreement: high inflation windows
hi = z["infl"] > 0.04
log(f"rolling 5y correlation when trailing inflation > 4%: mean {z.loc[hi,'corr'].mean():+.2f}; "
    f"when < 4%: {z.loc[~hi,'corr'].mean():+.2f}")

# consequence for a 60/40 portfolio: volatility with each era's correlation
for a, b in [("1966", "1999"), ("2000", "2021")]:
    x = d.loc[a:b]
    p = 0.6 * x["s"] + 0.4 * x["b"]
    ind = np.sqrt(0.36 * x["s"].var() + 0.16 * x["b"].var())
    log(f"   60/40 {a}-{b}: vol {100*np.sqrt(12)*p.std():.1f}% vs {100*np.sqrt(12)*ind:.1f}% if uncorrelated")

fig, axes = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True)
fig.subplots_adjust(hspace=0.3)
axes[0].plot(roll.index, roll, color=style.SERIES[0], lw=1.3)
axes[0].axhline(0, color=style.INK_2, lw=0.8)
axes[0].set_ylabel("5-year rolling correlation")
axes[0].set_title("Stock–bond correlation, US, 1876–2023")
axes[1].plot(infl_roll.index, 100 * infl_roll, color=style.SERIES[1], lw=1.3)
axes[1].axhline(4, color=style.NEUTRAL, lw=0.8, ls=":")
axes[1].set_ylabel("5-year average inflation (%)")
axes[1].set_title("Trailing inflation")
style.save(fig, "32_stock_bond.png")

open(style.os.path.join(style.FIG_DIR, "32_output.txt"), "w").write("\n".join(out) + "\n")
