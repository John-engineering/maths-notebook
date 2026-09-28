"""Note 29: calendar anomalies with multiple testing and publication dates.

1. Month-of-year effects, Shiller monthly S&P 1871-2026 (price returns; the
   averaging artefact of note 24 smears returns across adjacent months, so we
   also use Fama-French month-end market returns 1926-2018).
   12 simultaneous tests -> Holm correction; empirical-Bayes shrinkage of the
   12 monthly means (James-Stein toward the grand mean).
2. "Sell in May": Nov-Apr minus May-Oct, before and after Bouman & Jacobsen (2002).
3. Turn of the month (last trading day + first 3) and day of week, daily S&P
   1999-2018.
4. Overnight vs intraday, NASDAQ 1999-2018 (genuine opens; the S&P open field
   is unusable before ~2006).
"""
import numpy as np
import pandas as pd
from scipy.stats import t as tdist
import matplotlib.pyplot as plt
import data
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


def tstat(x):
    x = np.asarray(x)
    return x.mean() / (x.std(ddof=1) / np.sqrt(len(x)))


def holm(pvals):
    p = np.asarray(pvals)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(p) - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


# ------------------------------------------------------------ 1. month of year (FF month-end returns)
ff = data.french_monthly()
mk = ff["Mkt-RF"] + ff["RF"]
months = mk.index.month
rows = []
for m in range(1, 13):
    x = mk[months == m]
    rest = mk[months != m]
    diff = x.mean() - rest.mean()
    se = np.sqrt(x.var() / len(x) + rest.var() / len(rest))
    tt = diff / se
    rows.append((m, 100 * x.mean(), tt, 2 * tdist.sf(abs(tt), len(x) - 1)))
tab = pd.DataFrame(rows, columns=["month", "mean %", "t vs rest", "p"])
tab["Holm p"] = holm(tab["p"].values)
# James-Stein shrinkage of the 12 monthly means toward the grand mean
means = tab["mean %"].values / 100
se_m = np.array([mk[months == m].std() / np.sqrt((months == m).sum()) for m in range(1, 13)])
grand = mk.mean()
dev = means - grand
k = len(means)
shrink = max(0.0, 1 - (k - 3) * np.mean(se_m**2) / np.sum(dev**2))
tab["shrunk mean %"] = 100 * (grand + shrink * dev)
log(f"Fama-French total market return by calendar month, {mk.index[0]:%Y}-{mk.index[-1]:%Y}:")
for _, r_ in tab.iterrows():
    log(f"   month {int(r_['month']):2d}: mean {r_['mean %']:+5.2f}%  t vs rest {r_['t vs rest']:+5.2f}  "
        f"p {r_['p']:.3f}  Holm {r_['Holm p']:.3f}  James-Stein shrunk {r_['shrunk mean %']:+5.2f}%")
log(f"James-Stein shrinkage factor toward the grand mean: {shrink:.2f} (0 = all differences are noise)")

# ------------------------------------------------------------ 2. Sell in May
sh = data.sp500_monthly()["ret_price"].dropna()
def halfyear(series, a, b):
    s = series[(series.index.year >= a) & (series.index.year <= b)]
    win = s[s.index.month.isin([11, 12, 1, 2, 3, 4])]
    summ = s[s.index.month.isin([5, 6, 7, 8, 9, 10])]
    # annual differences: sum of Nov-Apr minus sum of May-Oct within each season year
    yr = np.where(s.index.month >= 11, s.index.year + 1, s.index.year)
    g = pd.DataFrame({"r": np.log1p(s.values), "win": s.index.month.isin([11, 12, 1, 2, 3, 4]), "y": yr})
    agg = g.groupby(["y", "win"])["r"].sum().unstack()
    d_ = (agg[True] - agg[False]).dropna()
    return 100 * d_.mean(), tstat(d_), len(d_)
for label, series, spans in [("Shiller S&P (monthly averages)", sh, [(1872, 2001), (2002, 2025)]),
                             ("Fama-French market (month-end)", mk, [(1927, 2001), (2002, 2018)])]:
    for a, b in spans:
        m_, tt, n = halfyear(series, a, b)
        log(f"Sell in May, {label} {a}-{b}: Nov-Apr minus May-Oct = {m_:+.2f}%/yr, t = {tt:+.2f} ({n} years)")

# ------------------------------------------------------------ 3. turn of month, day of week (daily S&P)
spx = data.sp500_daily()["Close"]
r = np.log(spx).diff().dropna()
idx = r.index
ym = idx.to_period("M")
pos = pd.Series(range(len(idx)), index=idx).groupby(ym).rank(method="first").astype(int)
nlast = pd.Series(1, index=idx).groupby(ym).transform("sum")
from_end = nlast - pos
tom = (pos <= 3) | (from_end == 0)
x_tom, x_rest = r[tom], r[~tom]
tt = (x_tom.mean() - x_rest.mean()) / np.sqrt(x_tom.var() / len(x_tom) + x_rest.var() / len(x_rest))
log(f"turn of month (last day + first 3), S&P 1999-2018: mean {1e4*x_tom.mean():+.1f} bp/day vs "
    f"{1e4*x_rest.mean():+.1f} bp/day otherwise, t = {tt:+.2f}; TOM days are {100*tom.mean():.0f}% of days and "
    f"account for {100*x_tom.sum()/r.sum():.0f}% of the total log return")
for half, (a, b) in [("1999-2008", ("1999", "2008")), ("2009-2018", ("2009", "2018"))]:
    rr = r[a:b]
    tm = tom[a:b]
    log(f"   {half}: TOM {1e4*rr[tm].mean():+.1f} bp/day vs rest {1e4*rr[~tm].mean():+.1f} bp/day")
dow = r.groupby(idx.dayofweek)
ps = []
for d_, x in dow:
    rest = r[idx.dayofweek != d_]
    t_ = (x.mean() - rest.mean()) / np.sqrt(x.var() / len(x) + rest.var() / len(rest))
    ps.append(2 * tdist.sf(abs(t_), len(x) - 1))
    log(f"   weekday {['Mon','Tue','Wed','Thu','Fri'][d_]}: {1e4*x.mean():+.1f} bp/day, t vs rest {t_:+.2f}")
log(f"   day-of-week Holm-adjusted p-values: " + ", ".join(f"{p:.2f}" for p in holm(ps)))

# ------------------------------------------------------------ 4. overnight vs intraday (NASDAQ)
nq = data.nasdaq_daily()
on = np.log(nq["Open"] / nq["Close"].shift(1)).dropna()
intra = np.log(nq["Close"] / nq["Open"]).reindex(on.index)
log(f"NASDAQ 1999-2018: cumulative overnight log return {on.sum():+.2f}, intraday {intra.sum():+.2f}, "
    f"close-to-close {on.sum()+intra.sum():+.2f}")
log(f"   overnight mean {1e4*on.mean():+.2f} bp (t {tstat(on):+.2f}, sd {1e4*on.std():.0f} bp); "
    f"intraday mean {1e4*intra.mean():+.2f} bp (t {tstat(intra):+.2f}, sd {1e4*intra.std():.0f} bp)")
for a, b in [("1999", "2000"), ("2001", "2005"), ("2006", "2010"), ("2011", "2018")]:
    o_, i_ = on[a:b], intra[a:b]
    log(f"   {a}-{b}: overnight {1e4*o_.mean():+.1f} bp/day, intraday {1e4*i_.mean():+.1f} bp/day, "
        f"corr(overnight, intraday) {np.corrcoef(o_, i_)[0,1]:+.2f}")
spo = data.sp500_daily()
on2 = np.log(spo["Open"] / spo["Close"].shift(1))
in2 = np.log(spo["Close"] / spo["Open"])
d2 = pd.DataFrame({"on": on2, "in": in2}).dropna()["2008":"2018"]
d2 = d2[np.abs(d2["on"]) > 1e-9]                     # days with a genuine open
log(f"S&P 500 2008-2018, {len(d2)} days with genuine opens: overnight {d2['on'].sum():+.2f}, "
    f"intraday {d2['in'].sum():+.2f} (overnight t {tstat(d2['on']):+.2f})")

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2))
fig.subplots_adjust(wspace=0.25)
axes[0].bar(tab["month"] - 0.18, tab["mean %"], width=0.36, color=style.SERIES[0], label="raw monthly mean")
axes[0].bar(tab["month"] + 0.18, tab["shrunk mean %"], width=0.36, color=style.SERIES[2], label="James–Stein shrunk")
axes[0].axhline(100 * grand, color=style.INK_2, lw=0.8, ls=":")
axes[0].set_xticks(range(1, 13))
axes[0].set_xticklabels(list("JFMAMJJASOND"))
axes[0].set_ylabel("mean total return (%/month)")
axes[0].set_title("Month-of-year effects, 1926–2018")
axes[0].legend(fontsize=8)
axes[1].plot(on.index, 100 * on.cumsum(), color=style.SERIES[0], label="overnight (close → open)")
axes[1].plot(intra.index, 100 * intra.cumsum(), color=style.SERIES[1], label="intraday (open → close)")
axes[1].axhline(0, color=style.INK_2, lw=0.8)
axes[1].set_ylabel("cumulative log return (%)")
axes[1].set_title("NASDAQ Composite: where the returns happened")
axes[1].legend(fontsize=8)
style.save(fig, "29_calendar.png")

open(style.os.path.join(style.FIG_DIR, "29_output.txt"), "w").write("\n".join(out) + "\n")
