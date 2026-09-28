"""Note 21 (real data): rolling LPPL bubble detection on daily S&P 500 and NASDAQ,
1999-2018. Windows of 500 trading days, advanced 10 days at a time. For each
window: fit, apply the loose and strict qualification filters from
21_lppl_false_positives.py, then record what happened over the NEXT 60 and 120
trading days (return and maximum drawdown).
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


# import the fitting functions without re-running the simulation study
src = open(style.os.path.join(style.HERE, "21_lppl_false_positives.py")).read()
head = src.split("N_DAYS = 500")[0]
ns = {}
exec(compile(head, "21_lppl_false_positives.py", "exec"), ns)
fit_lppl, qualifies, qualifies_loose = ns["fit_lppl"], ns["qualifies"], ns["qualifies_loose"]
_ols = ns["ols_sse"]


def _safe_ols(t, y, tc, m, w):
    # real price levels make X'X badly conditioned for some (tc, m, w):
    # treat those points as a failed fit rather than crashing
    try:
        return _ols(t, y, tc, m, w)
    except np.linalg.LinAlgError:
        return 1e10, np.zeros(4)


ns["ols_sse"] = _safe_ols

W, STEP = 500, 10
rows = []
for name, df in [("S&P 500", data.sp500_daily()), ("NASDAQ", data.nasdaq_daily())]:
    lp = np.log(df["Close"].values)
    dates = df.index
    for end in range(W, len(lp) - 120, STEP):
        y = lp[end - W:end]
        f = fit_lppl(y)
        fwd60 = lp[end + 59] - lp[end - 1]
        fwd120 = lp[end + 119] - lp[end - 1]
        path = lp[end - 1:end + 120]
        mdd = np.min(path - np.maximum.accumulate(path))
        rows.append(dict(index=name, date=dates[end - 1], strict=qualifies(f), loose=qualifies_loose(f),
                         runup=y[-1] - y[0], fwd60=fwd60, fwd120=fwd120, mdd120=mdd,
                         tc_days=(f["tc"] - 1.0) * W))
res = pd.DataFrame(rows)
for name, g in res.groupby("index"):
    for flag in ["strict", "loose"]:
        fl = g[g[flag]]
        nf = g[~g[flag]]
        up = g["runup"] > np.log(1.5)
        log(f"{name:8s} {flag:6s}: flagged {len(fl):3d}/{len(g)} windows ({100*len(fl)/len(g):4.1f}%); "
            f"next-120d return flagged {100*fl['fwd120'].mean():+5.1f}% vs unflagged {100*nf['fwd120'].mean():+5.1f}%; "
            f"next-120d max drawdown {100*fl['mdd120'].mean():5.1f}% vs {100*nf['mdd120'].mean():5.1f}%")
    fl = g[g["loose"]]
    if len(fl):
        spells = fl["date"].dt.to_period("Q").astype(str).value_counts().sort_index()
        log(f"   {name} loose-flag quarters: " + ", ".join(f"{q}({n})" for q, n in spells.items()))

# figure: price with flagged windows
fig, axes = plt.subplots(2, 1, figsize=(11, 6.0), sharex=True)
fig.subplots_adjust(hspace=0.3)
for ax, (name, df) in zip(axes, [("S&P 500", data.sp500_daily()), ("NASDAQ", data.nasdaq_daily())]):
    ax.plot(df.index, df["Close"], color=style.SERIES[0], lw=1)
    g = res[res["index"] == name]
    for _, r_ in g[g["loose"]].iterrows():
        ax.axvline(r_["date"], color=style.SERIES[1], lw=0.8, alpha=0.5)
    for _, r_ in g[g["strict"]].iterrows():
        ax.axvline(r_["date"], color=style.SERIES[6], lw=1.6, alpha=0.9)
    ax.set_yscale("log")
    ax.set_title(f"{name}: orange = loose LPPL bubble flag, violet = strict flag (at window end)", fontsize=10)
    style.plain_log(ax)
style.save(fig, "21_lppl_real.png")

res.to_csv(style.os.path.join(style.FIG_DIR, "21_lppl_real_windows.csv"), index=False)
open(style.os.path.join(style.FIG_DIR, "21_lppl_real_output.txt"), "w").write("\n".join(out) + "\n")
