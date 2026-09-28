"""Note 21: how often does a random walk look like a bubble?

LPPL: log p(t) = A + B f + C1 f cos(w ln(tc - t)) + C2 f sin(w ln(tc - t)),
f = (tc - t)^m  (Filimonov-Sornette linearisation: A, B, C1, C2 by OLS).
We fit windows of simulated geometric Brownian motion that happened to rise
strongly, apply Sornette-style bubble qualification filters, and count how
many pass. As a control, the same pipeline on synthetic LPPL bubbles + noise.
"""
import numpy as np
from numba import njit
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(21)
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


@njit(cache=True)
def ols_sse(t, y, tc, m, w):
    n = len(t)
    XtX = np.zeros((4, 4))
    Xty = np.zeros(4)
    for i in range(n):
        dt = tc - t[i]
        f = dt**m
        lg = np.log(dt)
        x = np.array([1.0, f, f * np.cos(w * lg), f * np.sin(w * lg)])
        for a in range(4):
            Xty[a] += x[a] * y[i]
            for b in range(4):
                XtX[a, b] += x[a] * x[b]
    beta = np.linalg.solve(XtX + 1e-12 * np.eye(4), Xty)
    sse = 0.0
    for i in range(n):
        dt = tc - t[i]
        f = dt**m
        lg = np.log(dt)
        pred = beta[0] + beta[1] * f + beta[2] * f * np.cos(w * lg) + beta[3] * f * np.sin(w * lg)
        sse += (y[i] - pred) ** 2
    return sse, beta


@njit(cache=True)
def grid_search(t, y, tcs, ms, ws):
    best = 1e300
    bi = (0.0, 0.0, 0.0)
    for tc in tcs:
        for m in ms:
            for w in ws:
                sse, _ = ols_sse(t, y, tc, m, w)
                if sse < best:
                    best = sse
                    bi = (tc, m, w)
    return best, bi


def fit_lppl(y):
    n = len(y)
    t = np.arange(n, dtype=float) / n            # window mapped to [0, 1)
    tcs = 1.0 + np.linspace(0.002, 0.3, 25)
    ms = np.linspace(0.05, 0.95, 13)
    ws = np.linspace(2.0, 16.0, 22)
    _, (tc0, m0, w0) = grid_search(t, y, tcs, ms, ws)

    def obj(z):
        tc, m, w = z
        if tc <= t[-1] + 1e-4 or not (0.01 < m < 0.99) or not (1.0 < w < 25.0):
            return 1e10
        return ols_sse(t, y, tc, m, w)[0]

    res = minimize(obj, np.array([tc0, m0, w0]), method="Nelder-Mead",
                   options=dict(xatol=1e-5, fatol=1e-10, maxiter=600))
    tc, m, w = res.x
    sse, beta = ols_sse(t, y, tc, m, w)
    A, B, C1, C2 = beta
    C = np.hypot(C1, C2)
    # rescale time to trading days for the oscillation count / tc checks
    return dict(tc=tc, m=m, w=w, A=A, B=B, C=C, sse=sse,
                damping=m * abs(B) / (w * C + 1e-300),
                osc=w / (2 * np.pi) * np.log((tc - t[0]) / (tc - t[-1])),
                r2=1 - sse / np.sum((y - y.mean()) ** 2))


def qualifies(f):
    return (0.1 <= f["m"] <= 0.9 and 6.0 <= f["w"] <= 13.0 and f["B"] < 0
            and f["tc"] <= 1.1 and f["damping"] >= 0.8 and f["osc"] >= 2.5)


def qualifies_loose(f):
    return (0.01 <= f["m"] <= 0.99 and 2.0 <= f["w"] <= 15.0 and f["B"] < 0
            and f["tc"] <= 1.2 and f["damping"] >= 0.5)


N_DAYS = 500
MU, SIG = 0.08, 0.20


def gbm_window(r_):
    x = np.cumsum((MU - 0.5 * SIG**2) / 252 + SIG / np.sqrt(252) * r_.standard_normal(N_DAYS))
    return x


def lppl_bubble(r_, noise=0.01):
    t = np.arange(N_DAYS) / N_DAYS
    tc = 1.0 + r_.uniform(0.01, 0.08)
    m = r_.uniform(0.3, 0.7)
    w = r_.uniform(7, 11)
    B = -0.6
    C = 0.06 * abs(B)             # small enough to satisfy the damping filter
    phi = r_.uniform(0, 2 * np.pi)
    f = (tc - t) ** m
    y = 4.0 + B * f + C * f * np.cos(w * np.log(tc - t) - phi)
    # noise: an AR(1) in log price, so residuals look like real data
    e = np.zeros(N_DAYS)
    for i in range(1, N_DAYS):
        e[i] = 0.95 * e[i - 1] + noise * r_.standard_normal()
    return y + e


# random walks, conditioned on a strong run-up (when people look for bubbles)
results = {"all GBM windows": [], "GBM windows up > 50%": [],
           "LPPL bubble, noise 0.005": [], "LPPL bubble, noise 0.01": []}
n_target = 300
while len(results["all GBM windows"]) < n_target or len(results["GBM windows up > 50%"]) < n_target:
    y = gbm_window(rng)
    up = np.exp(y[-1] - y[0]) - 1
    need_all = len(results["all GBM windows"]) < n_target
    need_up = up > 0.5 and len(results["GBM windows up > 50%"]) < n_target
    if not (need_all or need_up):
        continue
    f = fit_lppl(y)
    f["runup"] = up
    if need_all:
        results["all GBM windows"].append(f)
    if need_up:
        results["GBM windows up > 50%"].append(f)
tc_err = {}
for lvl in (0.005, 0.01):
    key = f"LPPL bubble, noise {lvl:g}"
    errs = []
    for _ in range(n_target):
        r_ = np.random.default_rng(rng.integers(1 << 31))
        tc_true = 1.0 + r_.uniform(0.01, 0.08)
        r_ = np.random.default_rng(rng.integers(1 << 31))
        y = lppl_bubble(r_, noise=lvl)
        f = fit_lppl(y)
        results[key].append(f)
    tc_err[key] = None

log("sample                      n    strict filters   loose filters   median R²")
rates = {}
for k, fs in results.items():
    strict = np.mean([qualifies(f) for f in fs])
    loose = np.mean([qualifies_loose(f) for f in fs])
    rates[k] = (strict, loose)
    log(f"{k:26s} {len(fs):4d}   {100*strict:8.1f}%       {100*loose:8.1f}%      {np.median([f['r2'] for f in fs]):.3f}")

# the predicted critical time for flagged random walks
flag = [f for f in results["GBM windows up > 50%"] if qualifies_loose(f)]
if flag:
    tcd = np.array([(f["tc"] - 1.0) * N_DAYS for f in flag])
    log(f"flagged up-trend random walks: predicted crash {np.median(tcd):.0f} trading days after "
        f"window end (IQR {np.percentile(tcd,25):.0f}-{np.percentile(tcd,75):.0f}); no crash follows "
        f"by construction")

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3))
fig.subplots_adjust(wspace=0.25)
# an example flagged random walk
rr = np.random.default_rng(5)
example = None
for _ in range(400):
    y = gbm_window(rr)
    if np.exp(y[-1] - y[0]) - 1 > 0.5:
        f = fit_lppl(y)
        if qualifies(f) or (example is None and qualifies_loose(f)):
            example = (y, f)
            if qualifies(f):
                break
if example is not None:
    y, f = example
    t = np.arange(N_DAYS) / N_DAYS
    tt = np.linspace(0, f["tc"] - 1e-3, 800)
    _, beta = ols_sse(t, y, f["tc"], f["m"], f["w"])
    ff = (f["tc"] - tt) ** f["m"]
    lg = np.log(f["tc"] - tt)
    fitc = beta[0] + beta[1] * ff + beta[2] * ff * np.cos(f["w"] * lg) + beta[3] * ff * np.sin(f["w"] * lg)
    axes[0].plot(np.arange(N_DAYS), np.exp(y - y[0]), color=style.SERIES[0], lw=1.2, label="geometric random walk")
    axes[0].plot(tt * N_DAYS, np.exp(fitc - y[0]), color=style.SERIES[1], lw=1.5, label="LPPL fit")
    axes[0].axvline(f["tc"] * N_DAYS, color=style.SERIES[1], ls=":", lw=1)
    axes[0].annotate("predicted\ncritical time", (f["tc"] * N_DAYS, np.exp(y - y[0]).min()),
                     fontsize=8, color=style.INK, ha="right")
    axes[0].set_xlabel("trading day")
    axes[0].set_ylabel("price (start = 1)")
    axes[0].set_title(f"A random walk flagged as a bubble (ω = {f['w']:.1f})", fontsize=11)
    axes[0].legend(fontsize=8, loc="upper left")
labels = list(rates.keys())
x = np.arange(len(labels))
axes[1].bar(x - 0.18, [100 * rates[k][0] for k in labels], width=0.34, color=style.SERIES[0], label="strict filters")
axes[1].bar(x + 0.18, [100 * rates[k][1] for k in labels], width=0.34, color=style.SERIES[2], label="loose filters")
axes[1].set_xticks(x)
axes[1].set_xticklabels(["random walk,\nany window", "random walk,\nup > 50%", "LPPL bubble,\nlow noise",
                         "LPPL bubble,\nhigher noise"], fontsize=9)
axes[1].set_ylabel("% of windows flagged as a bubble")
axes[1].set_title("Flag rates: random walks vs LPPL bubbles", fontsize=11)
axes[1].legend(fontsize=8)
style.save(fig, "21_lppl.png")

open(style.os.path.join(style.FIG_DIR, "21_output.txt"), "w").write("\n".join(out) + "\n")
