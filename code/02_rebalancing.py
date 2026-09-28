"""Note 02: where the rebalancing premium comes from, and what it costs.

A. An exact, model-free identity: log(EW / market) = change in log-diversity
   + sum of log(AM/GM) of the period's gross returns. Checked on fat-tailed paths.
B. Two worlds: independent equal-growth stocks (diversity collapses) versus a
   rank-based "Atlas" market (diversity is stationary).
C. Saturation: the per-period gain E log(AM/GM) is s^2/2 only while
   s^2 << 2 ln n; after that it grows like s*sqrt(2 ln n).
D. Transaction costs: band rebalancing, delta* = (3c/16)^(1/3) independent of vol,
   loss ~ c^(2/3). Calendar rebalancing has the same exponent, worse constant.
E. With volatility regimes the band rule adapts on its own; the calendar rule
   doesn't.
"""
import numpy as np
from scipy.special import logsumexp
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(2)
out = []


def log(msg):
    print(msg)
    out.append(msg)


# ============================================================ A. exact identity
n, T = 20, 2000
vols = rng.uniform(0.03, 0.15, n)                       # per-period vols
logR = vols * rng.standard_t(3, (T, n)) / np.sqrt(3) + rng.normal(0, 0.01, n)
R = np.exp(logR)
mu0 = rng.dirichlet(np.ones(n) * 0.7)                   # initial market weights

log_ew = np.log(R.mean(1)).sum()                        # equal weight, rebalanced
caps = mu0 * np.exp(logR.cumsum(0))                     # buy and hold
log_mkt = np.log(caps[-1].sum())                        # market started at 1
muT = caps[-1] / caps[-1].sum()
logG = lambda m: np.log(m).mean()                       # log geometric mean weight
amgm = (np.log(R.mean(1)) - logR.mean(1)).sum()
lhs = log_ew - log_mkt
rhs = logG(muT) - logG(mu0) + amgm
log(f"A. identity check with Student-t(3) returns: lhs={lhs:.12f} rhs={rhs:.12f} "
    f"diff={lhs-rhs:.2e}")

# ============================================================ B. two worlds
def atlas_drift(x, g):
    # every stock drifts down by g except the current smallest, which is
    # pushed up by (n-1)g: total drift zero, ranks mean-revert
    d = np.full(len(x), -g)
    d[np.argmin(x)] = (len(x) - 1) * g
    return d


def world(kind, n=50, years=150, steps_per_year=12, sigma=0.35, g_atlas=0.05,
          seed=0, burn_in_years=600):
    r = np.random.default_rng(seed)
    dt = 1 / steps_per_year
    N = years * steps_per_year
    x = r.normal(0, 1.0, n)                             # log caps
    if kind == "atlas":                                 # start from stationarity
        for _ in range(burn_in_years * steps_per_year):
            x += atlas_drift(x, g_atlas) * dt + sigma * np.sqrt(dt) * r.standard_normal(n)
    mu_start = np.exp(x - logsumexp(x))
    cum_amgm = np.zeros(N)
    cum_G = np.zeros(N)
    rel = np.zeros(N)
    log_ew = 0.0
    log_mkt_cap0 = logsumexp(x)
    acc = 0.0
    for t in range(N):
        drift = atlas_drift(x, g_atlas) if kind == "atlas" else np.zeros(n)
        dx = drift * dt + sigma * np.sqrt(dt) * r.standard_normal(n)
        R_ = np.exp(dx)
        acc += np.log(R_.mean()) - dx.mean()
        log_ew += np.log(R_.mean())
        x = x + dx
        mu = np.exp(x - logsumexp(x))
        cum_amgm[t] = acc
        cum_G[t] = logG(mu) - logG(mu_start)
        rel[t] = log_ew - (logsumexp(x) - log_mkt_cap0)
    return dict(t=np.arange(1, N + 1) * dt, amgm=cum_amgm, G=cum_G, rel=rel,
                final_mu=np.sort(mu)[::-1])


w_ind = world("independent")
w_atl = world("atlas")
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
fig.subplots_adjust(wspace=0.08)
for ax, w, title in [(axes[0], w_ind, "Independent stocks: diversity collapses"),
                     (axes[1], w_atl, "Atlas market: diversity is stationary")]:
    ax.plot(w["t"], w["amgm"], color=style.SERIES[0], label="Σ log(AM/GM): dispersion harvested")
    ax.plot(w["t"], w["G"], color=style.SERIES[1], label="Δ log diversity (geo-mean weight)")
    ax.plot(w["t"], w["rel"], color=style.INK, lw=1.3, label="log(EW / market)")
    ax.axhline(0, color=style.NEUTRAL, lw=1)
    ax.set_title(title)
    ax.set_xlabel("years")
axes[0].set_ylabel("cumulative log points")
axes[0].legend(loc="upper left")
style.save(fig, "02_two_worlds.png")
for name, w in [("independent", w_ind), ("atlas", w_atl)]:
    log(f"B. {name:11s}: harvested {w['amgm'][-1]:+.2f}, diversity {w['G'][-1]:+.2f}, "
        f"EW-vs-market {w['rel'][-1]:+.2f} log pts over 150y")

# capital distribution curves at the end
fig, ax = plt.subplots(figsize=(6.2, 4.2))
ranks = np.arange(1, 51)
ax.loglog(ranks, w_ind["final_mu"], color=style.SERIES[1], marker="o", ms=3,
          label="independent (after 150y)")
ax.loglog(ranks, w_atl["final_mu"], color=style.SERIES[0], marker="o", ms=3,
          label="Atlas (after 150y)")
ax.set_xlabel("rank")
ax.set_ylabel("market weight")
ax.set_title("Capital distribution curves")
style.plain_log(ax)
style.plain_log(ax, "x")
ax.legend()
style.save(fig, "02_capital_distribution.png")

# ============================================================ C. saturation
s2_grid = np.logspace(-2, 1.6, 40)
fig, ax = plt.subplots(figsize=(6.8, 4.4))
for i, nn in enumerate([2, 10, 100, 1000]):
    vals = []
    for s2 in s2_grid:
        s = np.sqrt(s2)
        Y = s * rng.standard_normal((4000 if nn <= 100 else 800, nn))
        vals.append(np.mean(logsumexp(Y, axis=1) - np.log(nn) - Y.mean(1)))
    ax.loglog(s2_grid, vals, color=style.SERIES[i], label=f"n = {nn}")
ax.loglog(s2_grid, s2_grid / 2 * 1.0, color=style.NEUTRAL, ls="--", lw=1.2,
          label="s²/2 (n → ∞)")
ax.set_xlabel("cross-sectional log-return variance per rebalance, s²")
ax.set_ylabel("E log(AM/GM) per rebalance")
ax.set_title("Rebalancing gain saturates once s² ≳ 2 ln n")
style.plain_log(ax)
style.plain_log(ax, "x")
ax.legend()
style.save(fig, "02_saturation.png")

# ============================================================ D. costs
from numba import njit

SIG = 0.30
dt = 1 / 252
YEARS = 20000
steps = int(YEARS / dt)


@njit(cache=True)
def _simulate(cost, deltas, intervals, sig_path, dt, seed):
    np.random.seed(seed)
    nd, ni = len(deltas), len(intervals)
    w_band = np.full(nd, 0.5)
    w_cal = np.full(ni, 0.5)
    lg_band = np.zeros(nd)
    lg_cal = np.zeros(ni)
    cv_band = np.zeros(nd)
    cv_cal = np.zeros(ni)
    lg_ref = 0.0
    for t in range(len(sig_path)):
        s = sig_path[t]
        dx1 = s * np.sqrt(dt) * np.random.randn() - 0.5 * s * s * dt
        dx2 = s * np.sqrt(dt) * np.random.randn() - 0.5 * s * s * dt
        r1, r2 = np.exp(dx1), np.exp(dx2)
        diff = dx1 - dx2
        lg_ref += np.log(0.5 * r1 + 0.5 * r2)
        for i in range(nd):
            w = w_band[i]
            # control variate: (w - 1/2)(dX1 - dX2) has mean zero and removes
            # the martingale noise that otherwise swamps a few-bp/yr signal
            cv_band[i] += (w - 0.5) * diff
            g = w * r1 + (1 - w) * r2
            lg_band[i] += np.log(g)
            w = w * r1 / g
            d = deltas[i]
            if w > 0.5 + d:
                lg_band[i] += np.log1p(-2 * cost * (w - 0.5 - d))
                w = 0.5 + d
            elif w < 0.5 - d:
                lg_band[i] += np.log1p(-2 * cost * (0.5 - d - w))
                w = 0.5 - d
            w_band[i] = w
        for i in range(ni):
            w = w_cal[i]
            cv_cal[i] += (w - 0.5) * diff
            g = w * r1 + (1 - w) * r2
            lg_cal[i] += np.log(g)
            w = w * r1 / g
            if (t + 1) % intervals[i] == 0:
                lg_cal[i] += np.log1p(-2 * cost * abs(w - 0.5))
                w = 0.5
            w_cal[i] = w
    return lg_ref, lg_band, lg_cal, cv_band, cv_cal


def run_rules(cost, deltas, intervals, sig_path, seed=11):
    """Two uncorrelated equal-growth assets, target 50/50. Returns the annual
    growth shortfall versus the frictionless daily-rebalanced portfolio."""
    lg_ref, lb, lc, cb, cc = _simulate(cost, np.asarray(deltas, float),
                                       np.asarray(intervals, np.int64),
                                       np.asarray(sig_path, float), dt, seed)
    years = len(sig_path) * dt
    return (lg_ref - lb + cb) / years, (lg_ref - lc + cc) / years


def cal_loss(D, c, sig):
    """Calendar rule with interval D years: tracking loss sig^4 D / 16, plus
    cost 2c E|x_D| / D with x_D ~ N(0, sig^2 D / 8)."""
    sw = sig / (2 * np.sqrt(2))
    return sig**4 * D / 16 + 2 * c * sw * np.sqrt(2 / np.pi) / np.sqrt(D)


def pred_cal(c, sig):
    return (8 * c / (np.sqrt(np.pi) * sig**3)) ** (2 / 3)


costs = [0.0005, 0.002, 0.008]
deltas = np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.13,
                   0.17, 0.22, 0.3])
intervals = np.array([5, 10, 21, 42, 63, 126, 252, 504, 1008, 2016])
const_sig = np.full(steps, SIG)
loss_band, loss_cal = {}, {}
for c in costs:
    loss_band[c], loss_cal[c] = run_rules(c, deltas, intervals, const_sig)
    b = np.argmin(loss_band[c])
    k = np.argmin(loss_cal[c])
    pred_d = (3 * c / 16) ** (1 / 3)
    pred_loss = SIG ** 2 * pred_d ** 2
    log(f"D. c={c:.4f}: best band ±{deltas[b]:.3f} (theory ±{pred_d:.3f}), "
        f"loss {1e4*loss_band[c][b]:.2f}bp/yr (theory {1e4*pred_loss:.2f}); "
        f"best calendar {intervals[k]}d (theory {252*pred_cal(c, SIG):.0f}d), "
        f"loss {1e4*loss_cal[c][k]:.2f}bp/yr "
        f"(theory {1e4*cal_loss(pred_cal(c, SIG), c, SIG):.2f})")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
fig.subplots_adjust(wspace=0.08)
for i, c in enumerate(costs):
    axes[0].plot(deltas, 1e4 * loss_band[c], marker="o", ms=4, lw=0,
                 color=style.SERIES[i], label=f"c = {1e4*c:g} bp")
    d = np.geomspace(deltas[0], deltas[-1], 200)
    axes[0].plot(d, 1e4 * SIG**2 * (d**2 / 3 + c / (8 * d)), color=style.SERIES[i],
                 lw=1.5)
    axes[1].plot(intervals / 252, 1e4 * loss_cal[c], marker="o", ms=4, lw=0,
                 color=style.SERIES[i], label=f"c = {1e4*c:g} bp")
    D = np.geomspace(intervals[0], intervals[-1], 300) / 252
    axes[1].plot(D, 1e4 * cal_loss(D, c, SIG), color=style.SERIES[i], lw=1.5)
axes[0].set_xlabel("no-trade band half-width δ (portfolio weight)")
axes[0].set_ylabel("growth lost vs frictionless (bp/yr)")
axes[0].set_title("Band rule: dots simulated, lines theory")
axes[1].set_xlabel("rebalancing interval (years)")
axes[1].set_title("Calendar rule: dots simulated, lines theory")
for ax in axes:
    ax.set_xscale("log")
    ax.set_yscale("log")
    style.plain_log(ax)
    style.plain_log(ax, "x")
    ax.legend(loc="upper center")
style.save(fig, "02_costs.png")

# ============================================================ E. vol regimes
# Two-state Markov volatility: 15% calm / 60% stressed, stressed ~ 1/5 of time.
@njit(cache=True)
def _regimes(u, p_enter, p_exit):
    state = np.zeros(len(u), np.int64)
    s_ = 0
    for t in range(len(u)):
        if s_ == 0 and u[t] < p_enter:
            s_ = 1
        elif s_ == 1 and u[t] < p_exit:
            s_ = 0
        state[t] = s_
    return state


rs = np.random.default_rng(5)
state = _regimes(rs.random(steps), 1 / (252 * 4), 1 / 252)
sig_regime = np.where(state == 1, 0.60, 0.15)
rms = np.sqrt(np.mean(sig_regime**2))
log(f"E. regime vol: stressed fraction {state.mean():.2f}, rms vol {rms:.3f}")
c = 0.002
lb_const, lc_const = run_rules(c, deltas, intervals, np.full(steps, rms), seed=21)
lb_reg, lc_reg = run_rules(c, deltas, intervals, sig_regime, seed=21)
log(f"E. constant vol {rms:.2f}: best band {1e4*lb_const.min():.2f}bp, best calendar "
    f"{1e4*lc_const.min():.2f}bp  (ratio {lc_const.min()/lb_const.min():.2f})")
log(f"E. regime vol:            best band {1e4*lb_reg.min():.2f}bp, best calendar "
    f"{1e4*lc_reg.min():.2f}bp  (ratio {lc_reg.min()/lb_reg.min():.2f})")
log(f"E. regime vol: argmin band {deltas[np.argmin(lb_reg)]}, const {deltas[np.argmin(lb_const)]}; "
    f"argmin calendar {intervals[np.argmin(lc_reg)]}d vs const {intervals[np.argmin(lc_const)]}d")

open(style.os.path.join(style.FIG_DIR, "02_output.txt"), "w").write("\n".join(out) + "\n")
