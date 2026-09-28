"""Note 18: how long until you notice your edge is gone?

Daily strategy returns (in units of their volatility) are N(d, 1) while the edge
is alive (d = SR / sqrt(252)) and N(d_after, 1) after it dies. CUSUM for the
change d -> d_after:
    W_t = max(0, W_{t-1} + (d_after - d)(x_t - (d + d_after)/2)),   alarm at W >= h.
For d_after = 0 this is d * (drawdown of the cumulative P&L minus half its
expected drift). Lorden/Moustakides: CUSUM is optimal; delay ~ log(ARL0) / KL,
KL = (d - d_after)^2 / 2 per day.
"""
import numpy as np
from numba import njit
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


@njit(cache=True)
def drawdown_run_length(mu, h, n_paths, max_steps, seed):
    """Days until the plain drawdown of cumulative P&L first reaches h."""
    np.random.seed(seed)
    res = np.empty(n_paths)
    for p in range(n_paths):
        cum = 0.0
        peak = 0.0
        t = 0
        while t < max_steps:
            cum += mu + np.random.randn()
            t += 1
            if cum > peak:
                peak = cum
            if peak - cum >= h:
                break
        res[p] = t
    return res


@njit(cache=True)
def llr_cusum(mu, d, d_after, h, n_paths, max_steps, seed):
    np.random.seed(seed)
    res = np.empty(n_paths)
    for p in range(n_paths):
        W = 0.0
        t = 0
        while t < max_steps:
            x = mu + np.random.randn()
            t += 1
            l = (d_after - d) * x - 0.5 * (d_after**2 - d**2)
            W = max(0.0, W + l)
            if W >= h:
                break
        res[p] = t
    return res


def arl_and_delay(sr, sr_after, thresholds, rule, n=800):
    d, da = sr / np.sqrt(252), sr_after / np.sqrt(252)
    arl, delay = [], []
    for k, h in enumerate(thresholds):
        if rule == 0:
            rl0 = llr_cusum(d, d, da, h, n, 252 * 400, 11 + k)
            rl1 = llr_cusum(da, d, da, h, n, 252 * 400, 911 + k)   # change at time 0 (worst case start W=0)
        else:
            rl0 = drawdown_run_length(d, h, n, 252 * 400, 11 + k)
            rl1 = drawdown_run_length(da, h, n, 252 * 400, 911 + k)
        arl.append(np.mean(rl0) / 252)
        delay.append(np.mean(rl1) / 252)
    return np.array(arl), np.array(delay)


def interp_delay(arl, delay, target):
    """log-log interpolation of delay at a given ARL (years)."""
    o = np.argsort(arl)
    return float(np.exp(np.interp(np.log(target), np.log(arl[o]), np.log(delay[o]))))


def brownian(h, tau):
    """Continuous-time CUSUM: ARL0 and detection delay in units of tau = 1/KL."""
    return tau * (np.exp(h) - h - 1), tau * (h - 1 + np.exp(-h))


fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
hs = np.linspace(0.5, 6.0, 12)
for i, sr in enumerate([0.5, 1.0, 2.0]):
    d = sr / np.sqrt(252)
    tau = 1 / (d**2 / 2) / 252                     # years
    arl, delay = arl_and_delay(sr, 0.0, hs, 0)
    axes[0].plot(arl, delay, marker="o", ms=4, lw=0, color=style.SERIES[i], label=f"Sharpe {sr} → 0")
    hh = np.linspace(0.3, 7, 200)
    A, D = brownian(hh, tau)
    axes[0].plot(A, D, color=style.SERIES[i], lw=1.5)
    arl_n, delay_n = arl_and_delay(sr, -sr, hs, 0)
    axes[0].plot(arl_n, delay_n, marker="s", ms=4, lw=0, color=style.SERIES[i], alpha=0.5,
                 label=f"Sharpe {sr} → −{sr}")
    A2, D2 = brownian(hh, tau / 4)
    axes[0].plot(A2, D2, color=style.SERIES[i], lw=1, ls="--", alpha=0.7)
    for target in (5, 20):
        log(f"SR {sr}: timescale 2/SR^2 = {tau:.1f}y; at 1 false alarm per {target}y, detecting death "
            f"to 0 takes {interp_delay(arl, delay, target):.1f}y (Brownian theory "
            f"{np.interp(target, A, D):.1f}y); death to -SR takes {interp_delay(arl_n, delay_n, target):.1f}y")
axes[0].set_xscale("log")
axes[0].set_yscale("log")
axes[0].set_xlim(0.5, 300)
axes[0].set_ylim(0.05, 40)
axes[0].set_xlabel("average years between false alarms (edge still alive)")
axes[0].set_ylabel("average years to detect a dead edge")
axes[0].set_title("CUSUM: dots simulated, lines Brownian theory")
axes[0].legend(fontsize=7.5, ncol=2, loc="upper left")
style.plain_log(axes[0])
style.plain_log(axes[0], "x")

# plain drawdown rule vs CUSUM for SR = 1 -> 0
sr = 1.0
arl_c, del_c = arl_and_delay(sr, 0.0, hs, 0)
dd_levels = np.linspace(8, 80, 12)             # drawdown in daily-sd units
arl_d, del_d = arl_and_delay(sr, 0.0, dd_levels, 1)
axes[1].plot(arl_c, del_c, marker="o", ms=4, color=style.SERIES[0], label="CUSUM = drawdown of (P&L − ½ expected drift)")
axes[1].plot(arl_d, del_d, marker="o", ms=4, color=style.SERIES[1], label="plain drawdown of P&L")
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlabel("average years between false alarms")
axes[1].set_ylabel("average years to detect")
axes[1].set_title("Sharpe 1 → 0: CUSUM vs a plain drawdown stop")
axes[1].legend(fontsize=8)
style.plain_log(axes[1])
style.plain_log(axes[1], "x")
style.save(fig, "18_dead_edge.png")
for target in (5, 20, 50):
    dd_needed = np.exp(np.interp(np.log(target), np.log(arl_d), np.log(dd_levels))) / np.sqrt(252)
    log(f"SR 1 -> 0 at 1 false alarm per {target}y: CUSUM delay {interp_delay(arl_c, del_c, target):.1f}y, "
        f"plain drawdown delay {interp_delay(arl_d, del_d, target):.1f}y "
        f"(drawdown trigger = {dd_needed:.2f} annual vols)")

open(style.os.path.join(style.FIG_DIR, "18_output.txt"), "w").write("\n".join(out) + "\n")
