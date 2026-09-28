"""Note 06: the next price move as a harmonic measure.

Bid and ask queue sizes (b, a) move like a correlated 2-D random walk in the
positive quadrant. The mid-price ticks up when the ask queue empties first,
down when the bid does. Whitening the correlation maps the quadrant to a wedge
of angle alpha = arccos(-rho), where

    P(up | a, b) = psi / alpha,  cos psi = (a - rho b) / sqrt(a^2 - 2 rho a b + b^2)

and the waiting time for the next price change has a power-law tail
P(tau > t) ~ t^(-pi / (2 alpha)): finite mean only if rho < 0.
"""
import numpy as np
from numba import njit
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg)
    out.append(msg)


def p_up(a, b, rho):
    alpha = np.arccos(-rho)
    c = (a - rho * b) / np.sqrt(a * a - 2 * rho * a * b + b * b)
    return np.arccos(np.clip(c, -1, 1)) / alpha


def expected_tau(a, b, rho, var_rate):
    """E[time to exit] for rho < 0, in the diffusion limit. var_rate is the
    variance per unit time of each queue."""
    alpha = np.arccos(-rho)
    if np.cos(alpha) <= 0:
        return np.inf
    # whitened polar coordinates: r^2 = x' Sigma^{-1} x, angle from the a-axis
    r2 = (a * a - 2 * rho * a * b + b * b) / (1 - rho**2) / var_rate
    phi = np.arccos((a - rho * b) / np.sqrt(a * a - 2 * rho * a * b + b * b))
    return 0.5 * r2 * (np.cos(2 * phi - alpha) / np.cos(alpha) - 1)


@njit(cache=True)
def race(a0, b0, c, sign, n_paths, max_events, seed):
    """Queues with unit-rate independent +/-1 moves on each side plus correlated
    events at rate c each: (+1,+1),(-1,-1) if sign>0, (+1,-1),(-1,+1) if sign<0.
    Returns up flags (1 up, 0 down, -1 undecided) and exit times."""
    np.random.seed(seed)
    total = 4.0 + 2.0 * c
    up = np.empty(n_paths, np.int64)
    tau = np.empty(n_paths)
    for p in range(n_paths):
        a, b, t = a0, b0, 0.0
        res = -1
        for _ in range(max_events):
            t += np.random.exponential(1.0 / total)
            u = np.random.random() * total
            if u < 1.0:
                a += 1
            elif u < 2.0:
                a -= 1
            elif u < 3.0:
                b += 1
            elif u < 4.0:
                b -= 1
            elif u < 4.0 + c:
                a += 1
                b += sign
            else:
                a -= 1
                b -= sign
            if a <= 0 and b <= 0:          # simultaneous depletion: coin flip
                res = 1 if np.random.random() < 0.5 else 0
                break
            if a <= 0:
                res = 1
                break
            if b <= 0:
                res = 0
                break
        up[p] = res
        tau[p] = t
    return up, tau


def params(rho):
    c = abs(rho) / (1 - abs(rho))
    sign = 1 if rho >= 0 else -1
    var_rate = 2 * (1 + c)
    return c, sign, var_rate


# ------------------------------------------------------------- 1. P(up) check
sizes = [1, 2, 3, 5, 8, 13, 21]
rhos = [-0.5, 0.0, 0.5]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.0), sharey=True)
fig.subplots_adjust(wspace=0.08)
worst = {}
for ax, rho in zip(axes, rhos):
    c, sign, _ = params(rho)
    ratios, sim, errs = [], [], []
    for a in sizes:
        for b in sizes:
            up, _ = race(a, b, c, sign, 3000, 2_000_000, 17 + 31 * a + b)
            dec = up >= 0
            sim.append(up[dec].mean())
            ratios.append(np.arctan2(b, a))
            errs.append(up[dec].mean() - p_up(a, b, rho))
    worst[rho] = np.max(np.abs(errs))
    th = np.linspace(0.001, np.pi / 2 - 0.001, 200)
    ax.plot(th * 180 / np.pi, p_up(np.cos(th), np.sin(th), rho), color=style.SERIES[0],
            label="harmonic measure of the wedge")
    ax.plot(np.array(ratios) * 180 / np.pi, sim, marker="o", ms=4, lw=0,
            color=style.SERIES[1], label="simulated discrete queues")
    ax.plot(th * 180 / np.pi, np.sin(th) / (np.sin(th) + np.cos(th)), color=style.NEUTRAL,
            ls="--", lw=1.2, label="linear imbalance b/(a+b)")
    ax.set_title(f"ρ = {rho:+.1f}   (wedge angle {np.degrees(np.arccos(-rho)):.0f}°)")
    ax.set_xlabel("queue angle arctan(bid / ask), degrees")
    log(f"1  rho={rho:+.1f}: max |sim - theory| over 49 (a,b) pairs = {worst[rho]:.3f}")
axes[0].set_ylabel("P(next mid move is up)")
axes[0].legend(fontsize=8, loc="upper left")
style.save(fig, "06_p_up.png")

# small queues are where discreteness bites: show error by min(a,b)
for rho in rhos:
    c, sign, _ = params(rho)
    line = []
    for a, b in [(1, 1), (1, 3), (3, 1), (2, 5), (10, 10), (5, 20)]:
        up, _ = race(a, b, c, sign, 20000, 2_000_000, 7 * a + b)
        dec = up >= 0
        line.append(f"({a},{b}) sim {up[dec].mean():.3f} th {p_up(a, b, rho):.3f}")
    log(f"1b rho={rho:+.1f}: " + "; ".join(line))

# ------------------------------------------------------------- 2. waiting times
fig, ax = plt.subplots(figsize=(7.2, 4.4))
taus = {}
for i, rho in enumerate(rhos):
    c, sign, var_rate = params(rho)
    up, tau = race(10, 10, c, sign, 40000, 20_000_000, 101 + i)
    taus[rho] = tau
    t_sorted = np.sort(tau)
    surv = 1 - np.arange(len(t_sorted)) / len(t_sorted)
    ax.loglog(t_sorted, surv, color=style.SERIES[i], label=f"ρ = {rho:+.1f}")
    alpha = np.arccos(-rho)
    k = np.pi / (2 * alpha)
    # reference slope anchored at the 90th percentile
    t0 = np.quantile(tau, 0.9)
    tt = np.geomspace(t0, t_sorted[-1], 50)
    ax.loglog(tt, 0.1 * (tt / t0) ** (-k), color=style.SERIES[i], ls="--", lw=1)
    # fitted tail exponent via Hill estimator on the top 2%
    top = t_sorted[-int(0.02 * len(t_sorted)):]
    hill = 1 / np.mean(np.log(top / top[0]))
    msg = (f"2  rho={rho:+.1f}: wedge {np.degrees(alpha):.0f} deg, tail exponent theory "
           f"{k:.3f}, Hill {hill:.3f}; undecided {np.mean(up < 0):.4f}")
    if rho < 0:
        msg += (f"; E[tau] sim {tau.mean():.1f} vs diffusion theory "
                f"{expected_tau(10, 10, rho, var_rate):.1f} "
                f"(with half-tick boundary correction {expected_tau(10.5, 10.5, rho, var_rate):.1f})")
    log(msg)
ax.set_xlabel("time until the next price change (event-rate units, log)")
ax.set_ylabel("P(wait > t) (log)")
ax.set_title("Waiting times: power laws with exponent π / (2α)")
ax.legend()
style.plain_log(ax)
style.save(fig, "06_waiting_times.png")

# ------------------------------------------------------------- 2b. subdiffusion
# If the queues regenerate at (10,10) after every price change, the number of
# price changes N(t) is a renewal process. With tail exponent k < 1 renewal
# theory gives E N(t) ~ t^k, so the price variance grows like t^k: subdiffusive.
fig, ax = plt.subplots(figsize=(7.2, 4.2))
grid_t = np.geomspace(1e2, 3e5, 25)
r_ = np.random.default_rng(8)
for i, rho in enumerate(rhos):
    tau = taus[rho]
    counts = np.zeros(len(grid_t))
    n_paths = 3000
    for _ in range(n_paths):
        arr = np.cumsum(r_.choice(tau, 4000))
        counts += np.searchsorted(arr, grid_t, side="right")
    counts /= n_paths
    slope = np.polyfit(np.log(grid_t[10:]), np.log(counts[10:]), 1)[0]
    k = min(1.0, np.pi / (2 * np.arccos(-rho)))
    ax.loglog(grid_t, counts, color=style.SERIES[i], marker="o", ms=3,
              label=f"ρ = {rho:+.1f}: fitted slope {slope:.2f}, tail exponent κ = {np.pi / (2 * np.arccos(-rho)):.2f}")
    log(f"2b rho={rho:+.1f}: E N(t) log-log slope {slope:.3f}, theory min(1, pi/2alpha) = {k:.3f}")
ax.set_xlabel("elapsed time t (log)")
ax.set_ylabel("expected number of price changes (log)")
ax.set_title("Price variance ∝ E N(t): subdiffusive unless ρ < 0")
ax.legend()
style.plain_log(ax)
style.plain_log(ax, "x")
style.save(fig, "06_subdiffusion.png")

# ------------------------------------------------------------- 3. geometry figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
fig.subplots_adjust(wspace=0.25)
rho = 0.5
alpha = np.arccos(-rho)
levels = np.linspace(0.1, 0.9, 9)
th = np.linspace(0.0005, np.pi / 2 - 0.0005, 4000)
pv = p_up(np.cos(th), np.sin(th), rho)
for lev in levels:
    col = style.BLUE_RAMP[1 + int(round(lev * 5))]
    t_l = th[np.argmin(np.abs(pv - lev))]
    axes[0].plot([0, 3 * np.cos(t_l)], [0, 3 * np.sin(t_l)], color=col, lw=1.3)
    # in whitened coordinates the level-lev ray sits at angle lev * alpha
    axes[1].plot([0, 3 * np.cos(lev * alpha)], [0, 3 * np.sin(lev * alpha)], color=col, lw=1.3)
    axes[1].annotate(f"{lev:.1f}", (3.15 * np.cos(lev * alpha), 3.15 * np.sin(lev * alpha)),
                     fontsize=7, color=style.INK_2, ha="center", va="center")
axes[0].plot([0, 3], [0, 0], color=style.INK, lw=2.2)
axes[0].plot([0, 0], [0, 3], color=style.INK, lw=2.2)
axes[0].annotate("bid empties → down", (1.6, 0.06), fontsize=8, color=style.INK_2)
axes[0].annotate("ask empties → up", (0.06, 3.02), fontsize=8, color=style.INK_2)
axes[1].plot([0, 3], [0, 0], color=style.INK, lw=2.2)
axes[1].plot([0, 3 * np.cos(alpha)], [0, 3 * np.sin(alpha)], color=style.INK, lw=2.2)
axes[1].annotate("bid empties → down", (1.6, -0.25), fontsize=8, color=style.INK_2)
axes[1].annotate("ask empties → up", (3 * np.cos(alpha) - 0.3, 3 * np.sin(alpha) + 0.15),
                 fontsize=8, color=style.INK_2)
axes[0].set_xlabel("ask queue a")
axes[0].set_ylabel("bid queue b")
axes[0].set_title("Quadrant: rays of constant P(up), ρ = +0.5")
axes[1].set_title(f"Whitened: a {np.degrees(alpha):.0f}° wedge, P(up) = angle / α")
axes[0].set_xlim(-0.2, 3.4)
axes[0].set_ylim(-0.2, 3.4)
axes[1].set_xlim(-2.0, 3.5)
axes[1].set_ylim(-0.5, 3.4)
for ax in axes:
    ax.set_aspect("equal")
    ax.grid(False)
axes[1].set_xticks([])
axes[1].set_yticks([])
axes[1].spines["left"].set_visible(False)
axes[1].spines["bottom"].set_visible(False)
style.save(fig, "06_geometry.png")

open(style.os.path.join(style.FIG_DIR, "06_output.txt"), "w").write("\n".join(out) + "\n")
