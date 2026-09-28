"""Note 07: luck, track records and continuous monitoring.

1. Arcsine law: fraction of time a zero-skill manager is ahead of benchmark.
2. Sampling to a foregone conclusion: P(t-stat ever exceeds 2 under monitoring),
   with a Pickands/OU approximation.
3. Always-valid boundary from Robbins' normal-mixture martingale + Ville's
   inequality; false-alarm rate and detection times for real skill.
"""
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(7)
out = []


def log(msg):
    print(msg)
    out.append(msg)


MONTHS = 360
N_MGR = 100_000

# ------------------------------------------------------------------ 1. arcsine
x = rng.standard_normal((N_MGR, MONTHS)).astype(np.float32)
S = np.cumsum(x, axis=1)
frac_ahead = (S > 0).mean(1)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
fig.subplots_adjust(wspace=0.25)
axes[0].hist(frac_ahead, bins=50, density=True, color=style.SERIES[0], alpha=0.85,
             label="100,000 zero-skill managers")
u = np.linspace(0.002, 0.998, 400)
axes[0].plot(u, 1 / (np.pi * np.sqrt(u * (1 - u))), color=style.SERIES[1],
             label="arcsine density 1/(π√(x(1−x)))")
axes[0].set_ylim(0, 6)
axes[0].set_xlabel("fraction of 30 years spent ahead of the benchmark")
axes[0].set_ylabel("density")
axes[0].set_title("Zero skill looks like consistency")
axes[0].legend(fontsize=8, loc="upper center")
p_ext = np.mean((frac_ahead > 0.95) | (frac_ahead < 0.05))
p_ext90 = np.mean((frac_ahead > 0.9) | (frac_ahead < 0.1))
th_ext = 2 * (2 / np.pi) * np.arcsin(np.sqrt(0.05))
th_ext90 = 2 * (2 / np.pi) * np.arcsin(np.sqrt(0.1))
log(f"1  P(ahead >95% or <5% of the time) sim {p_ext:.3f} theory {th_ext:.3f}; "
    f">90%/<10%: sim {p_ext90:.3f} theory {th_ext90:.3f}")
# last time the cumulative relative return changed sign
sign = np.sign(S)
changed = sign[:, 1:] != sign[:, :-1]
last = np.where(changed.any(1), MONTHS - 1 - np.argmax(changed[:, ::-1], axis=1), 0)
p_last = np.mean((last < 0.1 * MONTHS) & (S[:, -1] > 0))
log(f"1  P(continuously ahead for the last 27 of 30 years) sim {p_last:.3f}, "
    f"theory {(1/np.pi)*np.arcsin(np.sqrt(0.1)):.3f}")
del x

# ------------------------------------------------------------------ 2. monitoring
start = 12
t = np.arange(1, MONTHS + 1)
Z = S / np.sqrt(t)
ever2 = np.maximum.accumulate(Z[:, start - 1:] > 2.0, axis=1).mean(0)
ever2_abs = np.maximum.accumulate(np.abs(Z[:, start - 1:]) > 2.0, axis=1).mean(0)
horizon = t[start - 1:] / 12
# continuous-time approximation: stationary OU in log-time, Pickands constant
L = np.log(t[start - 1:] / start)
c = 2.0
pick = 1 - (1 - norm.sf(c)) * np.exp(-L * c * norm.pdf(c) / 2)
log(f"2  P(t>2 at some monthly check, years 1..30): sim {ever2[-1]:.3f}; "
    f"continuous-monitoring Pickands approx {pick[-1]:.3f}; |t|>2 two-sided {ever2_abs[-1]:.3f}")
for yrs in [3, 5, 10, 20, 30]:
    i = yrs * 12 - start
    log(f"   by year {yrs:2d}: one-sided {ever2[i]:.3f}   two-sided {ever2_abs[i]:.3f}")

# ------------------------------------------------------------------ 3. always-valid
def mixture_boundary(tt, alpha=0.05, rho=12.0):
    """|S_t| >= sqrt((t+rho) (log((t+rho)/rho) + 2 log(1/alpha))): two-sided,
    P(ever crossed | no skill) <= alpha by Ville's inequality."""
    return np.sqrt((tt + rho) * (np.log((tt + rho) / rho) + 2 * np.log(1 / alpha)))


bnd = mixture_boundary(t)
crossed = np.abs(S) >= bnd
fa = np.maximum.accumulate(crossed, axis=1).mean(0)
log(f"3  always-valid boundary: false alarms by year 30 = {fa[-1]:.4f} (guarantee <= 0.05)")
for yrs in [1, 5, 10, 30]:
    log(f"   t-stat threshold at year {yrs:2d}: {bnd[yrs*12-1]/np.sqrt(yrs*12):.2f}")
del S, Z, crossed

fig2_x = horizon
axes[1].plot(fig2_x, 100 * ever2, color=style.SERIES[0], label="naive: stop when t > 2 (one-sided)")
axes[1].plot(fig2_x, 100 * ever2_abs, color=style.SERIES[1], label="naive: stop when |t| > 2")
axes[1].plot(fig2_x, 100 * pick, color=style.SERIES[0], ls="--", lw=1.2,
             label="continuous-time approximation")
axes[1].plot(fig2_x, 100 * fa[start - 1:], color=style.SERIES[2],
             label="always-valid mixture boundary")
axes[1].axhline(5, color=style.NEUTRAL, lw=1, ls=":")
axes[1].set_xlabel("years of monthly monitoring")
axes[1].set_ylabel("% of zero-skill managers flagged")
axes[1].set_title("Watching a t-stat 'finds' skill that isn't there")
axes[1].legend(fontsize=8, loc="upper left")
style.save(fig, "07_arcsine_monitoring.png")

# detection times for real skill
fig, ax = plt.subplots(figsize=(7.4, 4.2))
te = 0.05                                             # 5% tracking error, monthly units below
for i, sr in enumerate([0.25, 0.5, 1.0]):
    mu = sr / np.sqrt(12)                              # monthly mean in tracking-error units
    y = rng.standard_normal((20_000, MONTHS)).astype(np.float32) + mu
    Sy = np.cumsum(y, axis=1)
    hit = Sy >= bnd                                     # upper side only for detection
    first = np.where(hit.any(1), np.argmax(hit, axis=1) + 1, np.inf)
    detected = np.array([np.mean(first <= m) for m in t])
    ax.plot(t / 12, 100 * detected, color=style.SERIES[i], label=f"information ratio {sr}")
    med = np.median(first) / 12
    med_s = f"{med:.1f}y" if np.isfinite(med) else "beyond 30y"
    fixed_power = norm.cdf(sr * np.sqrt(30) - norm.isf(0.025))
    log(f"3  IR {sr}: always-valid detection within 30y {100*detected[-1]:.1f}%, median time "
        f"{med_s}; a single test at year 30 (t>1.96) would have power {100*fixed_power:.1f}%")
ax.set_xlabel("years")
ax.set_ylabel("% of skilled managers detected")
ax.set_title("Always-valid detection: how long real skill takes to show")
ax.legend()
style.save(fig, "07_detection.png")

open(style.os.path.join(style.FIG_DIR, "07_output.txt"), "w").write("\n".join(out) + "\n")
