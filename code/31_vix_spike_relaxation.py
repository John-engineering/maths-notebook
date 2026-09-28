"""Note 31: how do volatility shocks decay? Event study of VIX spikes, 1990-2026.

Events: VIX rises by more than 40% within 10 trading days to a level above 25,
at least 120 days after the previous event; aligned at the local peak.
Relaxation of the excess log-VIX over its pre-event baseline (median of the 60
days before the run-up), normalised to 1 at the peak, over 250 days.
Fits: exponential exp(-k/tau)  vs  power law (1 + k/k0)^(-beta).
Rough volatility (H ~ 0.1) suggests a power law with beta ~ 1/2 - H ~ 0.4.
"""
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import data
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


v = data.vix_daily()["close"].dropna()
lv = np.log(v.values)
dates = v.index
K = 250
events = []
last = -10_000
i = 70
while i < len(lv) - K - 1:
    run = lv[i] - lv[i - 10]
    if run > np.log(1.4) and v.values[i] > 25 and i - last > 120:
        # peak within the next 30 days (2008 and 2020 kept climbing for weeks)
        j = i + int(np.argmax(lv[i:i + 31]))
        base = np.median(lv[i - 70:i - 10])
        events.append((j, base))
        last = j
        i = j + 1
    else:
        i += 1
log(f"{len(events)} spike events: " + ", ".join(f"{dates[j]:%Y-%m-%d} ({v.values[j]:.0f})" for j, _ in events))

curves = []
pre = []
for j, base in events:
    exc = lv[j:j + K + 1] - base
    curves.append(exc / exc[0])
    pre.append((lv[j - 20:j + 1] - base) / (lv[j] - base))
curves = np.array(curves)
pre = np.array(pre)
mean_c = curves.mean(0)
med_c = np.median(curves, 0)
k = np.arange(K + 1)


def expo(k, tau, c):
    return (1 - c) * np.exp(-k / tau) + c


def power(k, k0, beta):
    return (1 + k / k0) ** (-beta)


def power_c(k, k0, beta, c):
    return (1 - c) * (1 + k / k0) ** (-beta) + c


sel = (k >= 1) & (k <= 120)
target = med_c
pe, _ = curve_fit(expo, k[sel], target[sel], p0=[20, 0.1], bounds=([0.5, -0.5], [2000, 0.9]))
pp, _ = curve_fit(power, k[sel], target[sel], p0=[5, 0.5], bounds=([0.01, 0.01], [500, 5]))
ppc, _ = curve_fit(power_c, k[sel], target[sel], p0=[5, 0.5, 0.0], bounds=([0.01, 0.01, -0.5], [500, 5, 0.9]))
sse = lambda f, p: np.sum((target[sel] - f(k[sel], *p)) ** 2)
# leave-one-event-out: which functional form predicts a held-out spike better?
wins = 0
for e in range(len(curves)):
    others = np.median(np.delete(curves, e, axis=0), 0)
    a, _ = curve_fit(expo, k[sel], others[sel], p0=[20, 0.1], bounds=([0.5, -0.5], [2000, 0.9]))
    b_, _ = curve_fit(power, k[sel], others[sel], p0=[5, 0.5], bounds=([0.01, 0.01], [500, 5]))
    err_e = np.sum((curves[e][sel] - expo(k[sel], *a)) ** 2)
    err_p = np.sum((curves[e][sel] - power(k[sel], *b_)) ** 2)
    wins += err_p < err_e
log(f"leave-one-spike-out: power law predicts the held-out spike better in {wins} of {len(curves)} cases")
log(f"exponential (+ floor): tau = {pe[0]:.1f} days, floor {pe[1]:.2f}; SSE {sse(expo, pe):.4f}")
log(f"power law: k0 = {pp[0]:.2f} days, beta = {pp[1]:.3f}; SSE {sse(power, pp):.4f}")
log(f"power law (+ floor): k0 = {ppc[0]:.2f}, beta = {ppc[1]:.3f}, floor {ppc[2]:.2f}; SSE {sse(power_c, ppc):.4f}")
for kk in [5, 20, 60, 120, 250]:
    log(f"   remaining excess after {kk:3d} days: mean {mean_c[kk]:.2f}, median {med_c[kk]:.2f}")
# asymmetry: rise vs decay time
half_rise = 20 - np.argmax(pre.mean(0) >= 0.5)
half_decay = np.argmax(med_c <= 0.5)
log(f"time-asymmetry: average spike rises from half to peak in {half_rise} days, decays to half in {half_decay} days")

# local slope of log excess vs log time (a power law is a straight line)
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
for c in curves:
    axes[0].plot(k, c, color=style.NEUTRAL, lw=0.5, alpha=0.5)
axes[0].plot(np.arange(-20, 1), pre.mean(0), color=style.SERIES[2], lw=2, label="average run-up")
axes[0].plot(k, med_c, color=style.SERIES[0], lw=2, label=f"median decay ({len(events)} spikes)")
axes[0].axhline(0, color=style.INK_2, lw=0.8)
axes[0].set_xlim(-20, 250)
axes[0].set_ylim(-0.6, 1.3)
axes[0].set_xlabel("trading days from the VIX peak")
axes[0].set_ylabel("excess log-VIX (peak = 1)")
axes[0].set_title("VIX spikes: fast up, slow down")
axes[0].legend(fontsize=8)
kk = k[sel]
axes[1].loglog(kk, med_c[sel], marker="o", ms=2.5, lw=0, color=style.SERIES[0], label="median decay")
axes[1].loglog(kk, expo(kk, *pe), color=style.SERIES[1], lw=1.5, label=f"exponential, τ = {pe[0]:.0f} d")
axes[1].loglog(kk, power(kk, *pp), color=style.SERIES[2], lw=1.5, label=f"power law, β = {pp[1]:.2f}")
axes[1].set_ylim(0.05, 1.2)
axes[1].set_xlabel("days after the peak (log)")
axes[1].set_ylabel("remaining excess (log)")
axes[1].set_title("Relaxation: power law or exponential?")
axes[1].legend(fontsize=8)
style.plain_log(axes[1])
style.plain_log(axes[1], "x")
style.save(fig, "31_vix_relaxation.png")

open(style.os.path.join(style.FIG_DIR, "31_output.txt"), "w").write("\n".join(out) + "\n")
