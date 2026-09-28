"""Note 09: when does volatility targeting raise the Sharpe ratio?

Model: r_t = mu_t + sigma_t eps_t with mu_t = kappa sigma_t^p and log sigma_t
Gaussian with standard deviation s. Exposure w_t = sigma_hat_t^(-m).
Closed form (lognormal algebra):
    SR_m / SR_opt = exp(-s^2 (m b + p - 2)^2 / 2 - m^2 e^2 / 2)
where log sigma_hat = b log sigma + eta, eta ~ N(0, e^2). With a perfect
forecast (b=1, e=0) vol targeting (m=1) beats constant exposure (m=0) iff p < 3/2.
"""
import numpy as np
from scipy.signal import lfilter
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(9)
out = []


def log(msg):
    print(msg)
    out.append(msg)


def ratio(m, p, s, b=1.0, e=0.0):
    return np.exp(-0.5 * s**2 * (m * b + p - 2) ** 2 - 0.5 * m**2 * e**2)


# ------------------------------------------------------- simulate a SV world
# daily log-vol AR(1), persistence 0.985, stationary sd s; kappa set so that
# the unconditional Sharpe is modest; returns Gaussian given vol.
N = 4_000_000
phi = 0.985
s = 0.45
eta = rng.standard_normal(N) * s * np.sqrt(1 - phi**2)
logsig = lfilter([1], [1, -phi], eta) + np.log(0.16 / np.sqrt(252))
sig = np.exp(logsig)
eps = rng.standard_normal(N)
s_emp = np.std(logsig)


def sharpe(pos, r):
    pnl = pos * r
    return pnl.mean() / pnl.std() * np.sqrt(252)


ms = np.linspace(-0.5, 3.0, 29)
fig, ax = plt.subplots(figsize=(7.6, 4.4))
for i, p in enumerate([0.0, 1.0, 1.5, 2.0]):
    mu = 0.02 * (sig / sig.mean()) ** p * sig.mean()      # kappa sigma^p, SR ~ 0.3
    r = mu + sig * eps
    sr_opt = sharpe(mu / sig**2, r)
    sims = np.array([sharpe(sig ** (-m), r) for m in ms])
    ax.plot(ms, ratio(ms, p, s_emp), color=style.SERIES[i], label=f"p = {p:g}")
    ax.plot(ms, sims / sr_opt, color=style.SERIES[i], lw=0, marker="o", ms=3.5)
    log(f"p={p:3.1f}: SR const {sharpe(np.ones(N), r):.3f}, vol-target {sharpe(1/sig, r):.3f}, "
        f"var-target {sharpe(1/sig**2, r):.3f}, optimal {sr_opt:.3f}; "
        f"theory ratio vol/const = {ratio(1,p,s_emp)/ratio(0,p,s_emp):.3f}, sim "
        f"{sharpe(1/sig, r)/sharpe(np.ones(N), r):.3f}")
for m_, lab in [(0, "constant\nexposure"), (1, "vol\ntarget"), (2, "variance\ntarget")]:
    ax.axvline(m_, color=style.NEUTRAL, lw=0.8, ls=":")
    ax.annotate(lab, (m_ + 0.04, 0.25), fontsize=7.5, color=style.INK_2)
ax.set_xlabel("exposure ∝ σ^(−m):  m")
ax.set_ylabel("Sharpe / best achievable Sharpe")
ax.set_title(f"Lines: exp(−s²(m+p−2)²/2), s = {s_emp:.2f};  dots: simulated")
ax.legend(title="expected return ∝ σ^p", fontsize=8, title_fontsize=8, loc="lower right")
ax.set_ylim(0.2, 1.05)
style.save(fig, "09_vol_targeting.png")

# ------------------------------------------------------- realistic forecasts
# forecast sigma by an EWMA of squared returns (RiskMetrics-style, lambda 0.94)
p = 0.0
mu = 0.02 * sig.mean() * np.ones(N)
r = mu + sig * eps
ew = lfilter([0.06], [1, -0.94], r**2)
sig_hat = np.sqrt(np.concatenate([[ew[0]], ew[:-1]]))  # known at t-1
ls_hat = np.log(sig_hat[1000:])
ls = logsig[1000:]
b_est, a_est = np.polyfit(ls, ls_hat, 1)
e_est = np.std(ls_hat - (a_est + b_est * ls))
log(f"EWMA forecast: log sigma_hat = {a_est:.2f} + {b_est:.3f} log sigma + noise(sd {e_est:.3f})")
r_ = r[1000:]
sr_opt = sharpe(mu[1000:] / sig[1000:] ** 2, r_)
sims = np.array([sharpe(sig_hat[1000:] ** (-m), r_) for m in ms])
m_star = b_est * s_emp**2 * (2 - p) / (b_est**2 * s_emp**2 + e_est**2)
log(f"with EWMA forecast, p=0: best m simulated {ms[np.argmax(sims)]:.2f}, "
    f"theory m* = b s^2 (2-p) / (b^2 s^2 + e^2) = {m_star:.2f}")
log(f"   SR at m=0 {sims[np.argmin(np.abs(ms))]:.3f}, m=1 {sims[np.argmin(np.abs(ms-1))]:.3f}, "
    f"m=2 {sims[np.argmin(np.abs(ms-2))]:.3f}, optimal-with-true-vol {sr_opt:.3f}")
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(ms, ratio(ms, p, s_emp, b_est, e_est), color=style.SERIES[0],
        label="theory with forecast error")
ax.plot(ms, sims / sr_opt, color=style.SERIES[0], lw=0, marker="o", ms=3.5,
        label="simulated, EWMA(0.94) forecast")
ax.plot(ms, ratio(ms, p, s_emp), color=style.NEUTRAL, ls="--", lw=1.2,
        label="perfect forecast")
ax.axvline(m_star, color=style.SERIES[1], lw=1)
ax.annotate(f"m* = {m_star:.2f}", (m_star + 0.05, 1.02), fontsize=8, color=style.INK)
ax.set_xlabel("exposure ∝ σ̂^(−m):  m")
ax.set_ylabel("Sharpe / best achievable Sharpe")
ax.set_title("Constant expected return (p = 0): noise shrinks the optimal m")
ax.legend(fontsize=8, loc="lower center")
ax.set_ylim(0.2, 1.05)
style.save(fig, "09_forecast_error.png")

open(style.os.path.join(style.FIG_DIR, "09_output.txt"), "w").write("\n".join(out) + "\n")
