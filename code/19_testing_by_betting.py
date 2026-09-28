"""Note 19: the Kelly bettor is a statistician.

1. Kelly wealth (leverage f = mu/sigma^2 on the asset) equals the likelihood
   ratio dP_mu/dP_0 of the observed path, path by path.
2. A Gaussian mixture of Kelly bettors (over leverage) has wealth equal to
   Robbins' normal-mixture martingale; "reject no-edge when wealth >= 1/alpha"
   is exactly note 07's always-valid boundary.
3. CUSUM (note 18) = log-wealth of a bettor wagering that the edge died,
   restarted at 1 whenever it would fall below 1.
"""
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(19)
out = []


def log(msg):
    print(msg)
    out.append(msg)


# monthly excess returns in units of monthly vol, true Sharpe 0.8/yr
sr = 0.8
mu = sr / np.sqrt(12)
T = 360
x = mu + rng.standard_normal(T)
S = np.cumsum(x)
t = np.arange(1, T + 1)

# 1. Kelly wealth (continuous-time idealisation per step) vs likelihood ratio
f = mu                                                   # Kelly leverage in vol units
logW_kelly = f * S - 0.5 * f**2 * t
logLR = np.cumsum(norm.logpdf(x, mu, 1) - norm.logpdf(x, 0, 1))
log(f"1  max |log Kelly wealth - log likelihood ratio| over {T} months: "
    f"{np.max(np.abs(logW_kelly - logLR)):.2e}")

# 2. mixture of Kelly bettors over leverage lambda ~ N(0, 1/rho): numerical vs closed form
rho = 12.0
lam = np.linspace(-3, 3, 4001)
wts = norm.pdf(lam, 0, 1 / np.sqrt(rho))
wts /= np.trapezoid(wts, lam)
mix_num = np.array([np.trapezoid(wts * np.exp(lam * S[k] - 0.5 * lam**2 * t[k]), lam) for k in range(T)])
mix_closed = np.sqrt(rho / (t + rho)) * np.exp(S**2 / (2 * (t + rho)))
log(f"2  mixture wealth: numerical integral vs Robbins closed form, max rel. error "
    f"{np.max(np.abs(mix_num / mix_closed - 1)):.2e}")

# rejection when wealth >= 20 (alpha = 5%) under the null, many paths
N0 = 20000
X0 = rng.standard_normal((N0, T))
S0 = np.cumsum(X0, axis=1)
M0 = np.sqrt(rho / (t + rho)) * np.exp(S0**2 / (2 * (t + rho)))
rej = (M0 >= 20).any(axis=1).mean()
log(f"2  zero-skill bettors whose mixture wealth ever reaches 20x in 30y: {rej:.4f} (Ville bound 0.05)")

# 3. CUSUM as the log-wealth of a bettor on "edge has died" (mean 0 vs mu), floored
x2 = np.concatenate([mu + rng.standard_normal(180), rng.standard_normal(180)])   # dies at month 180
l = norm.logpdf(x2, 0, 1) - norm.logpdf(x2, mu, 1)          # log-LR increments "dead" vs "alive"
W = np.zeros(T)
w = 0.0
for k in range(T):
    w = max(0.0, w + l[k])
    W[k] = w
# the same thing as a betting strategy: bet (with Kelly leverage mu) that the
# strategy UNDER-performs its alive expectation, i.e. hold -mu units of (x - mu);
# fair (a martingale) if the edge is alive. Reset wealth to 1 when it falls below 1.
logw_bet = np.zeros(T)
lw = 0.0
for k in range(T):
    lw = max(0.0, lw + (-mu) * (x2[k] - mu) - 0.5 * mu**2)
    logw_bet[k] = lw
log(f"3  max |CUSUM - log-wealth of restarted 'edge died' bettor|: {np.max(np.abs(W - logw_bet)):.2e}")

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))
fig.subplots_adjust(wspace=0.28)
axes[0].plot(t / 12, logW_kelly, color=style.SERIES[0], label="log Kelly wealth")
axes[0].plot(t / 12, logLR, color=style.SERIES[1], ls="--", label="log likelihood ratio (edge vs none)")
axes[0].set_xlabel("years")
axes[0].set_ylabel("log")
axes[0].set_title("Kelly wealth = likelihood ratio")
axes[0].legend(fontsize=8)
axes[1].plot(t / 12, np.log(mix_closed), color=style.SERIES[0], label="mixture of Kelly bettors (true Sharpe 0.8)")
for k in range(30):
    axes[1].plot(t / 12, np.log(M0[k]), color=style.NEUTRAL, lw=0.6, alpha=0.6,
                 label="zero-skill paths" if k == 0 else None)
axes[1].legend(fontsize=8, loc="upper left")
axes[1].axhline(np.log(20), color=style.SERIES[1], lw=1.2)
axes[1].annotate("reject 'no edge' at 20× (α = 5%)", (0.5, np.log(20) + 0.15), fontsize=8, color=style.INK)
axes[1].set_xlabel("years")
axes[1].set_ylabel("log wealth")
axes[1].set_title("Always-valid test = betting")
axes[1].set_ylim(-3, 8)
axes[2].plot(t / 12, W, color=style.SERIES[0], label="CUSUM statistic")
axes[2].axvline(15, color=style.SERIES[1], lw=1, ls=":")
axes[2].annotate("edge dies", (15.2, W.max() * 0.9), fontsize=8, color=style.INK)
axes[2].set_xlabel("years")
axes[2].set_ylabel("log wealth of the 'it died' bettor")
axes[2].set_title("CUSUM = restarted bettor")
style.save(fig, "19_testing_by_betting.png")

open(style.os.path.join(style.FIG_DIR, "19_output.txt"), "w").write("\n".join(out) + "\n")
