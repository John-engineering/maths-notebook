"""Note 04: trend following is the discrete Ito formula.

1. Exact identities: position = cumulative move -> P&L = 1/2[(total)^2 - sum r^2];
   EMA position -> P&L = sum [(1-l^2) m_t^2 - r_t^2] / (2 l) + boundary terms.
2. Spectral filter H(w) = (cos w - l) / (1 - 2 l cos w + l^2); integrates to 0,
   crossover period ~ pi sqrt(2N).
3. Expected P&L = sum_k l^k gamma(k+1) = integral H S; checked on four
   autocorrelation structures.
4. Drift-plus-noise world: optimal EMA decay l* = rho (persistence of drift),
   Sharpe ~ SNR * rho / sqrt(1 - rho^2).
5. Convexity: P&L vs total move is a parabola (long straddle).
6. Mirror: rebalancer's relative log performance ~ -1/4 x (trend P&L on log ratio).
"""
import numpy as np
from scipy.signal import lfilter
from scipy.integrate import quad
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(4)
out = []


def log(msg):
    print(msg)
    out.append(msg)


def ema_signal(r, lam):
    """m_t = lam m_{t-1} + r_t (position held over the NEXT step)."""
    return lfilter([1.0], [1.0, -lam], r)


def pnl(r, m):
    return np.sum(m[:-1] * r[1:])


# ------------------------------------------------------------------ 1. identity
r = 0.01 * rng.standard_t(4, 5000)
S = np.concatenate([[0], np.cumsum(r)])
theta = S[:-1]                                   # position = S_t - S_0
lhs = np.sum(theta * r)
rhs = 0.5 * (S[-1] ** 2 - np.sum(r**2))
log(f"1a cumulative-position identity: {lhs:.12f} vs {rhs:.12f}")
lam = 0.95
m = ema_signal(r, lam)
lhs = pnl(r, m)
T = len(r)
rhs = ((1 - lam**2) * np.sum(m[:-1] ** 2) + m[-1] ** 2 - lam**2 * m[0] ** 2
       - np.sum(r[1:] ** 2) - (1 - lam**2) * 0) / (2 * lam)
# exact bookkeeping: sum_{t=0}^{T-2} m_t r_{t+1}
#   = [m_{T-1}^2 - lam^2 m_0^2 + (1-lam^2) sum_{t=1}^{T-2} m_t^2 - sum_{t=1}^{T-1} r_t^2]/(2 lam)
rhs = (m[-1] ** 2 - lam**2 * m[0] ** 2 + (1 - lam**2) * np.sum(m[1:-1] ** 2)
       - np.sum(r[1:] ** 2)) / (2 * lam)
log(f"1b EMA identity (lam={lam}):          {lhs:.12f} vs {rhs:.12f}")

# ------------------------------------------------------------------ 2. filters
w = np.linspace(1e-4, np.pi, 4000)


def H(w, lam):
    return (np.cos(w) - lam) / (1 - 2 * lam * np.cos(w) + lam**2)


fig, ax = plt.subplots(figsize=(7.5, 4.2))
for i, N in enumerate([10, 50, 250]):
    lam_ = 1 - 1 / N
    freq = w / (2 * np.pi)
    ax.plot(freq, H(w, lam_), color=style.SERIES[i], label=f"EMA N = {N}")
    pc = 2 * np.pi / np.arccos(lam_)
    ax.axvline(1 / pc, color=style.SERIES[i], lw=0.8, ls=":")
    log(f"2  N={N:4d}: crossover period {pc:6.1f} steps; pi*sqrt(2N) = {np.pi*np.sqrt(2*N):6.1f}; "
        f"integral of H = {quad(H, 0, np.pi, args=(lam_,), limit=500, points=[np.arccos(lam_)])[0]:+.1e}")
ax.axhline(0, color=style.NEUTRAL, lw=1)
ax.set_xscale("log")
ax.set_yscale("symlog", linthresh=1.0)
ax.set_xlabel("frequency (cycles per step, log); dotted: crossover 1/(π√2N)")
ax.set_ylabel("H(ω): P&L weight on spectral power")
ax.set_title("Long slow cycles, short fast ones (∫H dω = 0)")
ax.legend(loc="upper right")
style.save(fig, "04_filters.png")

# ------------------------------------------------------------------ 3. worlds
def arma_returns(n, kind, r_):
    e = r_.standard_normal(n + 2000)
    if kind == "iid":
        x = e
    elif kind == "momentum AR(1) +0.05":
        x = lfilter([1], [1, -0.05], e)
    elif kind == "reversal AR(1) -0.05":
        x = lfilter([1], [1, 0.05], e)
    elif kind == "slow drift + fast reversal":
        mu = lfilter([1], [1, -0.99], 0.02 * r_.standard_normal(n + 2000))
        noise = e + (-0.15) * np.concatenate([[0], e[:-1]])     # MA(1) reversal
        x = mu + noise
    return x[2000:]


def autocov(kind, K=4000):
    """Theoretical autocovariances gamma(1..K) for each world."""
    k = np.arange(1, K + 1)
    if kind == "iid":
        return np.zeros(K)
    if kind == "momentum AR(1) +0.05":
        phi = 0.05
        return phi**k / (1 - phi**2)
    if kind == "reversal AR(1) -0.05":
        phi = -0.05
        return phi**k / (1 - phi**2)
    if kind == "slow drift + fast reversal":
        rho, s = 0.99, 0.02
        g = s**2 / (1 - rho**2) * rho**k
        g[0] += -0.15
        return g


kinds = ["iid", "momentum AR(1) +0.05", "reversal AR(1) -0.05", "slow drift + fast reversal"]
Ns = np.array([2, 3, 5, 8, 12, 20, 35, 60, 100, 170, 300])
n_sim = 2_000_000
res = {}
for kind in kinds:
    x = arma_returns(n_sim, kind, np.random.default_rng(7))
    sims, theo = [], []
    g = autocov(kind)
    for N in Ns:
        lam_ = 1 - 1 / N
        m_ = ema_signal(x, lam_)
        # normalise: position scaled so the signal has unit variance
        scale = np.std(m_)
        sims.append(pnl(x, m_) / (len(x) - 1) / scale)
        k = np.arange(len(g))
        theo.append(np.sum(lam_**k * g) / scale)
    res[kind] = (np.array(sims), np.array(theo))
    log(f"3  {kind:28s} sim/theory at N=20: {sims[5]:+.5f} / {theo[5]:+.5f}")

fig, ax = plt.subplots(figsize=(7.5, 4.2))
for i, kind in enumerate(kinds):
    sims, theo = res[kind]
    ax.plot(Ns, theo, color=style.SERIES[i], label=kind)
    ax.plot(Ns, sims, color=style.SERIES[i], marker="o", ms=4, lw=0)
ax.axhline(0, color=style.NEUTRAL, lw=1)
ax.set_xscale("log")
ax.set_xlabel("EMA length N (steps, log)")
ax.set_ylabel("expected P&L per step (unit-vol position)")
ax.set_title("E[P&L] = Σ λᵏ γ(k+1): lines theory, dots simulated")
ax.legend(fontsize=8)
style.plain_log(ax, "x")
style.save(fig, "04_worlds.png")

# ------------------------------------------------------------------ 4. optimal lambda
# r_t = mu_t + eps_t, mu AR(1) with persistence rho. Sharpe(lam) prop. to
# sqrt(1-lam^2)/(1-lam rho): maximised at lam = rho.
for tau in [20, 60, 250]:
    rho = 1 - 1 / tau
    lams = 1 - 1 / np.geomspace(2, 3000, 60)
    f = np.sqrt(1 - lams**2) / (1 - lams * rho)
    log(f"4  drift persistence tau={tau}: argmax lambda -> N = {1/(1-lams[np.argmax(f)]):.0f} "
        f"(theory N = tau = {tau})")
# simulated check for tau = 60
tau = 60
rho = 1 - 1 / tau
snr_mu = 0.05
nn = 3_000_000
mu = lfilter([1], [1, -rho], snr_mu * np.sqrt(1 - rho**2) * rng.standard_normal(nn))
x = mu + rng.standard_normal(nn)
sr = []
for N in Ns:
    lam_ = 1 - 1 / N
    m_ = ema_signal(x, lam_)
    p = m_[:-1] * x[1:]
    sr.append(p.mean() / p.std())
sr = np.array(sr)
theory_sr = (snr_mu**2 * rho * np.sqrt(1 - (1 - 1 / Ns) ** 2)
             / (1 - (1 - 1 / Ns) * rho) / 1.0)
fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(Ns, sr * np.sqrt(252), marker="o", ms=4, lw=0, color=style.SERIES[0], label="simulated")
ax.plot(Ns, theory_sr * np.sqrt(252), color=style.SERIES[0], label="theory")
ax.axvline(tau, color=style.SERIES[1], lw=1, ls=":")
ax.annotate(f"drift persistence τ = {tau}", (tau, 0.02), xytext=(6, 0),
            textcoords="offset points", fontsize=8, color=style.INK_2)
ax.set_xscale("log")
ax.set_xlabel("EMA length N (days, log)")
ax.set_ylabel("annualised Sharpe")
ax.set_title("Best trend filter matches the drift's persistence")
ax.legend()
style.plain_log(ax, "x")
style.save(fig, "04_optimal_lambda.png")
log(f"4  tau=60 sim: best N = {Ns[np.argmax(sr)]}, max annual SR {np.sqrt(252)*sr.max():.3f}, "
    f"theory max {np.sqrt(252)*snr_mu**2*rho/np.sqrt(1-rho**2):.3f}")

# ------------------------------------------------------------------ 5. convexity
paths, steps = 4000, 250
r5 = 0.01 * rng.standard_normal((paths, steps))
total = r5.sum(1)
pn = np.zeros(paths)
for j in range(paths):
    m_ = ema_signal(r5[j], 1 - 1 / 50)
    pn[j] = pnl(r5[j], m_)
fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.scatter(total, pn, s=6, alpha=0.35, color=style.SERIES[0], edgecolors="none")
ax.axhline(0, color=style.NEUTRAL, lw=1)
ax.set_xlabel("total move over the year")
ax.set_ylabel("trend-follower P&L (EMA 50)")
ax.set_title("Zero-drift random walk: P&L is a straddle on the move")
style.save(fig, "04_convexity.png")
log(f"5  convexity: corr(P&L, total^2) = {np.corrcoef(pn, total**2)[0,1]:.3f}, mean P&L {pn.mean():+.2e}")

# ------------------------------------------------------------------ 6. mirror
# two assets; log ratio L. Rebalanced vs buy-and-hold relative log performance
# versus (1/4) x trend P&L with position L_t - L_0.
T6 = 2000
dL = 0.02 * rng.standard_normal(T6)
L = np.concatenate([[0], np.cumsum(dL)])
X1, X2 = L / 2, -L / 2                           # log prices, common part removed
R1, R2 = np.exp(np.diff(X1)), np.exp(np.diff(X2))
log_rb = np.sum(np.log(0.5 * R1 + 0.5 * R2))
log_bh = np.log(0.5 * np.exp(X1[-1]) + 0.5 * np.exp(X2[-1]))
trend = np.sum(L[:-1] * dL)
log(f"6  rebalanced - buy&hold = {log_rb - log_bh:+.4f};  -(1/4) trend P&L on L = {-trend/4:+.4f}")

open(style.os.path.join(style.FIG_DIR, "04_output.txt"), "w").write("\n".join(out) + "\n")
