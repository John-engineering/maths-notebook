"""Note 20: square-root impact from a diffusing latent order book.

Latent buy (A) and sell (B) intentions diffuse in price space with coefficient D
and annihilate on contact. The net density phi = rho_B - rho_A then solves the
heat equation exactly; in the stationary state phi(x) = L x (linear book). A
metaorder buying at rate m for time T removes sell liquidity at the current
price p_t, so

    L p_t = m int_0^{min(t,T)} K(p_t - p_s, t - s) ds,   K(y, u) = exp(-y^2/4Du)/sqrt(4 pi D u).

Solved step by step with exact product integration of the heat kernel on each
interval: int u^-1/2 e^{-b/u} du = 2 sqrt(u) e^{-b/u} - 2 sqrt(pi b) erfc(sqrt(b/u)).
Units: D = L = 1, so the natural flux J = D L = 1 and r = m / J.
Regimes: r << 1: I = (m/L) sqrt(T/(pi D));  r >> 1: I = sqrt(2 Q / L), Q = m T.
"""
import numpy as np
from scipy.special import erfc
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


D = 1.0
L = 1.0


def F(u, b):
    u = np.maximum(u, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        term1 = 2 * np.sqrt(u) * np.exp(-np.where(u > 0, b / np.where(u > 0, u, 1), np.inf))
        term2 = 2 * np.sqrt(np.pi * b) * erfc(np.sqrt(np.where(u > 0, b / np.where(u > 0, u, 1), np.inf)))
    return np.where(u > 0, term1 - term2, 0.0)


def solve(m, T=1.0, t_end=10.0, n_exec=400):
    dt = T / n_exec
    n = int(round(t_end / dt))
    t = np.arange(n + 1) * dt
    p = np.zeros(n + 1)
    active = t[:-1] < T - 1e-12                    # intervals during which the metaorder trades
    for k in range(1, n + 1):
        tk = t[k]
        lo = t[:k]
        hi = t[1:k + 1]
        act = active[:k]
        pmid_prev = 0.5 * (p[:k - 1] + p[1:k])     # midpoints of completed intervals

        def resid(x):
            pm = np.concatenate([pmid_prev, [0.5 * (p[k - 1] + x)]])
            b = (x - pm) ** 2 / (4 * D)
            u1 = tk - hi
            u2 = tk - lo
            integ = (F(u2, b) - F(u1, b)) / np.sqrt(4 * np.pi * D)
            return L * x - m * np.sum(integ[act])

        if not act.any() and k > 1:
            pass
        a, bnd = p[k - 1] - 5.0, p[k - 1] + 5.0
        ra, rb = resid(a), resid(bnd)
        while ra > 0:
            a -= 5.0
            ra = resid(a)
        while rb < 0:
            bnd += 5.0
            rb = resid(bnd)
        p[k] = brentq(resid, a, bnd, xtol=1e-10)
    return t, p


rs = np.geomspace(0.01, 100, 13)
peaks = []
paths = {}
for r in rs:
    t, p = solve(r, T=1.0, t_end=1.0 if r not in (0.1, 1.0, 10.0) else 10.0)
    peak = p[np.argmin(np.abs(t - 1.0))]
    peaks.append(peak)
    log(f"r = m/J = {r:8.3f}: peak impact I/sqrt(DT) = {peak:.4f};  slow-regime r/sqrt(pi) = "
        f"{r/np.sqrt(np.pi):.4f};  fast-regime sqrt(2r) = {np.sqrt(2*r):.4f}")
peaks = np.array(peaks)
for r in (0.1, 1.0, 10.0):
    t, p = solve(r, T=1.0, t_end=10.0)
    paths[r] = (t, p)

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
axes[0].loglog(rs, peaks, marker="o", ms=5, lw=0, color=style.SERIES[0], label="numerical solution")
rr = np.geomspace(0.005, 200, 100)
axes[0].loglog(rr, rr / np.sqrt(np.pi), color=style.SERIES[1], ls="--", lw=1.2, label="slow: (m/L)√(T/πD), linear in Q")
axes[0].loglog(rr, np.sqrt(2 * rr), color=style.SERIES[2], ls="--", lw=1.2, label="fast: √(2Q/L), square root")
axes[0].set_xlabel("execution rate relative to latent liquidity flux, m / (D L)")
axes[0].set_ylabel("peak impact / √(DT)")
axes[0].set_title("From linear to square-root impact")
axes[0].legend(fontsize=8)
style.plain_log(axes[0])
style.plain_log(axes[0], "x")
for i, r in enumerate((0.1, 1.0, 10.0)):
    t, p = paths[r]
    pk = p[np.argmin(np.abs(t - 1.0))]
    axes[1].plot(t, p / pk, color=style.SERIES[i], label=f"m/J = {r:g}")
tt = np.linspace(1.0001, 10, 300)
axes[1].plot(tt, np.sqrt(tt) - np.sqrt(tt - 1), color=style.NEUTRAL, ls="--", lw=1.2,
             label="√t − √(t−1) (slow-regime decay)")
axes[1].axvline(1.0, color=style.NEUTRAL, lw=0.8, ls=":")
axes[1].set_xlabel("time / execution duration T")
axes[1].set_ylabel("impact / peak impact")
axes[1].set_title("During and after the metaorder")
axes[1].legend(fontsize=8)
style.save(fig, "20_latent_book.png")

# decay after completion: impact at t = 2T and 10T relative to the peak
for r in (0.1, 1.0, 10.0):
    t, p = paths[r]
    pk = p[np.argmin(np.abs(t - 1.0))]
    v2 = p[np.argmin(np.abs(t - 2.0))] / pk
    v10 = p[np.argmin(np.abs(t - 10.0))] / pk
    log(f"relaxation m/J={r:g}: I(2T)/I(T) = {v2:.3f}, I(10T)/I(T) = {v10:.3f}; slow-regime formula "
        f"{np.sqrt(2)-1:.3f}, {np.sqrt(10)-np.sqrt(9):.3f}")

open(style.os.path.join(style.FIG_DIR, "20_output.txt"), "w").write("\n".join(out) + "\n")
