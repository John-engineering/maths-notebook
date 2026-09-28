"""Note 10: market impact kernels, efficiency and Bochner's theorem.

Propagator model: price = sum_s G(t-s) * (signed volume at s) + noise.
A. Round-trip cost is (1/2) v' G v; no manipulation <=> G positive definite
   <=> Fourier transform of G(|t|) >= 0 (Bochner).
B. Optimal buy schedules v* ~ G^{-1} 1 for different kernels: positive-definite
   but non-convex kernels make the optimal *buy* program sell at times.
C. With long-memory order signs (autocorrelation ~ l^-gamma), prices are
   diffusive only if G(l) ~ l^-beta with beta = (1 - gamma)/2.
"""
import numpy as np
from scipy.linalg import toeplitz
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(10)
out = []


def log(msg):
    print(msg)
    out.append(msg)


n = 240
t = np.arange(n, dtype=float)
TAU = 20.0
kernels = {
    "exponential e^(−t/τ)": np.exp(-t / TAU),
    "power law (1+t/τ)^(−0.5)": (1 + t / TAU) ** -0.5,
    "triangle (1−t/2τ)⁺": np.clip(1 - t / (2 * TAU), 0, None),
    "Gaussian e^(−t²/2τ²)": np.exp(-t**2 / (2 * TAU**2)),
    "rectangle 1{t<τ}": (t < TAU).astype(float),
}

# ------------------------------------------------------------------ A. Bochner
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
fig.subplots_adjust(wspace=0.25)
min_eigs = {}
for i, (name, g) in enumerate(kernels.items()):
    Gm = toeplitz(g)
    w, V = np.linalg.eigh(Gm)
    min_eigs[name] = w[0]
    # symbol of the Toeplitz operator (discrete-time Fourier transform of
    # G(|k|)), from a long sequence so truncation ripples are negligible
    K = 1 << 18
    kk = np.arange(K, dtype=float)
    gl = {"exponential e^(−t/τ)": np.exp(-kk / TAU),
          "power law (1+t/τ)^(−0.5)": (1 + kk / TAU) ** -0.5,
          "triangle (1−t/2τ)⁺": np.clip(1 - kk / (2 * TAU), 0, None),
          "Gaussian e^(−t²/2τ²)": np.exp(-kk**2 / (2 * TAU**2)),
          "rectangle 1{t<τ}": (kk < TAU).astype(float)}[name]
    sym = 2 * np.fft.rfft(gl).real - gl[0]
    om = 2 * np.pi * np.arange(len(sym)) / K
    sel = (om > 0.02 / TAU) & (om < 1.0)
    ref = sym[np.argmin(np.abs(om - 0.02 / TAU))]
    ft = sym[sel]
    axes[0].plot(om[sel] * TAU, ft / ref, color=style.SERIES[i], label=name)
    log(f"A  {name:28s} min eigenvalue of Toeplitz G: {w[0]:+.2e}   min of symbol: {sym.min():+.3e}")
axes[0].axhline(0, color=style.NEUTRAL, lw=1)
axes[0].set_xlabel("frequency ω·τ")
axes[0].set_ylabel("Ĝ(ω), normalised at low frequency")
axes[0].set_xscale("log")
axes[0].set_ylim(-0.4, 1.05)
axes[0].set_title("Bochner: manipulation-free ⇔ Ĝ ≥ 0")
axes[0].legend(fontsize=7.5)
# the most profitable round trip under the rectangular kernel
g = kernels["rectangle 1{t<τ}"]
w, V = np.linalg.eigh(toeplitz(g))
v = V[:, 0]
v = v / np.abs(v).sum() * 2                   # total traded volume 2 (1 bought, 1 sold)
cost = 0.5 * v @ toeplitz(g) @ v
axes[1].bar(t, v, width=1.0, color=np.where(v > 0, style.SERIES[0], style.SERIES[1]))
axes[1].axhline(0, color=style.NEUTRAL, lw=1)
axes[1].set_xlabel("time step")
axes[1].set_ylabel("signed trade size (blue buy, orange sell)")
axes[1].set_title(f"Rectangular decay: a round trip that earns {-cost:.3f}·η per unit² traded")
style.save(fig, "10_bochner.png")
# the period of the manipulation equals 2 pi / omega at the most negative FT
sign_changes = np.sum(np.diff(np.sign(v[np.abs(v) > 1e-6])) != 0)
log(f"A  rectangle: most profitable round trip has {sign_changes} sign changes over {n} steps "
    f"(period ~ {2*n/max(sign_changes,1):.1f} steps; theory 2π/(4.493/τ) = {2*np.pi*TAU/4.493:.1f}); "
    f"expected cost {cost:+.4f} (negative = profit)")

# ------------------------------------------------------------------ B. optimal buys
# add a small instantaneous (temporary) impact lam * v^2, e.g. half-spread and
# immediate book depletion; without it the Gaussian problem is ill-posed
LAM = 0.05
fig, axes = plt.subplots(1, 4, figsize=(15, 3.6), sharey=False)
fig.subplots_adjust(wspace=0.3)
m = 100                                          # execution horizon in steps
for ax, (i, (name, g)) in zip(axes, enumerate(list(kernels.items())[:4])):
    Gm = toeplitz(g[:m]) + LAM * np.eye(m)
    x = np.linalg.solve(Gm, np.ones(m))
    v = x / x.sum()                              # buy 1 unit in total
    ax.bar(np.arange(m), v, width=1.0, color=np.where(v >= 0, style.SERIES[0], style.SERIES[1]))
    ax.axhline(0, color=style.NEUTRAL, lw=1)
    ax.set_title(name, fontsize=10)
    ax.set_xlabel("time step")
    sells = -v[v < 0].sum()
    log(f"B  {name:28s} optimal buy program: first/last trade {v[0]:.3f}/{v[-1]:.3f}, "
        f"interior mean {v[5:-5].mean():.4f}, total sold along the way {sells:.3f}, "
        f"cost {0.5*v@Gm@v:.4f} vs uniform {0.5*np.ones(m)@Gm@np.ones(m)/m**2:.4f}")
axes[0].set_ylabel("trade size (orange = sells)")
style.save(fig, "10_optimal_schedules.png")

# ------------------------------------------------------------------ C. efficiency
def fgn(N, H, r_):
    """Fractional Gaussian noise by circulant embedding (Davies-Harte)."""
    k = np.arange(N + 1)
    gam = 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))
    row = np.concatenate([gam, gam[-2:0:-1]])
    lam = np.fft.fft(row).real
    lam = np.clip(lam, 0, None)
    M = len(row)
    Z = r_.standard_normal(M) + 1j * r_.standard_normal(M)
    X = np.fft.fft(np.sqrt(lam / M) * Z)
    return X.real[:N]


N = 2**21
H_signs = 0.75
eps = np.sign(fgn(N, H_signs, rng))           # order signs with long memory
gamma = 2 - 2 * H_signs
lags_ac = np.unique(np.geomspace(1, 2000, 25).astype(int))
e0 = eps - eps.mean()
ac = np.array([np.mean(e0[l:] * e0[:-l]) for l in lags_ac]) / np.var(eps)
gamma_fit = -np.polyfit(np.log(lags_ac[5:]), np.log(ac[5:]), 1)[0]
log(f"C  order signs: target gamma {gamma:.2f}, fitted autocorrelation exponent {gamma_fit:.3f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
fig.subplots_adjust(wspace=0.25)
axes[0].loglog(lags_ac, ac, marker="o", ms=4, color=style.SERIES[0], lw=0, label="order-sign autocorrelation")
axes[0].loglog(lags_ac, ac[3] * (lags_ac / lags_ac[3]) ** (-gamma), color=style.NEUTRAL, ls="--",
               label=f"ℓ^−{gamma:.1f}")
axes[0].set_xlabel("lag ℓ (trades)")
axes[0].set_ylabel("C(ℓ)")
axes[0].set_title("Order flow is long-memory (metaorders)")
axes[0].legend(fontsize=8)
style.plain_log(axes[0])
style.plain_log(axes[0], "x")
Lags = np.unique(np.geomspace(1, 20000, 30).astype(int))
beta_star = (1 - gamma) / 2
for i, beta in enumerate([0.0, 0.1, beta_star, 0.4, 0.5]):
    ell = np.arange(1, N + 1, dtype=float)
    Gk = ell ** (-beta)
    L2 = 1 << int(np.ceil(np.log2(2 * N)))
    p = np.fft.irfft(np.fft.rfft(eps, L2) * np.fft.rfft(Gk, L2), L2)[:N]
    p = p[N // 4:]                                # discard start-up
    D = np.array([np.var(p[l:] - p[:-l]) for l in Lags])
    slope = np.polyfit(np.log(Lags[8:]), np.log(D[8:]), 1)[0]
    lab = f"β = {beta:.2f}" + ("  = (1−γ)/2" if abs(beta - beta_star) < 1e-9 else "")
    axes[1].loglog(Lags, D / Lags / (D[0] / Lags[0]), color=style.SERIES[i], label=f"{lab}: D ∝ L^{slope:.2f}")
    log(f"C  beta={beta:.3f}: variogram exponent {slope:.3f} (theory 2 - 2beta - gamma = {2-2*beta-gamma:.3f}; "
        f"diffusive = 1)")
axes[1].set_xlabel("lag L (trades)")
axes[1].set_ylabel("D(L) / L  (normalised)")
axes[1].set_title("Only β = (1−γ)/2 gives a diffusive price")
axes[1].legend(fontsize=8)
style.plain_log(axes[1])
style.plain_log(axes[1], "x")
style.save(fig, "10_efficiency.png")

open(style.os.path.join(style.FIG_DIR, "10_output.txt"), "w").write("\n".join(out) + "\n")
