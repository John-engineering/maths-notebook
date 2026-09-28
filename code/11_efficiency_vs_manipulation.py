"""Note 11: does an efficient price imply no manipulation?

Efficient propagator: if order signs have best linear predictor
eps_hat_t = sum_j a_j eps_{t-j}, the kernel that makes the price a martingale is
G(k) = G(1) (1 - sum_{j<k} a_j), i.e. price = G(1) x cumulative surprise.
A strategic trader's expected cost is (1/2) v' T v with T(0) = G(1)(1 + s)
(s = extra immediate cost, e.g. half-spread, in units of G(1)), T(k) = G(k).

Symbol: T_hat(w) = G(1) [ s + sum_j a_j D_j(w) ],  D_j = Dirichlet kernel.
At w = pi: s - (a_1 - a_2 + a_3 - ...). So without a spread, persistent order
flow makes the efficient price manipulable by alternating trades.
"""
import numpy as np
from scipy.linalg import toeplitz, solve_toeplitz
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg)
    out.append(msg)


K = 3000


def fgn_sign_acf(H, K):
    k = np.arange(K + 1, dtype=float)
    rho = 0.5 * (np.abs(k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))
    return (2 / np.pi) * np.arcsin(rho)


def mixture_acf(K, phis=(0.3, 0.8, 0.97, 0.995), ws=(0.4, 0.25, 0.2, 0.15)):
    k = np.arange(K + 1, dtype=float)
    return sum(w * p**k for p, w in zip(phis, ws))


def predictor(acf, K):
    """Best linear predictor coefficients a_1..a_K (Yule-Walker)."""
    return solve_toeplitz(acf[:K], acf[1:K + 1])


def symbol_minus_s(a, om):
    j = np.arange(1, len(a) + 1)
    re = np.cos(np.outer(om, j)) @ a
    im = np.sin(np.outer(om, j)) @ a
    return re + im / np.tan(om / 2)


cases = {
    "fGn signs, H = 0.60 (γ = 0.8)": fgn_sign_acf(0.60, K + 1),
    "fGn signs, H = 0.75 (γ = 0.5)": fgn_sign_acf(0.75, K + 1),
    "fGn signs, H = 0.90 (γ = 0.2)": fgn_sign_acf(0.90, K + 1),
    "mixture of metaorder AR(1)s": mixture_acf(K + 1),
    "AR(1), φ = 0.5": 0.5 ** np.arange(K + 2, dtype=float),
}
om = np.linspace(0.01, np.pi, 1500)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
fig.subplots_adjust(wspace=0.25)
rows = []
for i, (name, acf) in enumerate(cases.items()):
    a = predictor(acf, K)
    G = 1 - np.concatenate([[0.0], np.cumsum(a)])          # G(1..K+1), G(1)=1
    ss = symbol_minus_s(a, om)
    s_star = max(0.0, -ss.min())
    w_min = om[np.argmin(ss)]
    alt = np.sum(a * (-1.0) ** (np.arange(1, K + 1) + 1))    # a1 - a2 + a3 - ...
    mono = np.all(np.diff(a[:200]) <= 1e-12)
    pos = np.all(a[:200] > 0)
    # check with a finite Toeplitz matrix just above and below the threshold
    n = 400
    def min_eig(s):
        col = np.concatenate([[1 + s], G[:n - 1]])
        return np.linalg.eigvalsh(toeplitz(col))[0]
    e_lo, e_hi = min_eig(s_star - 0.02), min_eig(s_star + 0.02)
    k = np.arange(1, K + 2)
    sel = (k > 30) & (k < 1000)
    # for long-memory flow the kernel is fully transient (G -> 0) and decays
    # like k^-(1-gamma)/2; fit log G directly
    beta_fit = -np.polyfit(np.log(k[sel]), np.log(np.clip(G[sel], 1e-15, None)), 1)[0]
    sig_innov = np.sqrt(1 - np.sum(a * acf[1:K + 1]))           # prediction error sd
    rows.append((name, a[0], alt, s_star, w_min, e_lo, e_hi))
    log(f"   -> innovation sd {sig_innov:.3f}; no-manipulation bound: "
        f"half-spread / volatility-per-trade >= s*/sigma_innov = {s_star/sig_innov:.3f}")
    log(f"{name:34s} a1={a[0]:.3f}  a1-a2+a3-...={alt:.3f}  min spread s*={s_star:.3f} "
        f"at ω={w_min:.3f} (π={np.pi:.3f}); a_j>0: {pos}, decreasing: {mono}; "
        f"Toeplitz min eig at s*∓0.02: {e_lo:+.4f} / {e_hi:+.4f}; kernel log-slope on [30,1000] {beta_fit:.2f}")
    show = k <= 1000                      # beyond ~K/3 the finite predictor order bites
    axes[0].loglog(k[show], G[show], color=style.SERIES[i], label=name)
    axes[1].plot(om, ss, color=style.SERIES[i], label=name)
axes[0].set_xlabel("lag k (trades)")
axes[0].set_ylabel("efficient impact kernel G(k) / G(1)")
axes[0].set_title("Martingale kernels: G(k) = 1 − Σ_{j<k} a_j")
axes[0].legend(fontsize=7.5, loc="lower left")
style.plain_log(axes[0])
style.plain_log(axes[0], "x")
axes[1].axhline(0, color=style.NEUTRAL, lw=1)
axes[1].set_xlabel("frequency ω")
axes[1].set_ylabel("Σ a_j D_j(ω)  (symbol without spread)")
axes[1].set_title("Negative at high frequency ⇒ need spread s ≥ −min")
axes[1].set_ylim(-1.0, 3.0)
style.save(fig, "11_efficient_kernels.png")

# ----------------------------------------------------------- a worked round trip
acf = cases["fGn signs, H = 0.75 (γ = 0.5)"]
a = predictor(acf, K)
n = 60
v = np.array([(-1) ** t for t in range(n)], float)
v[-1] = -v[:-1].sum()                                      # close out: sum v = 0
flow_pred = np.array([np.sum(a[:t][::-1] * v[:t]) for t in range(n)])
surprise = v - flow_pred
price = np.concatenate([[0.0], np.cumsum(surprise)])[:-1]    # price seen by trade t
for s in [0.0, 0.2, 0.5]:
    cost = np.sum(v * price) + 0.5 * (1 + s) * np.sum(v**2)
    log(f"alternating round trip ({n} trades), spread s={s:.1f}: expected cost {cost:+.3f} "
        f"({'profit' if cost < 0 else 'loss'} for the manipulator)")
fig, ax = plt.subplots(figsize=(8, 3.8))
ax.step(np.arange(n), price, where="post", color=style.SERIES[0], label="efficient price (cumulative surprise)")
ax.bar(np.arange(n), 0.25 * v, width=0.6, color=np.where(v > 0, style.SERIES[2], style.SERIES[1]),
       label="trades (scaled; green buy, orange sell)")
ax.axhline(0, color=style.NEUTRAL, lw=1)
ax.set_xlabel("trade number")
ax.set_ylabel("price (units of G(1))")
ax.set_title("A market that expects persistence overreacts to reversals")
ax.legend(fontsize=8, loc="upper right")
style.save(fig, "11_round_trip.png")

open(style.os.path.join(style.FIG_DIR, "11_output.txt"), "w").write("\n".join(out) + "\n")
