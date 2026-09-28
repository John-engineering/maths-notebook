"""Note 16: symmetric spreads vs skewed quotes as defences against manipulation.

Efficient mid: p_t = G(1) * cumulative surprise, eps_hat_t = sum_j a_j eps_{t-j}.
Unit trades v_t in {-1, 0, 1}.

(i) Fixed symmetric half-spread h: trader pays p_t + h sign(v_t).
    Expected cost = (1/2) v' M v, M(0) = 2h, M(k) = G(k)     (k >= 1)
    Break-even (average adverse selection) h* = G(1) sigma_innov^2.
    Manipulation-proof iff h >= G(1) (1 + s*)/2, s* = a1 - a2 + a3 - ...
(ii) Glosten-Milgrom skewed quotes: ask = p_t + G(1)(1 - eps_hat), bid = p_t - G(1)(1 + eps_hat)
    = the post-trade price. Expected cost = (1/2) v' Q v, Q(0) = 2G(1), Q(k) = G(k+1).
    Symbol = G(1)[1 + sum_j a_j D_{j-1}(w)] >= G(1) > 0 when a_j >= 0 decreasing
    (Abel summation into Fejer kernels).
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


def analyse(acf, n=300):
    a = solve_toeplitz(acf[:K], acf[1:K + 1])
    G = 1 - np.concatenate([[0.0], np.cumsum(a)])        # G(1..), G(1)=1
    s_star = max(0.0, np.sum(a * (-1.0) ** np.arange(K)))  # a1 - a2 + ...
    var_innov = 1 - np.sum(a * acf[1:K + 1])
    h_be = var_innov                                     # break-even half-spread (units G(1))
    h_min = (1 + s_star) / 2                             # manipulation-proof half-spread
    M = toeplitz(np.concatenate([[2 * h_be], G[:n - 1]]))
    Q = toeplitz(np.concatenate([[2.0], G[1:n]]))
    return dict(a=a, G=G, s_star=s_star, var_innov=var_innov, h_be=h_be, h_min=h_min,
                eig_fixed=np.linalg.eigvalsh(M)[0], eig_gm=np.linalg.eigvalsh(Q)[0], M=M, Q=Q)


# ------------------------------------------------------------ AR(1) sweep
phis = np.linspace(0.0, 0.95, 39)
hb, hm, ef, eg = [], [], [], []
for phi in phis:
    r = analyse(phi ** np.arange(K + 2, dtype=float))
    hb.append(r["h_be"]); hm.append(r["h_min"]); ef.append(r["eig_fixed"]); eg.append(r["eig_gm"])
hb, hm = np.array(hb), np.array(hm)
cross = phis[np.argmin(np.abs(hb - hm))]
log(f"AR(1) signs: break-even half-spread 1-phi^2 meets manipulation bound (1+phi)/2 at phi = {cross:.3f} (theory 0.5)")
fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
fig.subplots_adjust(wspace=0.25)
axes[0].plot(phis, hb, color=style.SERIES[0], label="competitive (break-even) half-spread = σ²_innov")
axes[0].plot(phis, hm, color=style.SERIES[1], label="half-spread needed to stop manipulation")
axes[0].fill_between(phis, hb, hm, where=hm > hb, color=style.SERIES[1], alpha=0.15, lw=0)
axes[0].axvline(0.5, color=style.NEUTRAL, ls=":", lw=1)
axes[0].annotate("manipulable\nwith symmetric\nspreads", (0.68, 0.62), fontsize=8, color=style.INK)
axes[0].set_xlabel("order-sign persistence φ (AR(1))")
axes[0].set_ylabel("half-spread, units of G(1)")
axes[0].set_title("Symmetric quotes: safe only if φ ≤ ½")
axes[0].legend(fontsize=8, loc="lower left")
axes[1].plot(phis, ef, color=style.SERIES[0], label="symmetric break-even spread")
axes[1].plot(phis, eg, color=style.SERIES[2], label="Glosten–Milgrom skewed quotes")
axes[1].axhline(0, color=style.NEUTRAL, lw=1)
axes[1].set_xlabel("order-sign persistence φ (AR(1))")
axes[1].set_ylabel("min eigenvalue of cost form")
axes[1].set_title("Negative ⇒ a profitable round trip exists")
axes[1].legend(fontsize=8)
style.save(fig, "16_ar1_threshold.png")

# ------------------------------------------------------------ realistic flows
cases = {
    "fGn signs, H = 0.60": fgn_sign_acf(0.60, K + 1),
    "fGn signs, H = 0.75": fgn_sign_acf(0.75, K + 1),
    "fGn signs, H = 0.90": fgn_sign_acf(0.90, K + 1),
    "mixture of metaorder AR(1)s": mixture_acf(K + 1),
    "AR(1), φ = 0.5": 0.5 ** np.arange(K + 2, dtype=float),
    "AR(1), φ = 0.7": 0.7 ** np.arange(K + 2, dtype=float),
}
log("case                           sigma2_innov  s*     break-even h   needed h   min eig (sym)  min eig (GM skew)")
for name, acf in cases.items():
    r = analyse(acf)
    log(f"{name:30s} {r['var_innov']:8.3f}   {r['s_star']:.3f}   {r['h_be']:8.3f}    {r['h_min']:8.3f}"
        f"   {r['eig_fixed']:+10.4f}    {r['eig_gm']:+10.4f}")
    # expected P&L of a 200-trade alternating round trip under each regime
    if name in ("mixture of metaorder AR(1)s", "AR(1), φ = 0.7"):
        v = np.array([(-1) ** t for t in range(200)], float)
        v[-1] = -v[:-1].sum()
        cost_sym = 0.5 * v @ r["M"][:200, :200] @ v
        cost_gm = 0.5 * v @ r["Q"][:200, :200] @ v
        log(f"   alternating 200-trade round trip: expected cost {cost_sym:+.2f} with symmetric "
            f"break-even spread, {cost_gm:+.2f} with skewed quotes (negative = manipulator profits)")

open(style.os.path.join(style.FIG_DIR, "16_output.txt"), "w").write("\n".join(out) + "\n")
