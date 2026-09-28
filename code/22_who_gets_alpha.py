"""Note 22: who gets the alpha? Fees, impact and crowding under the square-root law.

Gross alpha a per dollar; impact cost per dollar of AUM k S^delta (delta = 1/2
from the square-root law); investors are competitive (Berk-Green: net alpha 0).
  monopolist:        fee share delta/(1+delta), impact share 1/(1+delta)
  N-firm Cournot:    fee share delta/(N+delta), impact share N/(N+delta)
Checked numerically by best-response iteration.
"""
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg)
    out.append(msg)


alpha, k = 0.04, 0.01


def cournot(N, delta, iters=500):
    A = np.full(N, 0.1)
    for _ in range(iters):
        for i in range(N):
            others = A.sum() - A[i]

            def neg_rev(a):
                return -a * (alpha - k * (a + others) ** delta)

            A[i] = minimize_scalar(neg_rev, bounds=(0, 1e4), method="bounded",
                                   options=dict(xatol=1e-12)).x
    S = A.sum()
    impact = k * S**delta
    return S, impact / alpha, (alpha - impact) / alpha


log("delta   N    fee share (numeric / theory)   impact share (numeric / theory)")
fig, ax = plt.subplots(figsize=(7.4, 4.3))
Ns = np.arange(1, 21)
for i, delta in enumerate([0.5, 1.0]):
    th_fee = delta / (Ns + delta)
    ax.plot(Ns, 100 * th_fee, color=style.SERIES[i], label=f"fee share, impact ∝ AUM^{delta:g}")
    ax.plot(Ns, 100 * (1 - th_fee), color=style.SERIES[i], ls="--", lw=1.2,
            label=f"impact share, impact ∝ AUM^{delta:g}")
    for N in [1, 2, 5, 10, 20]:
        S, imp, fee = cournot(N, delta)
        ax.plot(N, 100 * fee, marker="o", ms=5, lw=0, color=style.SERIES[i])
        log(f"{delta:4.1f}  {N:3d}      {fee:.4f} / {delta/(N+delta):.4f}                {imp:.4f} / {N/(N+delta):.4f}")
ax.set_xlabel("number of managers running the same strategy")
ax.set_ylabel("% of gross alpha")
ax.set_title("Crowding hands the alpha to liquidity providers")
ax.legend(fontsize=8, loc="center right")
style.save(fig, "22_alpha_shares.png")

open(style.os.path.join(style.FIG_DIR, "22_output.txt"), "w").write("\n".join(out) + "\n")
