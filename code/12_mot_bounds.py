"""Note 12: model-free bounds on forward-start options (martingale optimal transport).

Given all vanilla prices at T1 and T2 (hence the risk-neutral marginals of S1
and S2), the price of a payoff c(S1, S2) is only pinned down to an interval:
  [min, max] of E_pi[c] over couplings pi with those marginals and
  E_pi[S2 | S1] = S1 (martingale).
That is a linear program; its dual is the cheapest semi-static super-hedge
u1(S1) + u2(S2) + Delta(S1) (S2 - S1) >= c(S1, S2).
"""
import numpy as np
from scipy.optimize import linprog
from scipy.stats import norm
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg)
    out.append(msg)


S0 = 100.0


def discretise_lognormal(grid, sig, T, fwd=S0):
    """Probabilities on a grid with exact total mass 1 and mean fwd."""
    edges = np.concatenate([[0], 0.5 * (grid[1:] + grid[:-1]), [np.inf]])
    m = np.log(fwd) - 0.5 * sig**2 * T
    cdf = norm.cdf((np.log(np.maximum(edges, 1e-300)) - m) / (sig * np.sqrt(T)))
    p = np.diff(cdf)
    p /= p.sum()
    return tilt_to_mean(grid, p, fwd)


def tilt_to_mean(grid, p, target):
    """Exponentially tilt p so that its mean is exactly target."""
    lo, hi = -1.0, 1.0
    for _ in range(200):
        th = 0.5 * (lo + hi)
        w = p * np.exp(th * (grid - target) / target)
        w /= w.sum()
        if w @ grid < target:
            lo = th
        else:
            hi = th
    return w


def model_coupling(x, px, y, kernel_sig, dt, mix=None):
    """A martingale transition x -> y: lognormal increments, tilted so that
    E[S2 | S1 = x] = x exactly. mix = (weight, sigma_high) gives a two-state
    stochastic-volatility-like kernel."""
    P = np.zeros((len(x), len(y)))
    for i, xi in enumerate(x):
        if mix is None:
            row = discretise_lognormal(y, kernel_sig, dt, fwd=xi)
        else:
            wgt, s_hi = mix
            row = (1 - wgt) * discretise_lognormal(y, kernel_sig, dt, fwd=xi) \
                + wgt * discretise_lognormal(y, s_hi, dt, fwd=xi)
            row = tilt_to_mean(y, row, xi)
        P[i] = px[i] * row
    return P


def mot_bounds(x, px, y, py, C):
    n1, n2 = len(x), len(y)
    A_eq, b_eq = [], []
    for i in range(n1):                        # row sums = mu1
        r = np.zeros((n1, n2)); r[i, :] = 1
        A_eq.append(r.ravel()); b_eq.append(px[i])
    for j in range(n2 - 1):                    # column sums = mu2 (one redundant)
        r = np.zeros((n1, n2)); r[:, j] = 1
        A_eq.append(r.ravel()); b_eq.append(py[j])
    for i in range(n1):                        # martingale: sum_j pi_ij (y_j - x_i) = 0
        r = np.zeros((n1, n2)); r[i, :] = y - x[i]
        A_eq.append(r.ravel()); b_eq.append(0.0)
    A_eq = np.array(A_eq)
    b_eq = np.array(b_eq)
    res = {}
    for sense, sgn in [("min", 1), ("max", -1)]:
        sol = linprog(sgn * C.ravel(), A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
        assert sol.status == 0, sol.message
        res[sense] = (sgn * sol.fun, sol.x.reshape(n1, n2))
    return res


T1, T2 = 0.5, 1.0
x = np.linspace(55, 165, 45)
y = np.linspace(40, 210, 70)

worlds = {
    "flat 20% vol (Black–Scholes)": dict(sig1=0.20, ksig=0.20, mix=None),
    "regime-switching forward vol": dict(sig1=0.20, ksig=0.12, mix=(0.25, 0.38)),
}
payoffs = {
    "forward-start straddle |S₂ − S₁|": lambda X, Y: np.abs(Y - X),
    "cliquet leg (S₂/S₁ − 1)⁺": lambda X, Y: np.maximum(Y / X - 1, 0) * S0,
    "forward-start digital 1{S₂ > S₁}": lambda X, Y: (Y > X).astype(float) * S0,
    "two-date variance (log S₂/S₁)²": lambda X, Y: np.log(Y / X) ** 2 * S0,
    "separable check: S₂² − S₁²": lambda X, Y: (Y**2 - X**2) / S0,
}
X, Y = np.meshgrid(x, y, indexing="ij")
results = {}
for wname, wp in worlds.items():
    px = discretise_lognormal(x, wp["sig1"], T1)
    P = model_coupling(x, px, y, wp["ksig"], T2 - T1, wp["mix"])
    py = P.sum(0)
    # implied T2 vol of the resulting marginal (ATM, from the model coupling)
    call_atm = np.sum(py * np.maximum(y - S0, 0))
    from scipy.optimize import brentq
    bs = lambda s: S0 * norm.cdf(0.5 * s * np.sqrt(T2)) - S0 * norm.cdf(-0.5 * s * np.sqrt(T2))
    iv2 = brentq(lambda s: bs(s) - call_atm, 1e-4, 2)
    log(f"[{wname}] T2 ATM implied vol of the discretised marginal: {iv2:.4f}")
    for pname, f in payoffs.items():
        C = f(X, Y)
        r = mot_bounds(x, px, y, py, C)
        model = np.sum(P * C)
        results[(wname, pname)] = (r["min"][0], model, r["max"][0], r)
        log(f"   {pname:34s} lower {r['min'][0]:8.3f}   model {model:8.3f}   upper {r['max'][0]:8.3f}"
            f"   width/model {(r['max'][0]-r['min'][0])/max(abs(model),1e-9):.2f}")

# ------------------------------------------------------------- figure: couplings
wname = "flat 20% vol (Black–Scholes)"
pname = "forward-start straddle |S₂ − S₁|"
lo, mdl, hi, r = results[(wname, pname)]
px = discretise_lognormal(x, 0.20, T1)
P_model = model_coupling(x, px, y, 0.20, T2 - T1)
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
fig.subplots_adjust(wspace=0.08)
for ax, (title, Pm) in zip(axes, [(f"cheapest coupling: {lo:.2f}", r["min"][1]),
                                   (f"Black–Scholes coupling: {mdl:.2f}", P_model),
                                   (f"most expensive coupling: {hi:.2f}", r["max"][1])]):
    Pn = Pm / Pm.max()
    ii, jj = np.nonzero(Pn > 1e-4)
    ax.scatter(x[ii], y[jj], s=60 * Pn[ii, jj] + 2, color=style.SERIES[0], alpha=0.7,
               edgecolors="none")
    ax.plot([40, 210], [40, 210], color=style.NEUTRAL, lw=1, ls="--")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("S₁ (at 6 months)")
    ax.set_xlim(50, 170)
axes[0].set_ylabel("S₂ (at 1 year)")
fig.suptitle("Forward-start straddle: same vanilla prices, three joint laws", x=0.07, ha="left",
             fontsize=12, fontweight="bold", y=1.02)
style.save(fig, "12_couplings.png")

# the Cauchy-Schwarz ceiling: E|S2-S1| <= sqrt(E[(S2-S1)^2]) = sqrt(E S2^2 - E S1^2)
for wname in worlds:
    fv = results[(wname, "separable check: S₂² − S₁²")][1] * S0
    lo, mdl, hi, _ = results[(wname, "forward-start straddle |S₂ − S₁|")]
    log(f"[{wname}] forward variance E(S2-S1)^2 = {fv:.1f}; sqrt = {np.sqrt(fv):.2f}; straddle "
        f"lower/model/upper as fraction of sqrt: {lo/np.sqrt(fv):.3f} / {mdl/np.sqrt(fv):.3f} / "
        f"{hi/np.sqrt(fv):.3f}   (normal increments: sqrt(2/pi) = {np.sqrt(2/np.pi):.3f})")

# support size per row
for sense in ["min", "max"]:
    Pm = r[sense][1]
    k = [(Pm[i] > 1e-9 * Pm[i].sum()).sum() for i in range(len(x)) if Pm[i].sum() > 1e-9]
    log(f"straddle {sense}: support points per S1 value: median {np.median(k):.0f}, max {max(k)}")

# ------------------------------------------------------------- figure: bounds
fig, ax = plt.subplots(figsize=(9, 4.4))
names = list(payoffs.keys())[:4]
for k_, wname in enumerate(worlds):
    for i, pname in enumerate(names):
        lo, mdl, hi, _ = results[(wname, pname)]
        xpos = i + (k_ - 0.5) * 0.28
        ax.plot([xpos, xpos], [lo / mdl, hi / mdl], color=style.SERIES[k_], lw=6,
                solid_capstyle="butt", label=wname if i == 0 else None)
        ax.plot(xpos, 1.0, marker="o", color=style.INK, ms=4)
ax.set_xticks(range(len(names)))
ax.set_xticklabels([n.split(" ", 1)[0] + "\n" + n.split(" ", 1)[1] for n in names], fontsize=8)
ax.set_ylabel("price / model price")
ax.set_title("Model-free price ranges given all vanilla prices (dot = model)")
ax.legend(fontsize=8, loc="upper left")
style.save(fig, "12_bounds.png")

open(style.os.path.join(style.FIG_DIR, "12_output.txt"), "w").write("\n".join(out) + "\n")
