"""Note 23: where the impact-decay exponent comes from.

Latent book = populations with renewal (cancel + re-deposit around the current
price) rates nu, liquidity shares rho(nu). Linear response to a small trade dq:
    L dp(t) = dq sum_nu rho_nu e^{-nu t} / sqrt(4 pi D t)
              + L sum_nu nu rho_nu int_0^t dp(s) e^{-nu (t-s)} ds.
Laplace: dp_hat(z) = (dq / (2 L sqrt D)) * sum rho_nu (z+nu)^-1/2 / sum rho_nu z/(z+nu).
If rho(nu) ~ nu^(a-1) at small nu with 1/2 < a < 1, then G(t) ~ t^-(1-a):
beta = 1 - a. Efficiency (note 10) needs beta = (1-gamma)/2, i.e. a = (1+gamma)/2.
Single population: G(t) = c [e^{-nu t}/sqrt(pi t) + sqrt(nu) erf(sqrt(nu t))].
"""
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


D = 1.0
L = 1.0
dq = 1.0


def u_hat(z, nus, rho):
    """Laplace transform of the linear response (dq = L = D = 1 units)."""
    z = np.asarray(z, dtype=complex)[..., None]
    num = np.sum(rho / np.sqrt(z + nus), axis=-1)
    den = np.sum(rho * z / (z + nus), axis=-1)
    return dq / (2 * L * np.sqrt(D)) * num / den


def talbot(Fhat, t, M=64):
    """Fixed-Talbot numerical inverse Laplace transform (Abate & Valko 2004).
    Valid when all singularities of Fhat lie on the non-positive real axis."""
    t = np.atleast_1d(t)
    out = np.empty(len(t))
    k = np.arange(1, M)
    theta = k * np.pi / M
    cot = 1 / np.tan(theta)
    for i, ti in enumerate(t):
        r = 2 * M / (5 * ti)
        S = r * theta * (cot + 1j)
        sig = theta + (theta * cot - 1) * cot
        val = 0.5 * np.real(Fhat(np.array([r + 0j]))[0] * np.exp(r * ti))
        val += np.sum(np.real(np.exp(ti * S) * Fhat(S) * (1 + 1j * sig)))
        out[i] = r / M * val
    return out


times = np.geomspace(1e-3, 1e7, 300)

# 1. single population: check against the closed form
nu = 1e-2
dp1 = talbot(lambda z: u_hat(z, np.array([nu]), np.array([1.0])), times)
exact = dq / (2 * L * np.sqrt(D)) * (np.exp(-nu * times) / np.sqrt(np.pi * times) + np.sqrt(nu) * erf(np.sqrt(nu * times)))
sel = times > 1e-3
log(f"single population nu={nu}: max relative error vs closed form (t > 1e-3): "
    f"{np.max(np.abs(dp1[sel] / exact[sel] - 1)):.2e}; permanent impact sqrt(nu)/2 = {np.sqrt(nu)/2:.4f}, "
    f"numerical at t=1e7: {dp1[-1]:.4f}")

# 2. power-law spectrum of renewal rates
nus = np.geomspace(1e-9, 1e2, 400)
dlog = np.log(nus[1] / nus[0])
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
fig.subplots_adjust(wspace=0.25)
axes[0].loglog(times, dp1, color=style.NEUTRAL, lw=1.5, label=f"single rate ν = {nu:g}")
res = {}
for i, a in enumerate([0.55, 0.65, 0.75, 0.9]):
    rho = nus**a * dlog                      # density nu^(a-1) d(nu) = nu^a d(log nu)
    rho /= rho.sum()
    dp = talbot(lambda z, r_=rho: u_hat(z, nus, r_), times)
    win = (times > 1e2) & (times < 1e6)       # well inside [1/nu_max, 1/nu_min]
    beta_fit = -np.polyfit(np.log(times[win]), np.log(dp[win]), 1)[0]
    res[a] = beta_fit
    axes[0].loglog(times, dp, color=style.SERIES[i], label=f"a = {a}: fitted β = {beta_fit:.3f} (theory {1-a:.2f})")
    log(f"rho(nu) ~ nu^(a-1), a = {a}: fitted decay exponent {beta_fit:.3f}, theory 1 - a = {1-a:.3f}")
rho = nus**0.3 * dlog
rho /= rho.sum()
dp_low = talbot(lambda z: u_hat(z, nus, rho), times)
win = (times > 1e2) & (times < 1e6)
b_low = -np.polyfit(np.log(times[win]), np.log(dp_low[win]), 1)[0]
log(f"a = 0.3 (< 1/2): fitted decay exponent {b_low:.3f}, theory 1/2 (the heat kernel wins)")
axes[0].set_xlim(1e-3, 1e7)
axes[0].set_xlabel("time since trade (log)")
axes[0].set_ylabel("price impact of a small trade (log)")
axes[0].set_title("Heterogeneous liquidity horizons set the decay")
axes[0].legend(fontsize=7.5, loc="lower left")
style.plain_log(axes[0], "x")

# 3. the efficiency map: beta = 1 - a  and  beta = (1 - gamma)/2
gam = np.linspace(0.0, 1.0, 50)
axes[1].plot(gam, (1 + gam) / 2, color=style.SERIES[0], label="required liquidity-horizon exponent a = (1+γ)/2")
axes[1].axhspan(0.5, 1.0, color=style.SERIES[2], alpha=0.08, lw=0)
axes[1].annotate("admissible range ½ < a < 1", (0.05, 0.93), fontsize=8, color=style.INK_2)
axes[1].plot(0.5, 0.75, marker="o", ms=7, color=style.SERIES[1])
axes[1].annotate("γ ≈ 0.5 (typical stocks) → a ≈ 0.75, β ≈ 0.25", (0.52, 0.72), fontsize=8, color=style.INK)
axes[1].set_xlabel("order-sign memory exponent γ  (C(ℓ) ~ ℓ^−γ)")
axes[1].set_ylabel("a  (share of liquidity with rate < ν ~ ν^a)")
axes[1].set_title("Efficiency ties takers' memory to makers' horizons")
axes[1].set_ylim(0.4, 1.02)
axes[1].legend(fontsize=8, loc="lower right")
style.save(fig, "23_multiscale_book.png")

open(style.os.path.join(style.FIG_DIR, "23_output.txt"), "w").write("\n".join(out) + "\n")
