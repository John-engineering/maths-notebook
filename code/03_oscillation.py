"""Note 03: level, slope and curvature are a theorem.

Gantmacher-Krein: if a symmetric matrix is oscillatory (totally non-negative
with some power totally positive), its k-th eigenvector has exactly k-1 sign
changes. Correlation matrices of ordered quantities whose correlation decays
with distance are usually of this kind, so the first three principal components
of a yield curve *must* look like level, slope and curvature.

Here we: (1) count sign changes over thousands of random matrices from
different families, (2) show the ordering is what matters (shuffle it and the
shapes vanish), (3) show that the *number* of factors needed is a smoothness
property of the kernel, which curve-fitting inflates.
"""
import numpy as np
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(3)
out = []


def log(msg):
    print(msg)
    out.append(msg)


mats = np.array([0.25, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30])   # years
x = np.log(mats)
m = len(mats)


def sign_changes(v, tol=1e-9):
    s = np.sign(v[np.abs(v) > tol * np.abs(v).max()])
    return int(np.sum(s[1:] != s[:-1]))


def top_eig(C, k=3):
    w, V = np.linalg.eigh(C)
    idx = np.argsort(w)[::-1]
    V = V[:, idx[:k]]
    V = V * np.sign(V.sum(0) + 1e-15)  # orient
    return w[idx], V


def to_corr(S):
    d = np.sqrt(np.diag(S))
    return S / d[:, None] / d[None, :]


# ------------------------------------------------------------- families
def fam_ou(r):
    ell = r.uniform(0.3, 6.0)
    return np.exp(-np.abs(x[:, None] - x[None, :]) / ell)


def fam_gauss(r):
    ell = r.uniform(0.5, 4.0)
    return np.exp(-((x[:, None] - x[None, :]) ** 2) / (2 * ell**2)) + 1e-10 * np.eye(m)


def fam_markov(r):
    # Gauss-Markov chain along maturity with random positive coefficients
    a = r.uniform(0.05, 1.3, m - 1)
    b = r.uniform(0.05, 1.0, m - 1)
    L = np.zeros((m, m))
    L[0, 0] = r.uniform(0.3, 2.0)
    for i in range(1, m):
        L[i] = a[i - 1] * L[i - 1]
        L[i, i] = b[i - 1]
    return to_corr(L @ L.T)


def fam_positive_factor(r):
    # all correlations positive, but no ordering structure
    k = 3
    B = np.abs(r.normal(0, 1, (m, k))) * r.uniform(0.2, 1.0, k)
    S = B @ B.T + np.diag(r.uniform(0.05, 0.5, m))
    return to_corr(S)


def fam_markov_shuffled(r):
    C = fam_markov(r)
    p = r.permutation(m)
    return C[np.ix_(p, p)]


def total_positivity_order2(C):
    """Check all 2x2 minors with i<j, k<l are >= 0 (a necessary condition)."""
    M = C[:, None, :, None] * C[None, :, None, :] - C[:, None, None, :] * C[None, :, :, None]
    i, j, k, l = np.meshgrid(*[np.arange(m)] * 4, indexing="ij")
    mask = (i < j) & (k < l)
    return bool(np.all(M[mask] >= -1e-12))


families = {
    "OU kernel in log-maturity": fam_ou,
    "Gaussian kernel in log-maturity": fam_gauss,
    "random Gauss–Markov chain": fam_markov,
    "positive factor model (no order)": fam_positive_factor,
    "Gauss–Markov, maturities shuffled": fam_markov_shuffled,
}
TRIALS = 3000
log(f"{'family':36s} P(v1 has 0 sc)  P(v2 has 1)  P(v3 has 2)  P(all three)  P(TP2)")
results = {}
for name, f in families.items():
    ok = np.zeros((TRIALS, 3), bool)
    tp2 = 0
    for t in range(TRIALS):
        C = f(rng)
        _, V = top_eig(C)
        ok[t] = [sign_changes(V[:, i]) == i for i in range(3)]
        if t < 300:
            tp2 += total_positivity_order2(C)
    results[name] = ok
    log(f"{name:36s} {ok[:,0].mean():14.3f}  {ok[:,1].mean():11.3f}  {ok[:,2].mean():11.3f}"
        f"  {ok.all(1).mean():12.3f}  {tp2/300:6.3f}")

# ------------------------------------------------------------- figure 1: shapes
fig, axes = plt.subplots(1, 4, figsize=(14, 3.6), sharey=True)
fig.subplots_adjust(wspace=0.08)
examples = [("OU kernel", fam_ou), ("random Gauss–Markov", fam_markov),
            ("positive factor model", fam_positive_factor),
            ("Gauss–Markov, shuffled", fam_markov_shuffled)]
labels = ["PC1", "PC2", "PC3"]
for ax, (title, f) in zip(axes, examples):
    C = f(np.random.default_rng(12))
    w, V = top_eig(C)
    for i in range(3):
        ax.plot(np.arange(m), V[:, i], marker="o", ms=4, color=style.SERIES[i],
                label=f"{labels[i]} ({100*w[i]/w.sum():.0f}%)")
    ax.axhline(0, color=style.NEUTRAL, lw=1)
    ax.set_xticks(np.arange(m))
    ax.set_xticklabels([f"{t:g}" for t in mats], fontsize=7)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("maturity (y)")
    ax.legend(fontsize=8, loc="lower left")
axes[0].set_ylabel("eigenvector loading")
style.save(fig, "03_shapes.png")

# ------------------------------------------------------------- figure 2: spectra
# Continuous-ish grid: 30 maturities. OU (rough) vs Gaussian (smooth) kernels
# calibrated so PC1 explains ~85%.
grid = np.log(np.geomspace(0.25, 30, 30))


def pc_share(C):
    w = np.sort(np.linalg.eigvalsh(C))[::-1]
    return w / w.sum()


def calibrate(kernel, target=0.85):
    lo, hi = 0.01, 100.0
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if pc_share(kernel(mid))[0] < target:
            lo = mid
        else:
            hi = mid
    return mid


ou = lambda ell: np.exp(-np.abs(grid[:, None] - grid[None, :]) / ell)
ga = lambda ell: np.exp(-((grid[:, None] - grid[None, :]) ** 2) / (2 * ell**2)) + 1e-12 * np.eye(30)
ell_ou, ell_ga = calibrate(ou), calibrate(ga)
sh_ou, sh_ga = pc_share(ou(ell_ou)), pc_share(ga(ell_ga))
log(f"calibrated to PC1=85%: OU ell={ell_ou:.2f}, top-3 share {sh_ou[:3].sum():.4f}; "
    f"Gaussian ell={ell_ga:.2f}, top-3 share {sh_ga[:3].sum():.4f}")

# the curve-fitting experiment: rough "true" curve changes + measurement noise,
# then smoothed by least squares onto a few Nelson-Siegel-like basis functions
tau = np.exp(grid)
T = 5000
Lc = np.linalg.cholesky(ou(ell_ou) + 1e-12 * np.eye(30))
dy_true = rng.standard_normal((T, 30)) @ Lc.T
dy_obs = dy_true + 0.15 * rng.standard_normal((T, 30))
lam = 0.6


def ns_basis(tau, k):
    f1 = np.ones_like(tau)
    f2 = (1 - np.exp(-lam * tau)) / (lam * tau)
    f3 = f2 - np.exp(-lam * tau)
    cols = [f1, f2, f3]
    if k >= 4:
        lam2 = 0.1
        f4 = (1 - np.exp(-lam2 * tau)) / (lam2 * tau) - np.exp(-lam2 * tau)
        cols.append(f4)
    return np.stack(cols[:k], 1)


def smooth(dy, k):
    B = ns_basis(tau, k)
    P = B @ np.linalg.pinv(B)
    return dy @ P.T


def emp_share(dy):
    w = np.sort(np.linalg.eigvalsh(np.cov(dy.T)))[::-1]
    return w / w.sum()


shares = {
    "true (OU kernel)": emp_share(dy_true),
    "observed (+ noise)": emp_share(dy_obs),
    "Svensson-fitted (4 basis fns)": emp_share(smooth(dy_obs, 4)),
    "Nelson–Siegel-fitted (3)": emp_share(smooth(dy_obs, 3)),
}
for k_, v in shares.items():
    log(f"  {k_:30s} PC1 {v[0]:.3f}  top-3 {v[:3].sum():.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
fig.subplots_adjust(wspace=0.25)
k = np.arange(1, 31)
axes[0].semilogy(k, sh_ou, marker="o", ms=4, color=style.SERIES[0], label="OU kernel (rough)")
axes[0].semilogy(k, sh_ga, marker="o", ms=4, color=style.SERIES[1], label="Gaussian kernel (smooth)")
axes[0].set_ylim(1e-7, 1.2)
axes[0].set_xlabel("principal component")
axes[0].set_ylabel("share of variance (log)")
axes[0].set_title("Same PC1 (85%), very different tails")
axes[0].legend()
for i, (k_, v) in enumerate(shares.items()):
    axes[1].plot(k[:8], 100 * np.cumsum(v[:8]), marker="o", ms=4, color=style.SERIES[i], label=k_)
axes[1].set_xlabel("number of components")
axes[1].set_ylabel("cumulative % of variance")
axes[1].set_title("Curve fitting manufactures 'three factors'")
axes[1].legend(loc="lower right")
style.save(fig, "03_spectra.png")

open(style.os.path.join(style.FIG_DIR, "03_output.txt"), "w").write("\n".join(out) + "\n")
