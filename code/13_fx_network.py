"""Note 13: cross rates, Kirchhoff and cohomology.

Quotes w_e = phi_j - phi_i + noise_e (log rates) on the edges of a quote graph.
* Weighted least squares for phi (the consistent prices) = a Laplacian solve;
  Var(implied log cross rate i/j) = effective resistance R_eff(i, j) with
  resistances sigma_e^2. The optimal weight on each route = current flow.
* Residual = curl part (triangular arbitrage) + harmonic part (arbitrage around
  cycles no triangle detects). Under pure noise, weighted residual sum of
  squares ~ chi^2 with (E - V + 1) dof, split as rank(curl) + beta_1.
"""
import itertools
import numpy as np
import matplotlib.pyplot as plt
import style

style.apply()
rng = np.random.default_rng(13)
out = []


def log(msg):
    print(msg)
    out.append(msg)


ccy = ["USD", "EUR", "JPY", "GBP", "CHF", "AUD", "CAD", "NZD", "SEK", "NOK", "MXN", "ZAR", "TRY"]
idx = {c: i for i, c in enumerate(ccy)}
quotes = [  # (base, quote, noise sd in bp): roughly half-spread-sized quote noise
    ("EUR", "USD", 0.5), ("USD", "JPY", 0.6), ("GBP", "USD", 0.8), ("USD", "CHF", 1.0),
    ("AUD", "USD", 1.0), ("USD", "CAD", 1.0), ("NZD", "USD", 1.5), ("USD", "SEK", 3.0),
    ("USD", "NOK", 3.0), ("USD", "MXN", 5.0), ("USD", "ZAR", 8.0), ("USD", "TRY", 15.0),
    ("EUR", "JPY", 1.0), ("EUR", "GBP", 1.0), ("EUR", "CHF", 1.2), ("EUR", "SEK", 2.5),
    ("EUR", "NOK", 2.5), ("EUR", "TRY", 20.0), ("EUR", "ZAR", 12.0), ("GBP", "JPY", 2.0),
    ("AUD", "JPY", 2.0), ("AUD", "NZD", 2.5), ("CAD", "JPY", 3.0), ("NOK", "SEK", 4.0),
    ("CHF", "JPY", 2.5), ("GBP", "CHF", 3.0),
]
V, E = len(ccy), len(quotes)
B = np.zeros((E, V))
sig = np.zeros(E)
for e, (a, b, s) in enumerate(quotes):
    B[e, idx[a]] = -1.0
    B[e, idx[b]] = 1.0
    sig[e] = s * 1e-4
W = np.diag(1 / sig**2)
L = B.T @ W @ B                                # weighted Laplacian (conductance 1/sigma^2)
Lp = np.linalg.pinv(L)


def r_eff(i, j):
    d = np.zeros(V); d[i], d[j] = 1, -1
    return d @ Lp @ d


# ---------------------------------------------- 1. estimation variance = R_eff
phi_true = rng.normal(0, 1, V)
n_sims = 20000
noise = rng.standard_normal((n_sims, E)) * sig
wq = (B @ phi_true)[None, :] + noise
phi_hat = (Lp @ B.T @ W @ wq.T).T
pairs = [("EUR", "JPY"), ("GBP", "JPY"), ("TRY", "JPY"), ("ZAR", "TRY"), ("MXN", "NZD"), ("SEK", "NOK")]
log("pair       direct quote sd   network sd (sim)   sqrt(R_eff)   improvement")
direct = {(a, b): s for a, b, s in quotes}
for a, b in pairs:
    i, j = idx[a], idx[b]
    est = phi_hat[:, j] - phi_hat[:, i]
    sd_sim = np.std(est - (phi_true[j] - phi_true[i])) * 1e4
    sd_th = np.sqrt(r_eff(i, j)) * 1e4
    d = direct.get((a, b), direct.get((b, a)))
    ds = f"{d:5.1f} bp" if d else "  (none)"
    imp = f"{d/sd_th:5.2f}x" if d else "   -"
    log(f"{a}/{b}    {ds}          {sd_sim:6.2f} bp          {sd_th:6.2f} bp     {imp}")

# ---------------------------------------------- 2. current flow = route weights
def current_flow(a, b):
    d = np.zeros(V); d[idx[a]], d[idx[b]] = 1, -1
    pot = Lp @ d
    return (B @ pot) / sig**2 * -1          # current along each quoted edge (a -> b positive)


pos = {}
angles = np.linspace(0, 2 * np.pi, V - 1, endpoint=False)
pos["USD"] = np.array([0.0, 0.0])
for k, c in enumerate(ccy[1:]):
    pos[c] = np.array([np.cos(angles[k]), np.sin(angles[k])]) * (1.0 if c in ("EUR", "JPY", "GBP", "CHF", "AUD", "CAD", "NZD") else 1.6)
fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))
for ax, (a, b) in zip(axes, [("EUR", "JPY"), ("TRY", "JPY")]):
    cur = current_flow(a, b)
    total = np.abs(cur).max()
    for e, (u, v_, s) in enumerate(quotes):
        p, q = pos[u], pos[v_]
        f = abs(cur[e])
        ax.plot([p[0], q[0]], [p[1], q[1]], color=style.NEUTRAL, lw=0.8, zorder=1)
        if f > 0.02:
            ax.plot([p[0], q[0]], [p[1], q[1]], color=style.SERIES[0], lw=1 + 9 * f, alpha=0.8, zorder=2)
            mid = 0.5 * (p + q)
            ax.annotate(f"{100*f:.0f}%", mid, fontsize=7, color=style.INK, ha="center", zorder=4)
    for c, p in pos.items():
        col = style.SERIES[1] if c in (a, b) else style.SURFACE
        ax.scatter(*p, s=420, color=col, edgecolors=style.INK_2, zorder=3)
        ax.annotate(c, p, ha="center", va="center", fontsize=8, zorder=5)
    ax.set_title(f"How to price {a}/{b}: current flow = weight on each quote", fontsize=11)
    ax.set_aspect("equal")
    ax.axis("off")
style.save(fig, "13_current_flow.png")
for a, b in [("EUR", "JPY"), ("TRY", "JPY")]:
    cur = current_flow(a, b)
    d_e = [e for e, (u, v_, s) in enumerate(quotes) if {u, v_} == {a, b}]
    via = sorted([(abs(cur[e]), quotes[e][:2]) for e in range(E)], reverse=True)[:5]
    log(f"flow {a}->{b}: direct quote carries {100*abs(cur[d_e[0]]) if d_e else 0:.0f}% ; top edges "
        + ", ".join(f"{u}/{v_} {100*f:.0f}%" for f, (u, v_) in via))


# ---------------------------------------------- 3. Hodge decomposition
def triangles(edges_list, nV):
    adj = {}
    for e, (u, v_) in enumerate(edges_list):
        adj[(u, v_)] = (e, 1.0)
        adj[(v_, u)] = (e, -1.0)
    tris = []
    for i, j, k in itertools.combinations(range(nV), 3):
        if (i, j) in adj and (j, k) in adj and (k, i) in adj:
            row = np.zeros(len(edges_list))
            for (u, v_) in [(i, j), (j, k), (k, i)]:
                e, sgn = adj[(u, v_)]
                row[e] += sgn
            tris.append(row)
    return np.array(tris)


def hodge(Bm, sigv, w, C):
    """Whitened orthogonal decomposition of w into gradient + curl + harmonic."""
    Wh = 1 / sigv
    wt = w * Wh
    Bt = Bm * Wh[:, None]
    Ct = C / Wh[None, :] if len(C) else np.zeros((0, len(w)))
    Pg = Bt @ np.linalg.pinv(Bt)
    grad = Pg @ wt
    if len(Ct):
        Pc = Ct.T @ np.linalg.pinv(Ct.T)
        curl = Pc @ wt
    else:
        curl = np.zeros_like(wt)
    harm = wt - grad - curl
    return grad, curl, harm


edges_ij = [(idx[a], idx[b]) for a, b, _ in quotes]
C = triangles(edges_ij, V)
rank_c = np.linalg.matrix_rank(C)
beta1 = (E - V + 1) - rank_c
log(f"quote graph: V={V}, E={E}, independent cycles E-V+1 = {E-V+1}, triangles {len(C)} "
    f"(rank {rank_c}), harmonic dimension beta_1 = {beta1}")
ss_curl, ss_harm = [], []
for k in range(4000):
    w = B @ phi_true + rng.standard_normal(E) * sig
    g, c, h = hodge(B, sig, w, C)
    ss_curl.append(c @ c)
    ss_harm.append(h @ h)
log(f"null: E[curl SS] = {np.mean(ss_curl):.2f} (chi2 dof {rank_c}); E[harmonic SS] = {np.mean(ss_harm):.3f} (dof {beta1})")

# a venue ring (crypto-style): 6 venues quoting BTC, transfers only between
# neighbours on a ring plus one hub link. Cycles here are NOT filled by triangles.
nv = 6
ring = [(i, (i + 1) % nv) for i in range(nv)] + [(0, 3)]
Bv = np.zeros((len(ring), nv))
for e, (u, v_) in enumerate(ring):
    Bv[e, u], Bv[e, v_] = -1, 1
sv = np.full(len(ring), 2e-4)
Cv = triangles(ring, nv)
bv = (len(ring) - nv + 1) - (np.linalg.matrix_rank(Cv) if len(Cv) else 0)
log(f"venue graph: V={nv}, E={len(ring)}, triangles {len(Cv)}, harmonic dimension {bv}")
# inject a 10 bp arbitrage around the 4-cycle 0-1-2-3-0 (no triangle sees it)
phi_v = rng.normal(0, 0.01, nv)
w = Bv @ phi_v
arb = np.zeros(len(ring))
for e in [0, 1, 2]:
    arb[e] += 10e-4 / 4
arb[6] -= 10e-4 / 4               # edge (0,3) traversed backwards closes the cycle
w_arb = w + arb
g, c, h = hodge(Bv, sv, w_arb, Cv)
log(f"venue ring with 10bp 4-cycle arbitrage: curl SS {c@c:.3f}, harmonic SS {h@h:.3f} "
    f"(in units of quote variance); triangle checks find nothing because there are no triangles")
# and a triangular arbitrage in the FX graph
w_fx = B @ phi_true
tri_edges = [e for e, (a, b, s) in enumerate(quotes) if {a, b} in ({"EUR", "GBP"}, {"GBP", "JPY"}, {"EUR", "JPY"})]
bump = np.zeros(E)
bump[tri_edges[0]] += 3e-4
g, c, h = hodge(B, sig, w_fx + bump, C)
log(f"FX graph with a 3bp mispricing on one EUR/GBP/JPY leg: curl SS {c@c:.2f}, harmonic SS {h@h:.4f} (null expectation of curl SS: {rank_c})")

open(style.os.path.join(style.FIG_DIR, "13_output.txt"), "w").write("\n".join(out) + "\n")
