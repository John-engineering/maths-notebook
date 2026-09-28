"""Note 08: rough volatility from nearly critical Hawkes order flow.

Order flow is a Hawkes process with branching ratio a (close to 1) and kernel
phi(t) = a * alpha * c^alpha / (c + t)^(1 + alpha). Near criticality the
integrated intensity behaves like rough volatility with H = alpha - 1/2
(Jaisson-Rosenbaum; El Euch-Fukasawa-Rosenbaum). We simulate it and measure H
with the Gatheral-Jaisson-Rosenbaum moment-scaling method.

The power-law kernel is written as a Laplace mixture of exponentials,
(c+t)^-(1+alpha) = (1/Gamma(1+alpha)) int x^alpha e^{-x c} e^{-x t} dx,
discretised on a geometric grid, so the process is Markov in M states.
"""
import numpy as np
from numba import njit
from scipy.special import gamma as G
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


def powerlaw_mixture(a, alpha, c=1.0, xmin=1e-8, xmax=1e2, M=48):
    x = np.geomspace(xmin, xmax, M)
    dlogx = np.log(x[1] / x[0])
    # phi(t) = sum_j w_j exp(-x_j t)
    w = a * alpha * c**alpha / G(1 + alpha) * x**alpha * np.exp(-x * c) * x * dlogx
    w *= a / np.sum(w / x)                  # exact branching ratio a
    return x, w


def exp_kernel(a, beta=0.05):
    return np.array([beta]), np.array([a * beta])


@njit(cache=True)
def simulate(mu, x, w, T, win, seed):
    np.random.seed(seed)
    M = len(x)
    Y = np.zeros(M)
    n_win = int(T / win)
    integ = np.zeros(n_win)
    counts = np.zeros(n_win, np.int64)
    t = 0.0
    k = 0
    B = win
    n_events = 0
    while k < n_win:
        lam_bar = mu
        for j in range(M):
            lam_bar += Y[j]
        dt = np.random.exponential(1.0) / lam_bar
        if t + dt >= B:
            h = B - t
            acc = mu * h
            for j in range(M):
                e = np.exp(-x[j] * h)
                acc += Y[j] * (1.0 - e) / x[j]
                Y[j] *= e
            integ[k] += acc
            t = B
            k += 1
            B += win
            continue
        acc = mu * dt
        lam_new = mu
        for j in range(M):
            e = np.exp(-x[j] * dt)
            acc += Y[j] * (1.0 - e) / x[j]
            Y[j] *= e
            lam_new += Y[j]
        integ[k] += acc
        t += dt
        if np.random.random() * lam_bar < lam_new:
            for j in range(M):
                Y[j] += w[j]
            counts[k] += 1
            n_events += 1
    return integ, counts, n_events


def roughness(v, lags, qs=(0.5, 1.0, 1.5, 2.0, 3.0)):
    lv = np.log(v)
    zetas, ms = [], {}
    for q in qs:
        m = np.array([np.mean(np.abs(lv[k:] - lv[:-k]) ** q) for k in lags])
        ms[q] = m
        zetas.append(np.polyfit(np.log(lags), np.log(m), 1)[0])
    zetas = np.array(zetas)
    H = np.polyfit(np.array(qs), zetas, 1)[0]       # zeta_q ~ q H
    return H, zetas, ms


def local_slopes(v, win, n_lags=60, max_lag=10000, half_width=0.25):
    """Local log-log slope of the q=1 structure function of log volatility."""
    lv = np.log(v)
    lags = np.unique(np.geomspace(1, max_lag, n_lags).astype(int))
    m1 = np.array([np.mean(np.abs(lv[k:] - lv[:-k])) for k in lags])
    x = np.log10(lags * win)
    y = np.log10(m1)
    centers, slopes = [], []
    for c in x:
        sel = np.abs(x - c) <= half_width
        if sel.sum() >= 4:
            centers.append(c)
            slopes.append(np.polyfit(x[sel], y[sel], 1)[0])
    return 10 ** np.array(centers), np.array(slopes), lags * win, m1


def main():
    from scipy.stats import skew
    T = 1_000_000.0
    WIN = 10.0
    A = 0.9995
    cases = [
        ("power-law kernel, α = 0.6", powerlaw_mixture(A, 0.6, xmin=1e-9), 0.05, 0.1, 0.6),
        ("power-law kernel, α = 0.8", powerlaw_mixture(A, 0.8, xmin=1e-9), 0.05, 0.3, 0.8),
        ("exponential kernel", exp_kernel(A, beta=0.05), 0.05, 0.5, None),
        ("power law α = 0.6, subcritical a = 0.7", powerlaw_mixture(0.7, 0.6), 30.0, None, 0.6),
    ]
    series = {}
    for name, (x, w), mu, H_th, alpha in cases:
        integ, counts, n_ev = simulate(mu, x, w, T, WIN, 5)
        # the process starts empty and takes several 10^4 time units to reach
        # its stationary level; discard the first 20% as burn-in
        v = np.maximum(integ[len(integ) // 5:], 1e-9)
        a_ = np.sum(w / x)
        tau_c = (1 - a_) ** (-1 / alpha) if alpha else 1 / (x[0] * (1 - a_))
        centers, slopes, lags_t, m1 = local_slopes(v, WIN)
        series[name] = dict(v=v, centers=centers, slopes=slopes, tau_c=tau_c, H_th=H_th)
        lv = np.log(v)
        sk = [skew(lv[k:] - lv[:-k]) for k in (1, 100)]
        # plateau: median local slope over the decade below the cutoff (or 10..300)
        lo, hi = (max(30.0, tau_c / 30), tau_c / 3) if H_th is not None else (30.0, 3000.0)
        sel = (centers >= lo) & (centers <= hi)
        plateau = np.median(slopes[sel]) if sel.any() else np.nan
        series[name]["plateau"] = (lo, hi, plateau)
        th = f"{H_th:.2f}" if H_th is not None else "n/a"
        log(f"{name:40s} {n_ev/1e6:5.1f}M events, cutoff tau_c ~ {tau_c:9.0f}; "
            f"local slope over [{lo:.0f}, {hi:.0f}] = {plateau:.3f} (theory H = {th}); "
            f"skew of dlogvol at lag 10 / 1000: {sk[0]:+.2f} / {sk[1]:+.2f}")

    # ------------------------------------------------------------ figures
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for i, (name, _, _, H_th, _) in enumerate(cases):
        d = series[name]
        ax.semilogx(d["centers"], d["slopes"], color=style.SERIES[i], label=name)
        if H_th is not None:
            ax.axhline(H_th, color=style.SERIES[i], ls=":", lw=1)
            ax.axvline(d["tau_c"], color=style.SERIES[i], ls="--", lw=0.8)
    ax.set_ylim(-0.05, 0.75)
    ax.set_xlabel("lag Δ (log)   ·   dashed: criticality cutoff τ_c   ·   dotted: α − ½")
    ax.set_ylabel("local roughness exponent")
    ax.set_title("Roughness appears between the micro scale and the cutoff")
    ax.legend(fontsize=8, loc="upper right")
    style.plain_log(ax, "x")
    style.save(fig, "08_local_H.png")

    fig, axes = plt.subplots(3, 1, figsize=(11, 7.2), sharex=True)
    fig.subplots_adjust(hspace=0.35)
    n_show = 20000
    for ax, (name, *_), col in zip(axes, cases[:3], style.SERIES[:3]):
        v = series[name]["v"][:n_show]
        ax.plot(np.arange(n_show) * WIN, np.log(v / v.mean()), color=col, lw=0.6)
        ax.set_title(name + "  (after burn-in)", fontsize=10)
        ax.set_ylabel("log vol")
    axes[-1].set_xlabel("time")
    style.save(fig, "08_paths.png")


if __name__ == "__main__":
    main()
    open(style.os.path.join(style.FIG_DIR, "08_output.txt"), "w").write("\n".join(out) + "\n")
