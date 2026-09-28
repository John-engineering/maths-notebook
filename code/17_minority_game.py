"""Note 17: the minority game: efficiency versus volatility.

N agents, each with S = 2 random lookup-table strategies over P information
states. Each round: info mu is drawn uniformly, agents play their best-scored
strategy, attendance A = sum of actions, the minority wins.
Scores: U_{i,s} += -a_{i,s}^mu (A - eta (a_i - a_{i,s}^mu)) / P
eta = 0: price-takers (ignore own impact); eta = 1: agents account for it.
Control parameter alpha = P / N. Measures: volatility sigma^2/N = <A^2>/N and
predictability H/N = (1/P) sum_mu <A | mu>^2 / N.
"""
import numpy as np
from numba import njit
import matplotlib.pyplot as plt
import style

style.apply()
out = []


def log(msg):
    print(msg, flush=True)
    out.append(msg)


@njit(cache=True)
def run_mg(N, P, eta, t_eq, t_meas, seed):
    np.random.seed(seed)
    strat = np.where(np.random.random((N, 2, P)) < 0.5, -1, 1)
    U = np.zeros((N, 2))
    # tiny random initial bias to break ties
    for i in range(N):
        U[i, 0] = 1e-6 * np.random.randn()
    sumA2 = 0.0
    condA = np.zeros(P)
    condN = np.zeros(P)
    for t in range(t_eq + t_meas):
        mu = np.random.randint(P)
        A = 0
        for i in range(N):
            s = 0 if U[i, 0] >= U[i, 1] else 1
            A += strat[i, s, mu]
        for i in range(N):
            s_used = 0 if U[i, 0] >= U[i, 1] else 1
            a_used = strat[i, s_used, mu]
            for s in range(2):
                a_s = strat[i, s, mu]
                U[i, s] -= a_s * (A - eta * (a_used - a_s)) / P
        if t >= t_eq:
            sumA2 += A * A
            condA[mu] += A
            condN[mu] += 1
    sigma2 = sumA2 / t_meas / N
    H = 0.0
    for mu in range(P):
        if condN[mu] > 0:
            H += (condA[mu] / condN[mu]) ** 2
    H = H / P / N
    return sigma2, H


def main():
    N = 301
    ms = np.arange(1, 13)
    res = {}
    for eta in (0.0, 1.0):
        rows = []
        for m in ms:
            P = 2**m
            t_eq = max(200 * P, 20000)
            t_meas = max(200 * P, 20000)
            s2, H = [], []
            for seed in range(4):
                a, b = run_mg(N, P, eta, t_eq, t_meas, 100 * m + seed)
                s2.append(a)
                H.append(b)
            rows.append((P / N, np.mean(s2), np.mean(H)))
            log(f"eta={eta:.0f}  m={m:2d}  alpha={P/N:7.4f}  sigma^2/N={np.mean(s2):8.4f}  H/N={np.mean(H):.4f}")
        res[eta] = np.array(rows)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    fig.subplots_adjust(wspace=0.25)
    for i, eta in enumerate((0.0, 1.0)):
        r = res[eta]
        lab = "price-takers (η = 0)" if eta == 0 else "impact-aware (η = 1)"
        axes[0].loglog(r[:, 0], r[:, 1], marker="o", ms=4, color=style.SERIES[i], label=lab)
        axes[1].semilogx(r[:, 0], r[:, 2], marker="o", ms=4, color=style.SERIES[i], label=lab)
    axes[0].axhline(1.0, color=style.NEUTRAL, ls="--", lw=1)
    axes[0].annotate("random agents (σ²/N = 1)", (0.02, 1.1), fontsize=8, color=style.INK_2)
    for ax in axes:
        ax.axvline(0.3374, color=style.NEUTRAL, ls=":", lw=1)
    axes[0].set_xlabel("α = P / N (information complexity per agent)")
    axes[0].set_ylabel("volatility σ² / N")
    axes[0].set_title("Volatility: crowding below α_c")
    axes[0].legend(fontsize=8)
    style.plain_log(axes[0])
    style.plain_log(axes[0], "x")
    axes[1].set_xlabel("α = P / N")
    axes[1].set_ylabel("predictability H / N")
    axes[1].set_title("Predictability: zero below α_c ≈ 0.34")
    axes[1].legend(fontsize=8)
    style.plain_log(axes[1], "x")
    style.save(fig, "17_minority_game.png")

    open(style.os.path.join(style.FIG_DIR, "17_output.txt"), "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
