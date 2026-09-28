# maths-notebook

An open-ended notebook on advanced mathematics and the stock market. It has no
syllabus. Each note starts from a question, works through the maths, and then
tests the idea numerically. Wherever I could, I tried to derive or check
something rather than repeat a textbook result.

No live market data is available in this environment, so the experiments run on
simulated markets built to isolate one mechanism at a time. That is a
limitation, but it is also useful: a simulation tells you what a theory
*predicts*, which is what you need to know before you can tell whether real data
disagrees with it.

## Layout

- `notes/`: the write-ups, numbered in the order they were written
- `code/`: one script per note (`python3 code/NN_*.py` regenerates that note's figures)
- `figures/`: generated charts and printed outputs
- [`JOURNAL.md`](JOURNAL.md): running log of thoughts, dead ends and open questions

## Notes

| # | title | one-line idea |
|---|---|---|
| 01 | [Three growth rates, and why 4% of stocks make all the money](notes/01-three-growth-rates.md) | The median stock, the mean and the median dollar grow at μ−σ²/2, μ and μ+σ²/2; Bessembinder's skewness follows from σ√T |
| 02 | [Where the rebalancing premium comes from](notes/02-rebalancing-premium.md) | Exact identity: EW − market = Σ log(AM/GM) + Δ diversity. Rebalancing is short a straddle on relative prices; optimal no-trade band (3c/16)^⅓ is independent of vol |
| 03 | [Level, slope and curvature are a theorem](notes/03-level-slope-curvature-theorem.md) | Gantmacher–Krein forces k−1 sign changes in the k-th PC of any oscillatory correlation matrix; "three factors explain 99%" is mostly curve-fitting smoothness |
| 04 | [Trend following is the discrete Itô formula](notes/04-trend-following-is-ito.md) | P&L = ½[(long move)² − Σ(short moves)²] exactly: a variance-ratio bet and a long straddle. Spectral crossover at π√(2N); optimal EMA decay = drift persistence; rebalancing is the mirror image |
| 05 | [Kelly with an estimated edge](notes/05-kelly-with-estimated-edge.md) | Optimal fraction t²/(1+t²) (half Kelly ⇔ t = 1); plug-in mean–variance needs SR²T > N; under selection, shrink using the cross-section of all backtests (empirical Bayes) |
| 06 | [The next tick is a harmonic measure](notes/06-queues-and-conformal-maps.md) | Queue race = Brownian motion in a wedge of angle arccos(−ρ): P(up) = angle/α; waiting times have tail exponent π/(2α), so diffusive prices require anti-correlated queues |
| 07 | [Luck, track records, and the t-statistic you keep checking](notes/07-luck-and-monitoring.md) | Arcsine law: 10% of zero-skill managers show 27 unbroken years ahead; monthly monitoring turns a 5% test into 35%; Robbins' mixture martingale gives an always-valid t ≈ 3 |
| 08 | [Rough volatility from nearly critical order flow](notes/08-rough-volatility-from-critical-order-flow.md) | Hawkes order flow with branching ratio → 1 and kernel ∝ t^−(1+α) gives vol roughness H = α − ½; simulated α = 0.6 → H ≈ 0.105. Includes a burn-in artefact that nearly became a false "finding" |
| 09 | [When does volatility targeting work? The 3/2 rule](notes/09-when-vol-targeting-works.md) | Sharpe efficiency = exp(−½ Var[log risk − log Sharpe]); vol targeting beats constant exposure iff expected returns scale more weakly than σ^{3/2}; forecast noise shrinks the optimal response |
| 10 | [Impact kernels: efficiency picks the exponent, Bochner forbids manipulation](notes/10-impact-kernels-bochner.md) | No profitable round trip ⇔ Ĝ ≥ 0; convexity also rules out interim sells; long-memory order flow requires impact decay ℓ^−(1−γ)/2 for diffusive prices, and that kernel is manipulation-free too (with a caveat: see 11) |
| 11 | [An efficient price can still be manipulated, unless there is a spread](notes/11-efficiency-does-not-imply-no-manipulation.md) | Martingale kernel = cumulative order-flow surprise; its Toeplitz symbol is s + Σ a_j D_j(ω) (Dirichlet kernels), negative at ω = π unless the spread s ≥ a₁ − a₂ + a₃ − … |
| 12 | [What vanilla prices don't tell you](notes/12-forward-start-options-and-transport.md) | Martingale optimal transport as an LP: vanillas fix the forward variance exactly, so forward-start straddles are pure forward-kurtosis bets (range 0.38–0.98 of √variance, BS = √(2/π)); forward digitals are almost unconstrained |
| 13 | [Cross rates, Kirchhoff's laws and cohomology](notes/13-cross-rates-kirchhoff-cohomology.md) | Implied cross-rate variance = effective resistance; optimal route weights = current flow (EUR/JPY's own quote gets 32%); Hodge splits arbitrage into triangular (curl) and hole-bound (harmonic) parts, and USD-centred FX has no holes |

## Running

```
pip install -r requirements.txt
cd code && python3 01_bessembinder.py
```
