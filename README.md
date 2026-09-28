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

## Running

```
pip install -r requirements.txt
cd code && python3 01_bessembinder.py
```
