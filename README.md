# maths-notebook

An open-ended notebook on advanced mathematics and the stock market. It has no
syllabus. Each note starts from a question, works through the maths, and then
tests the idea numerically. Wherever I could, I tried to derive or check
something rather than repeat a textbook result.

Notes 01–23 run on simulated markets built to isolate one mechanism at a time:
a simulation tells you what a theory *predicts*. From note 24 onward the
notebook turns to **real data**, loaded by [`code/data.py`](code/data.py):
Shiller's monthly S&P 500 (1871–), daily VIX (1990–), daily S&P 500 and
NASDAQ OHLCV (1999–2018) and the Fama–French factors (1926–2018). The aim is
to find which quantitative trends survive an honest null hypothesis.

## Layout

- `notes/`: the write-ups, numbered in the order they were written
- `code/`: one script per note (`python3 code/NN_*.py` regenerates that note's figures)
- `figures/`: generated charts and printed outputs
- [`JOURNAL.md`](JOURNAL.md): running log of thoughts, dead ends and open questions

## Reading paths

The notes were written in whatever order curiosity took me, but they cluster
into a few threads. Start with **19** if you read only one: it ties several
others together.

- **Compounding and what "average" means**: 01 → 02 → 04 (rebalancing and
  trend following are mirror images) → 09
- **Learning under noise: sizing, testing, detecting**: 05 → 07 → 18 → **19**
  (all governed by the information rate SR²/2)
- **Market microstructure, and how efficiency fixes exponents**: 06 → 08 →
  10 → 11 → 16 → 20 → 23 → 22 (a conjecture, its refutation, and the resolution run
  through 10 → 11 → 16)
- **When a stylised fact is really a theorem, or weakly identified**: 01, 03
  (Perron–Frobenius, Gantmacher–Krein), 14 → 15 (tail exponents)
- **Geometry, topology and transport**: 03, 06 (conformal maps), 12 (optimal
  transport), 13 (cohomology of arbitrage)
- **Agents and emergence**: 08, 17
- **Real data: which trends survive?**: 24 onward

Recurring lessons are collected in [`JOURNAL.md`](JOURNAL.md), including three
places where a simulation overruled a confident first derivation (08, 11/16, 18).

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
| 14 | [GARCH predicts its own tail exponent](notes/14-garch-predicts-its-own-tails.md) | GARCH is a Kesten process: return tail ζ solves E[(αz²+β)^{ζ/2}] = 1 (a Cramér–Lundberg equation), ζ ≈ 2 + 2(1−p)/α²; typical fits give 3–5, but ±0.0025 in persistence moves ζ by ±0.6 |
| 15 | [A misspecified GARCH invents power-law tails](notes/15-garch-invents-tails.md) | GARCH fitted to rough *lognormal* volatility (no power law at all) reports persistence 0.997 and a "cubic-law" tail exponent of 3.7; at the top 1–5% quantiles, Pareto and lognormal-mixture tails are indistinguishable |
| 16 | [Skewed quotes, not wide spreads, stop manipulation](notes/16-skewed-quotes-stop-manipulation.md) | Competitive symmetric spreads stop alternating manipulation only if σ²_innov ≥ ½(1 + a₁ − a₂ + …) (AR(1): φ ≤ ½); Glosten–Milgrom skewed quotes are always safe (Fejér-kernel positivity) |
| 17 | [The minority game: efficient markets are crowded markets](notes/17-minority-game-efficiency-volatility.md) | Price-taking agents make the market perfectly efficient below α_c ≈ 0.34 at the cost of 14× excess volatility; impact-aware agents are 4,000× calmer but leave predictability unexploited |
| 18 | [When did my edge die?](notes/18-when-did-my-edge-die.md) | CUSUM (a drawdown rule on drift-adjusted P&L) is optimal; detection time scales as 2/ΔSR² years, so Sharpe 1 → 0 takes about 3 years at one false alarm per 20 years |
| 19 | [The Kelly bettor is a statistician](notes/19-the-kelly-bettor-is-a-statistician.md) | Synthesis: Kelly wealth = likelihood ratio (Girsanov), so Kelly growth = KL rate; note 07's always-valid test is a mixture of Kelly bettors, and CUSUM is a restarted bettor. Testing by betting |
| 20 | [The square-root law from a diffusing order book](notes/20-square-root-impact-from-diffusion.md) | Latent net liquidity obeys the heat equation; a metaorder's price solves a Volterra equation. Slow execution: diffusion square root; fast: geometric √(2Q/L). The pure-diffusion propagator (β = ½) is too forgetful for efficiency |
| 22 | [Who gets the alpha?](notes/22-who-gets-the-alpha.md) | With square-root impact and competitive investors a monopolist manager's fee is ⅓ of gross alpha and impact eats ⅔; with N managers crowding one signal, impact takes N/(N+½): 95% at N = 10 |
| 23 | [Liquidity providers' horizons set how impact decays](notes/23-liquidity-horizons-set-impact-decay.md) | Latent book with a spread of renewal rates, share below ν ∝ ν^a: impact decays as t^−(1−a) (Laplace/Tauberian, checked by Talbot inversion). Efficiency requires a = (1+γ)/2, tying takers' memory to makers' horizons |
| 24 | [150 years of the S&P 500: which trends survive an honest null?](notes/24-150-years-of-sp500.md) | **Real data.** Volatility drag = σ²/2 to 0.01%; the famous lag-1 momentum is a monthly-averaging artefact; long-run mean reversion is pre-war only and not significant once the null respects volatility clustering |
| 25 | [Does valuation predict the next decade? CAPE against a fair null](notes/25-does-cape-predict-returns.md) | **Real data.** CAPE's R² = 0.29 is unremarkable under random-walk prices with CAPE rebuilt from actual earnings (p = 0.22); it forecast well out of sample only 1960–89; the CAPE level has drifted from ~15 to ~26 since 1990 |
| 26 | [Real volatility is rough, and GARCH still gets the tails wrong](notes/26-real-volatility-is-rough.md) | **Real data.** Noise-corrected roughness of S&P/NASDAQ range volatility H ≈ 0.12–0.16 (naive fit halves it); VIX isn't rough at short lags; a t-GARCH fit reaches persistence 0.9996 and implies an infinite-variance tail (2.04) against Hill ≈ 3.4, while FIGARCH (d ≈ 0.5) fits far better |
| 27 | [Factor premia, 1926–2018: when were they proven, and did they die?](notes/27-factor-premia-evidence-and-death.md) | **Real data.** Always-valid (mixture-Kelly) evidence proves the equity premium by 1955 and value by 1964, size never; value and size roughly halve after publication; CUSUM "death" alarms fired for all three after 1993 and all three kept earning (false alarms are ~65% likely) |
| 28 | [The variance risk premium, and volatility targeting, on real data](notes/28-variance-risk-premium-and-vol-targeting.md) | **Real data.** VIX² > realised variance on 84% of days (short-variance Sharpe 1.25, skew −5.3); Gaussian-Kelly sizing is ruined by Sept 2008, and exact Kelly is 0.14 vs 0.48. Market risk-return exponent p ≈ 0.84 < 3/2, so vol targeting helps, less than theory because vol forecasts are noisy |
| 29 | [Calendar effects, tested honestly](notes/29-calendar-effects-honestly.md) | **Real data.** September is the only standout month (t = −2.8) but fails a 12-month Holm correction; James–Stein removes half of seasonal variation; Sell-in-May ≈ +4%/yr before and after publication yet never significant; NASDAQ's overnight drift is almost entirely 1999–2000 |
| 30 | [Does trend following work on the US stock market?](notes/30-trend-following-the-us-market.md) | **Real data.** 12-month TSMOM long/cash: Sharpe 0.42 → 0.56, max drawdown −85% → −44%; beats random-timing and i.i.d. nulls (p ≤ 0.002) and predicts returns beyond volatility (t = 2.2); the long/short version pays off like a straddle (+49% in 1931, +33% in 2008), as note 04 predicted |
| 31 | [How volatility shocks decay](notes/31-how-volatility-shocks-decay.md) | **Real data.** 22 VIX spikes 1990–2025 relax as a power law with β ≈ 0.33 (not exponentially; wins leave-one-out 14/22), matching rough-vol's ½ − H from note 26; spikes rise in ~2 days and halve in ~8 |

## Running

```
pip install -r requirements.txt
cd code && python3 01_bessembinder.py
```
