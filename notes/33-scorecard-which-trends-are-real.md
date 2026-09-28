# 33 · Scorecard: which quantitative trends in the stock market are real?

A summary of the real-data notes (21, 24–32). Each "trend" was tested against a
null built to match the data's awkward features (monthly averaging, volatility
clustering, persistent predictors, multiple testing, publication dates).

## Survives honest testing

| trend | key numbers | note |
|---|---|---|
| **Volatility is rough** | noise-corrected H ≈ 0.12–0.16 (S&P, NASDAQ, every subperiod) | 26 |
| **Volatility shocks decay as a power law** | 22 VIX spikes: excess ∝ (1+k/k₀)^−0.33, beats exponential (14/22 held-out); matches ½ − H | 31 |
| **Long memory in volatility** | abs(r) autocorrelation 0.24 at 20 days, power-law decay exponent 0.41; FIGARCH d ≈ 0.5 beats GARCH by 11 log-lik points | 26 |
| **Leverage asymmetry** | corr(r_t, abs(r_t+k)) ≈ −0.1 for k = 1–5; the reverse ≈ 0; VIX up in ~2 days, halves in ~8 | 26, 31 |
| **Variance risk premium** | VIX² > realised variance on 84% of days; short-variance Sharpe 1.25, skew −5.3 | 28 |
| **12-month time-series momentum** | Sharpe 0.42 → 0.56, max DD −85% → −44% (1927–2018); random-timing null p = 0.001; predicts returns beyond vol (t = 2.2) | 30 |
| **Trend following is convex** | long/short TSMOM: +49% (1931), +40% (1974), +33% (2008); quadratic smile | 30 |
| **Volatility drag = σ²/2** | 1.01% measured vs 1.00% predicted, 152 years | 24 |
| **The equity premium** | always-valid evidence reached 20× by 1955, ~4,000× by 2018 | 27 |
| **Vol targeting helps the market** | risk–return exponent p ≈ 0.84 < 3/2; Sharpe 0.43 → 0.46 (monthly), 0.18 → 0.34 (daily, variance-targeted) | 28 |
| **Credit spreads predict volatility** | BAA−AAA explains 34% of next-year log-vol variation (p < 0.001); +1pp spread ≈ +43% vol | 35 |
| **Fat tails** | daily Hill ≈ 3–3.8, monthly ≈ 3 (but see note 15: not distinguishable from a lognormal mixture) | 24, 26 |

## Real but non-stationary (regimes, not constants)

| trend | key numbers | note |
|---|---|---|
| **Stock–bond correlation** | +0.28 (1966–99), −0.32 (2000–21), +0.46 (2022–23); positive when inflation > 4% | 32 |
| **Valuation level (CAPE)** | mean ≈ 15 for 1881–1989, 26.5 since 1990 | 25 |
| **Value and size premia** | Sharpe roughly halves after publication (HML 0.62 → 0.23; SMB 0.33 → 0.14) | 27 |
| **Turn-of-the-month concentration** | 19% of days gave 68% of 1999–2018 returns, mostly 1999–2008 | 29 |
| **Overnight drift** | NASDAQ +26 bp/night vs −24 bp/day in 1999–2000; roughly balanced since 2006 | 29 |

## Does not survive

| claim | why it fails | note |
|---|---|---|
| Monthly momentum (lag-1 autocorrelation 0.26) | pure monthly-averaging artefact (Working 1960) | 24 |
| Autocorrelation at 5–21 months | disappears once the null has the real volatility history (p = 0.11) | 24 |
| Long-horizon mean reversion | pre-war only; post-1946 variance ratios go the other way; never significant | 24 |
| CAPE predicts 10-year returns | R² 0.29 is typical of random-walk prices with rebuilt CAPE (p = 0.22); out-of-sample R² negative since 1990 | 25 |
| Month-of-year seasonality | only September stands out and it fails Holm (p = 0.07); James–Stein removes half the variation | 29 |
| Sell in May | stable at about +4%/yr but never individually significant | 29 |
| Day-of-week effects | all Holm p = 1.00 | 29 |
| Buy the dip | +1.5%/yr over 3–5y after 20% drawdowns, p ≈ 0.4 (only 20 episodes) | 35 |
| Credit spreads predict returns | R² 0.005, p = 0.44 | 35 |
| Size premium | never reaches always-valid significance in 90 years | 27 |
| LPPL bubble detection | flags 24% of random walks that rose 50%; flags didn't precede drawdowns 2001–18 | 21 |
| "The premium died" alarms | CUSUM fired for all three factors after 1993; ~65% false-alarm probability; all kept earning | 27 |

## Theory predictions confirmed by the data

| theory note | prediction | data note | result |
|---|---|---|---|
| 01 | geometric = arithmetic − σ²/2 | 24 | 1.01% vs 1.00% |
| 04 | trend following is a straddle; works via low-frequency dependence | 30 | quadratic smile; profitable despite insignificant individual autocorrelations |
| 05 | Gaussian Kelly overbets skewed payoffs | 28 | mean–variance Kelly 0.48 is ruined by Sept 2008; exact Kelly 0.14 |
| 07/19 | proving an edge of Sharpe ≈ 0.4 takes decades | 27 | equity premium: 27 years to 20× evidence |
| 08 | nearly critical order flow → H ≈ 0.1 | 26, 31 | H ≈ 0.12–0.16; power-law relaxation β ≈ ½ − H |
| 09 | vol targeting helps iff p < 3/2; forecast noise shrinks the gain | 28 | p ≈ 0.84; gain about half the perfect-forecast value |
| 14/15 | long memory pushes GARCH persistence → 1 and fattens implied tails | 26 | t-GARCH persistence 0.9996 → implied tail 2.04 vs Hill 3.4 |
| 18 | detecting a dead edge takes ~2/ΔSR² years; false alarms common | 27 | 3/3 alarms, all likely false |

## The pattern

1. **The robust trends are in volatility, not in returns.** Roughness,
   power-law relaxation, long memory, leverage asymmetry and the variance risk
   premium all survive every null. Return predictability is either weak (trend
   following, t ≈ 2) or doesn't survive.
2. **Most famous return anomalies are either artefacts of measurement or
   properties of one era.** The honest shrinkage factor is about one half.
3. **The important "trends" are regime variables**: valuation levels,
   stock–bond correlation, factor premia. They are real, but a model that
   treats them as constants will be wrong at the worst moment.
4. **The theory notes did well.** Everything that mathematics predicted about
   *volatility* was confirmed, often with matching exponents across
   independent measurements.
