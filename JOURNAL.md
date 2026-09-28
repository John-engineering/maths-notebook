# Journal

A running log: what I'm thinking about, what didn't work, and questions to come
back to. Newest entries at the bottom.

---

### Entry 1: starting point

Bessembinder's "4% of stocks" result looked like a good first test of a
suspicion: a lot of what gets reported as empirical surprises about stocks is
the arithmetic of multiplication in disguise. The volatility drag $-\sigma^2/2$
is taught as a technicality. I think it is closer to the central fact of
investing.

What came out of note 01:
- Median, mean and size-biased median grow at $m-\sigma^2/2$, $m$ and $m+\sigma^2/2$.
  I like this three-row table; I haven't seen it laid out this way.
- The share of stocks creating all the wealth is $\approx\bar\Phi(\sigma\sqrt T + c)$.
- A zero-alpha simulation comes out slightly *more* concentrated than the real
  data. Worth coming back to.
- Surprise: a buy-and-hold $k$-stock portfolio needs $k \gg e^{\sigma^2 T}$ to
  track the index, which is exponential in the horizon. Buy-and-hold
  de-diversifies.

Open questions from this:
1. The cap-weighted market is buy-and-hold, so why doesn't it concentrate without
   limit? What stabilises the capital distribution curve? (Candidates: new
   listings, deaths, mean reversion by rank as in Fernholz's Atlas model.)
2. Where exactly does the "rebalancing premium" come from, and is it a free lunch?
3. Is the real data less concentrated than the null because of long-horizon
   mean reversion in idiosyncratic returns?

---

### Entry 2: rebalancing and a theorem hiding in the yield curve

**Rebalancing (note 02).** The model-free identity
EW − market = Σ log(AM/GM) + Δ log diversity is simple and exact. I checked it
on Student-t paths to 1e-15. The picture I found most useful: rebalancing is
short a straddle on the relative price and collects theta. It doesn't raise
expected wealth at all, only the median. That connects directly to note 01:
rebalancing moves you from the median-stock economy toward the mean economy.

Dead end: my first transaction-cost simulation was pure noise. The signal is a
few bp/yr, while the martingale noise over 4,000 years was about 5 bp/yr. Fix:
a control variate $\sum (w_t - \tfrac12)(dX_1 - dX_2)$, which has mean zero
exactly. After that, theory and simulation agree to within a few percent.
Lesson: when measuring a small drift, subtract the known martingale part rather
than simulating longer.

Nice surprise: the optimal no-trade band doesn't depend on vol. Both terms of
the loss scale with $\sigma^2$. Calendar intervals do depend on vol, so
volatility clustering penalises calendar rules. Confirmed: band loss unchanged
under regime-switching vol, calendar loss +25%.

Process lesson: `pkill -f pattern` in a shell whose own command line contains
the pattern kills the shell itself. That cost me two runs.

**Oscillation theorem (note 03).** Level/slope/curvature is guaranteed by
Gantmacher–Krein for any oscillatory correlation matrix; the market mode is
Perron–Frobenius. Shuffling maturities destroys the shapes but not the
eigenvalues, which is a clean demonstration that the shapes live in the ordering.
The percentages ("3 factors = 99%") reflect kernel smoothness, which curve
fitting (Nelson–Siegel is literally level/slope/curvature) and yield averaging
both inflate.

Idea to test later: on raw bond quotes, count sign changes of PC2/PC3. A
violation of the $k-1$ pattern would be genuine evidence of segmentation.

Queue of ideas, roughly by how much I want to do them:
- Trend following = the discrete Itô formula; the P&L identity makes it a
  variance-ratio bet, the mirror image of rebalancing.
- Kelly under estimation error: optimal fraction t²/(1+t²), i.e. James–Stein.
- Queue imbalance and conformal maps: the probability the next price move is up
  is a harmonic measure of a wedge.
- Martingale optimal transport: model-free bounds on exotics as an LP.
- Arcsine laws and track records.
- GARCH as a Kesten process: fitted parameters imply their own tail exponent.
- Square-root impact from a latent order book (reaction–diffusion).
- Minority game phase transition: efficiency versus volatility.
- Rough volatility from nearly-critical Hawkes processes.
- Arbitrage as cohomology: Hodge decomposition of FX log-rates.

---

### Entry 3: trend following and Itô

Summation by parts gives the trend-follower identity, and the EMA version has
the same structure. What I didn't expect: the spectral crossover period is
$\pi\sqrt{2N}$, not $N$. A 250-day trend follower is short every cycle faster
than about 70 days. Also neat: the optimal EMA decay equals the drift's AR
coefficient; the $\lambda^*=\rho$ derivation is a few lines.

The biggest payoff is the mirror with note 02. Rebalancing and trend following
are $\pm$ the same quadratic form in the path. That resolves the apparent
paradox of rebalancing "winning" in a random walk: in arithmetic P&L it
doesn't; only the log/geometric view makes it look like a free lunch.

To test later: does vol-targeting help trend following because of the
$-\sum r^2$ term, or only through Sharpe (heteroskedasticity)? In expectation,
scaling a predictable position can't create P&L in a martingale world, so any
gain has to be in the second moment. Needs a proper look, possibly a note on
when volatility targeting raises Sharpe at all (Moreira–Muir style).

---

### Entry 4: Kelly and shrinkage

The single-asset algebra is short and gives $k^* = t^2/(1+t^2)$. What I like
about it: half Kelly gets a precise meaning ($t=1$), and "parameter uncertainty"
turns out to act only through the posterior mean, because growth is linear in
$\mu$.

Surprise: plug-in shrinkage for a single strategy ($1-1/\hat t^2$) is *worse*
than fixed half Kelly. I expected it to help. Of course: Stein's phenomenon
needs three or more dimensions. That pushed me to the cross-sectional
(empirical Bayes) version, which is the most practically useful thing so far:
bet size depends on how many things you tried.

Methods note: Tweedie's formula is elegant but fragile at the extreme order
statistic, which is precisely where selection puts you. NPMLE (g-modelling)
was far more robust. The first NPMLE implementation was too slow (EM on 2000×281
per iteration × 900 rounds); binning the data fixed it.

The one-line DeMiguel bound $T > N/(\mathrm{SR}^2_{\tan} - \mathrm{SR}^2_{1/N})$
reproduces their "~3000 months for 25 assets". Satisfying.

---

### Entry 5: queues, wedges and waiting times

I expected note 06 to be a small conformal-map exercise. The waiting-time part
turned out to be the interesting bit. The wedge exit-time exponent
$\pi/(2\alpha)$ with $\alpha=\arccos(-\rho)$ means independent queues give an
infinite mean waiting time, and a renewal argument turns that into subdiffusive
prices. So the observed diffusivity of prices constrains queue correlations or
drift. Real macro/micro consistency, from one angle.

Bug worth remembering: tie-breaking when both queues empty at once gave a
spurious 60/40 in a symmetric model. A good reminder that simultaneous events
in discrete simulations need explicit handling.

Checked that my wedge formula matches the Avellaneda–Reed–Stoikov form to 1e-14.

---

### Entry 6: monitoring

Note 07 is mostly classical (arcsine law, LIL, Robbins' mixture martingale).
What's new to me is putting them together for performance evaluation. The OU
time change makes the false-alarm rate depend on $\log(T_{\max}/T_{\min})$,
which I find a clean way to say it. The always-valid threshold landing at
t ≈ 3, the same number as Harvey–Liu–Zhu's multiple-testing correction, is a
coincidence, but a memorable one.

Sobering number: an IR-0.5 manager has about 50% odds of proving skill within
30 years under honest sequential testing.

Next up: rough volatility from nearly critical Hawkes processes. This is the
most computationally demanding idea on the list, so it's worth doing carefully.

---

### Entry 7: rough volatility, and a near miss

Hawkes simulation via the Laplace-mixture representation of the power-law
kernel works well: 48 exponentials, 0.03% kernel error, 10⁸ events in about a
minute with numba.

The near miss: I measured strongly positive skew in log-vol increments and
almost wrote it up as "Hawkes rough vol predicts the time-irreversibility of
real vol". It came entirely from the burn-in (the process starts empty and
climbs for ~5·10⁴ time units). Discarding 20% killed the skew and moved H from
0.13 to 0.105. Rule for myself: with long-memory processes, look at the raw
path before computing any statistic.

α = 0.8 remains unresolved (measured ~0.2 vs 0.3). I've written that up
honestly rather than tuning the lag window until it matched.

---

### Entry 8: vol targeting

The lognormal cosine identity $\mathbb E[XY]/\sqrt{\mathbb E X^2\mathbb E Y^2} = e^{-\frac12\mathrm{Var}(\log X-\log Y)}$
does all the work in note 09. The 3/2 break-even surprised me. I'd have
guessed 1.

A theme is emerging across notes 05, 07 and 09: **every noisy input to a
position-sizing rule should be shrunk, and the shrinkage factor has the form
signal²/(signal² + noise²).** Kelly under mean uncertainty ($t^2/(1+t^2)$),
James–Stein across assets, and vol-exponent shrinkage under forecast noise
($b s^2/(b^2 s^2+e^2)$) are the same regression-to-the-mean idea in different
places.
