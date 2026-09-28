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
