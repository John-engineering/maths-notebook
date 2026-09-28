"""Real market data used by the empirical notes (24 onwards).

Sources (fetched on demand, cached in data/raw/, not committed):
* Monthly S&P 500 composite 1871-: price, dividends, earnings, CPI, long rate,
  CAPE (Shiller's data as republished by github.com/datasets/s-and-p-500).
* Daily VIX 1990- (CBOE, via github.com/datasets/finance-vix).
* Daily S&P 500 and NASDAQ OHLCV 1999-2018 and monthly Fama-French factors
  1926-2018, bundled with the `arch` Python package.
"""
import os
import urllib.request
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
URLS = {
    "sp500_monthly.csv": "https://raw.githubusercontent.com/datasets/s-and-p-500/main/data/data.csv",
    "vix_daily.csv": "https://raw.githubusercontent.com/datasets/finance-vix/main/data/vix-daily.csv",
}


def _path(name):
    os.makedirs(RAW, exist_ok=True)
    p = os.path.join(RAW, name)
    if not os.path.exists(p):
        urllib.request.urlretrieve(URLS[name], p)
    return p


def sp500_monthly():
    """Monthly S&P composite. Adds nominal/real total returns where dividends
    are available (to mid-2023) and price returns throughout."""
    df = pd.read_csv(_path("sp500_monthly.csv"), parse_dates=["Date"]).set_index("Date")
    df = df.rename(columns={"SP500": "price", "Dividend": "div", "Earnings": "earn",
                            "Consumer Price Index": "cpi", "Long Interest Rate": "rate",
                            "Real Price": "real_price", "PE10": "cape"})
    for c in ["div", "earn", "cpi", "rate", "real_price", "cape"]:
        df.loc[df[c] <= 0, c] = np.nan
    df["ret_price"] = df["price"].pct_change()
    # Shiller convention: dividends are annualised, paid 1/12 per month
    df["ret_total"] = (df["price"] + df["div"] / 12) / df["price"].shift(1) - 1
    df["infl"] = df["cpi"].pct_change()
    df["ret_real"] = (1 + df["ret_total"]) / (1 + df["infl"]) - 1
    return df


def vix_daily():
    df = pd.read_csv(_path("vix_daily.csv"), parse_dates=["DATE"]).set_index("DATE")
    df.columns = [c.lower() for c in df.columns]
    return df


def sp500_daily():
    from arch.data import sp500
    return sp500.load()


def nasdaq_daily():
    from arch.data import nasdaq
    return nasdaq.load()


def french_monthly():
    """Fama-French 3 factors (percent per month), with a proper monthly index."""
    from arch.data import frenchdata
    f = frenchdata.load()
    ym = f.index.astype("int64")                  # stored as YYYYMM in nanoseconds
    idx = pd.to_datetime([f"{v // 100}-{v % 100:02d}-01" for v in ym])
    f.index = idx
    return f / 100.0
