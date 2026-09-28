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

## Running

```
pip install -r requirements.txt
cd code && python3 01_bessembinder.py
```
