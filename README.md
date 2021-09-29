# Code for "Robust Statistical Arbitrage Strategies"

## Eva Lutkebohmert, Julian Sester

# Abstract

We investigate statistical arbitrage strategies when there is ambiguity about the underlying time-discrete financial model. Pricing measures are assumed to be martingale measures calibrated to prices of liquidly traded options, whereas the set of admissible physical measures is not necessarily implied from market data. Our investigations rely on the mathematical characterization of statistical arbitrage, which was originally introduced by Bondarenko in 2003. In contrast to pure arbitrage strategies, statistical arbitrage strategies are not entirely risk-free, but the notion allows to identify strategies which are profitable on average, given the outcome of a specific sigma-algebra. Besides a characterization of robust statistical arbitrage, we also provide a super-/sub-replication theorem for the construction of statistical arbitrage strategies based on path-dependent options. In particular, we show that the range of statistical arbitrage-free prices is, in general, much tighter than the range of arbitrage-free prices. 

# Preprint

Can be found [here](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3436788)


# Content

This folder contains:

[functions.R](https://github.com/juliansester/statistical-arbitrage/blob/master/functions.R):, 
All Functions for Computation of No-Statistical Arbitrage Bounds and the corresponding replication strategies are defined here.


Further we provide the following jupyter notebooks:

1. [Example 4.7.ipynb](https://github.com/juliansester/statistical-arbitrage/blob/master/Example%204.7.ipynb), Contains Example 4.7.
2. [Example 5.1.ipynb](https://github.com/juliansester/statistical-arbitrage/blob/master/Example%205.1.ipynb), Contains Example 5.1. 
3. [Example 5.2.ipynb](https://github.com/juliansester/statistical-arbitrage/blob/master/Example%205.2.ipynb), Contains Example 5.2. 
4. [Section 5.2. Eurostoxx.ipynb](https://github.com/juliansester/statistical-arbitrage/blob/master/Section%205.2.%20Eurostoxx.ipynb), Contains the Example from Section 5.2.


## Data
Note that the data for the S&P 500 examples cannot be provided for legal reasons.

# License

MIT License

Copyright (c) 2021 Julian Sester

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




