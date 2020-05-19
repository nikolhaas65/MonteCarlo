# FE621 Notes

## Lattice-based Methods

Tree methods equations are derived by moment matching for mean and variance. Note that real world probabilities do not play any role in valuing an option in a tree method.

### Binomial Trees

* Additive tree (<img src="https://tex.s2cms.ru/svg/log%20S_t" alt="log S_t" />)

* Trigeogis tree: <img src="https://tex.s2cms.ru/svg/%5CDelta%20x_u%20%3D%20%5CDelta%20x_d" alt="\Delta x_u = \Delta x_d" />

<img src="https://tex.s2cms.ru/svg/%0A%20%20%20%20p_u%20%5CDelta%20x%20-%20p_d%20%5CDelta%20x%20%3D%20%5Cmu%20%5CDelta%20t%5C%5C%0A%20%20%20%20p_u%20%5CDelta%20x%5E2%20%2B%20p_d%20%5CDelta%20x%20%5E2%20-%20(%5Cmu%20%5CDelta%20t)%5E2%20%3D%20%5Csigma%5E2%20%5CDelta%20t%5C%5C%0A%20%20%20%20p_u%20%2B%20p_d%20%3D%201%5C%5C%0A%20%20%20%20%5Cmu%20%3D%20r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%0A" alt="
    p_u \Delta x - p_d \Delta x = \mu \Delta t\\
    p_u \Delta x^2 + p_d \Delta x ^2 - (\mu \Delta t)^2 = \sigma^2 \Delta t\\
    p_u + p_d = 1\\
    \mu = r-\frac{\sigma^2}{2}
" />

* Jarrow Rudd tree: <img src="https://tex.s2cms.ru/svg/p_u%20%3D%20p_d" alt="p_u = p_d" />

* Multiplicative tree (<img src="https://tex.s2cms.ru/svg/S_t" alt="S_t" />)

* Cox-Ross-Rubinstein (CRR) tree:  similar to Trigeogis tree, but different in the first moment  and involves approximation. <img src="https://tex.s2cms.ru/svg/u%20%3D%20e%5E%7B%5Csigma%20%5Csqrt%7B%5CDelta%20t%7D%7D" alt="u = e^{\sigma \sqrt{\Delta t}}" />

  
<img src="https://tex.s2cms.ru/svg/%0Ap_u%20u%20%2B%20p_d%20d%20%3D%20e%5E%7Br%20%5CDelta%20t%7D%5C%5C%0Ap_u%20%5CDelta%20x%5E2%20%2B%20p_d%20%5CDelta%20x%20%5E2%20-%20(%5Cmu%20%5CDelta%20t)%5E2%20%3D%20%5Csigma%5E2%20%5CDelta%20t%5C%5C%0Au%20%3D%20%5Cfrac%7B1%7D%7Bd%7D%3De%5E%7B%5CDelta%20x%7D%5C%5C%0Ap_u%20%2B%20p_d%20%3D%201%5C%5C%0A%5Cmu%20%3D%20r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%0A" alt="
p_u u + p_d d = e^{r \Delta t}\\
p_u \Delta x^2 + p_d \Delta x ^2 - (\mu \Delta t)^2 = \sigma^2 \Delta t\\
u = \frac{1}{d}=e^{\Delta x}\\
p_u + p_d = 1\\
\mu = r-\frac{\sigma^2}{2}
" />


### Trinomial Trees

<img src="https://tex.s2cms.ru/svg/%0Ap_u%20%5CDelta%20x%20-%20p_d%20%5CDelta%20x%20%3D%20%5Cmu%20%5CDelta%20t%5C%5C%0Ap_u%20%5CDelta%20x%5E2%20%2B%20p_d%20%5CDelta%20x%20%5E2%20-%20(%5Cmu%20%5CDelta%20t)%5E2%20%3D%20%5Csigma%5E2%20%5CDelta%20t%5C%5C%0Ap_u%20%2Bp_m%2B%20p_d%20%3D%201%5C%5C%0A%5Cmu%20%3D%20r-%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%0A" alt="
p_u \Delta x - p_d \Delta x = \mu \Delta t\\
p_u \Delta x^2 + p_d \Delta x ^2 - (\mu \Delta t)^2 = \sigma^2 \Delta t\\
p_u +p_m+ p_d = 1\\
\mu = r-\frac{\sigma^2}{2}
" />

In order to make probabilities to be numbers between 0 and 1, we have a sufficient condition: 
<img src="https://tex.s2cms.ru/svg/%0A%5CDelta%20x%20%3E%20%5Csigma%20%5Csqrt%7B3%20%5CDelta%20t%7D%0A" alt="
\Delta x &gt; \sigma \sqrt{3 \Delta t}
" />
Any <img src="https://tex.s2cms.ru/svg/%5CDelta%20x" alt="\Delta x" /> with this property produces a convergent tree.

### Finite Difference Method

More general to the addictive trinomial tree method.

### Convergence Comparison

| method             | rate of convergence                                          | convergence condition                    |
| :----------------- | ------------------------------------------------------------ | ---------------------------------------- |
| binomial tree      | <img src="https://tex.s2cms.ru/svg/O((%5CDelta%20x)%5E2%20%2B%20%5CDelta%20t)" alt="O((\Delta x)^2 + \Delta t)" />                                 | NA                                       |
| trinomial tree     | <img src="https://tex.s2cms.ru/svg/O((%5CDelta%20x)%5E2%20%2B%20%5CDelta%20t)" alt="O((\Delta x)^2 + \Delta t)" />                                 | <img src="https://tex.s2cms.ru/svg/%5CDelta%20x%20%5Cgeq%20%5Csigma%20%5Csqrt%7B3%20%5CDelta%20t%7D" alt="\Delta x \geq \sigma \sqrt{3 \Delta t}" /> |
| explicit FDM       | <img src="https://tex.s2cms.ru/svg/O((%5CDelta%20x)%5E2%20%2B%20%5CDelta%20t)" alt="O((\Delta x)^2 + \Delta t)" />                                 | <img src="https://tex.s2cms.ru/svg/%5CDelta%20x%20%5Cgeq%20%5Csigma%20%5Csqrt%7B3%20%5CDelta%20t%7D" alt="\Delta x \geq \sigma \sqrt{3 \Delta t}" /> |
| implicit FDM       | <img src="https://tex.s2cms.ru/svg/O((%5CDelta%20x)%5E2%20%2B%20%5CDelta%20t)" alt="O((\Delta x)^2 + \Delta t)" />                                 | stable                                   |
| Crank-Nicolson FDM | <img src="https://tex.s2cms.ru/svg/O((%5CDelta%20x)%5E2%20%2B%20(%5Cfrac%7B%5CDelta%20t%7D%7B2%7D)%5E2)" alt="O((\Delta x)^2 + (\frac{\Delta t}{2})^2)" />                   | stable                                   |
| Monte Carlo        | <img src="https://tex.s2cms.ru/svg/O%5Cleft(max%5Cleft(%5CDelta%20t%2C%20%5Cfrac%7B%5Csigma%7D%7B%5Csqrt%7BN_x%7D%7D%5Cright)%5Cright)" alt="O\left(max\left(\Delta t, \frac{\sigma}{\sqrt{N_x}}\right)\right)" /> | stable                                   |

## Variance Reduction

|      method        |              explanation                    |
| :-----------------:| :-----------------------------------------: |
| antithetic variates| the payoff of the antithetic pair <img src="https://tex.s2cms.ru/svg/(X_1%2CX_2)" alt="(X_1,X_2)" />, <img src="https://tex.s2cms.ru/svg/f(X_1)%2C%20f(X_2)" alt="f(X_1), f(X_2)" /> is negative correlated. A sufficient condition is to make payoff function monotone |
| delta-based control variates |                                                              |

## Risk-Neutral Measure

Risk-neutral measures make it easy to express the value of a derivative in a formula. 

<img src="https://tex.s2cms.ru/svg/H_0%20%3D%20P(0%2CT)%20E_Q%5BH_T%5D" alt="H_0 = P(0,T) E_Q[H_T]" />

where the risk-neutral measure is denoted by Q. This can be re-stated in terms of the physical measure P as

<img src="https://tex.s2cms.ru/svg/H_0%20%3D%20P(0%2CT)E_P%5B%5Cfrac%7BdQ%7D%7BdP%7D%20H_T%5D" alt="H_0 = P(0,T)E_P[\frac{dQ}{dP} H_T]" />

Another name for the risk-neutral measure is the equivalent martingale measure. If there is just one unique risk-neutral measure in the market, then there is a unique arbitrage-free price for each asset in the market. This is the fundamental theorem of arbitrage-free pricing.If there are more such measures, then in an interval of prices no arbitrage is possible. If no equivalent martingale measure exists, arbitrage opportunities do. 

Suppose our economy consists of 2 assets, a stock and a risk-free bond, and that we use Black-Scholes model. In the model the evolution of the stock price can be described by Geometric Brownian Motion:

<img src="https://tex.s2cms.ru/svg/dS_t%20%3D%20%5Calpha%20S_t%20dt%20%2B%20%5Csigma%20S_t%20dW_t" alt="dS_t = \alpha S_t dt + \sigma S_t dW_t" />

where <img src="https://tex.s2cms.ru/svg/W_t" alt="W_t" /> is a standard Brownian motion with respect to the physical measure. If we define

<img src="https://tex.s2cms.ru/svg/%5Ctilde%7BW%7D_t%20%3D%20W_t%20%2B%20%5Cfrac%7B%5Calpha-r%7D%7B%5Csigma%7Dt" alt="\tilde{W}_t = W_t + \frac{\alpha-r}{\sigma}t" />

Girsanov's theorem states that there exists a measure <img src="https://tex.s2cms.ru/svg/Q" alt="Q" /> under which <img src="https://tex.s2cms.ru/svg/%7B%5Cdisplaystyle%20%7B%5Ctilde%20%7BW%7D%7D_%7Bt%7D%7D" alt="{\displaystyle {\tilde {W}}_{t}}" /> is a Brownian motion. Put this back in the original equation:

<img src="https://tex.s2cms.ru/svg/dS_t%20%3D%20r%20S_t%20dt%20%2B%20%5Csigma%20S_t%20d%5Ctilde%7BW%7D_t" alt="dS_t = r S_t dt + \sigma S_t d\tilde{W}_t" />

Then, the discounted stock price <img src="https://tex.s2cms.ru/svg/%5Ctilde%7BS%7D_t" alt="\tilde{S}_t" />is a <img src="https://tex.s2cms.ru/svg/Q" alt="Q" />-martingale.

<img src="https://tex.s2cms.ru/svg/d%20%5Ctilde%7BS%7D_t%20%3D%20%5Csigma%20%5Ctilde%7BS%7D_t%20d%5Ctilde%7BW%7D_t" alt="d \tilde{S}_t = \sigma \tilde{S}_t d\tilde{W}_t" />

Note that risk neutral measure is powerful because you don't need to replicate a portfolio in order to be arbitrage-free compared to risk-free bond. Under such measure, expected value(first moment) of securities are equal to rolling it into a deposit account (inflate it to the maturity date at the riskless rate of interest). 
