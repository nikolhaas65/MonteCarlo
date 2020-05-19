# Read Me

Author: Jerry Xia

Date: 2018/06/19

*Note: The advanced Marckdown features such as math expression may not be compatible in GitHub, please see README.pdf instead if you want more details*

## Implementation

Please feel free to see the Monte Carlo engine: MonteCarlo.py

### Classification

* regular Monte Carlo simulation
* optimal hedged Monte Carlo simulation
* delta-based Monte Carlo simulation
* Monte Carlo with antithetic variates
* Least square method of Longstaff and Schiwatz (LSM)
* Hedeged Least Square method (HLSM)

### Underlying Process

* geometric Brownian motion
* CIR model
* Heston model

### Boundary Scheme (CIR model)

* absorption
* reflection
* Higham and Mao
* partial truncation
* full truncation



# Optimal Hedged Monte Carlo

Model Inventors: Marc Potters, Jean-Philippe Bouchaud, Dragan Sestovic



* ## 1 Introduction

  This is a Python Notebook about variance reduction Monte Carlo simulations. In this script, I implemented the following variance reduction methods as well as their antithetic variates' version:

  * regular Monte Carlo
  * Monte Carlo with delta-based control variates
  * optimal hedged Monte Carlo

  Due to the significance and robustness, I mainly focus on the optimal hedged Monte Carlo (OHMC) in option pricing. We invoke this method to price European options and make comparison with other methods.

  ### 1.1 Facts
  * The option price is not simply the average value of the discounted future pay-off over the objective (or historical) probability distribution
  * The requirement of absence of arbitrage opportunities is equivalent to the existence of "risk-neutral measure", such that the price is indeed its average discounted future pay-off.
  * Risk in option trading cannot be eliminated

  ### 1.2 Objective
  * It would be satisfactory to have an option theory where the objective stochastic process of the underlying is used to calculate the option price, the hedge strategy and the *residual risk*.

  ### 1.3 Advantages
  * It is a versatile methods to price complicated path-dependent options.
  * Considerable variance reduction scheme for Monte Carlo
  * It provide not only a numerical estimate of the option price, but also of the optimal hedge strategy and of the residual risk.
  * This method does not rely on the notion of risk-neutral measure, and can be used to any model of the true dynamics of the underlying

## 2 Underlying dynamics

### Black-Scholes Model
<img src="https://tex.s2cms.ru/svg/dS%20%3D%20%5Cmu%20S%20dt%20%2B%20%5Csigma%20S%20dW_t" alt="dS = \mu S dt + \sigma S dW_t" />
<img src="https://tex.s2cms.ru/svg/log%20S_%7Bt%2B1%7D%20%3D%20log%20S_t%20%2B(%5Cmu%20-%20%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D)%5CDelta%20t%20%2B%20%5Csigma%20%5Csqrt%7B%5CDelta%20t%7D%20%5Cepsilon" alt="log S_{t+1} = log S_t +(\mu - \frac{\sigma^2}{2})\Delta t + \sigma \sqrt{\Delta t} \epsilon" />

where

<img src="https://tex.s2cms.ru/svg/%5Cepsilon%20%5Csim%20N(0%2C1)" alt="\epsilon \sim N(0,1)" />

In risk neutral measure, <img src="https://tex.s2cms.ru/svg/%5Cmu%20%3D%20r%20-%20q" alt="\mu = r - q" />. 

### Heston Model
The basic Heston model assumes that $S_t$, the price of the asset, is determined by a stochastic process:

<img src="https://tex.s2cms.ru/svg/%0AdS_t%20%3D%20%5Cmu%20S_t%20dt%20%2B%20%5Csqrt%7Bv_t%7D%20S_t%20d%20W_t%5ES%5C%5C%0Adv_t%20%3D%20%5Ckappa%20(%5Ctheta%20-%20v_t)%20dt%20%2B%20%5Cxi%20%5Csqrt%7Bv_t%7D%20d%20W_t%5Ev%0A" alt="
dS_t = \mu S_t dt + \sqrt{v_t} S_t d W_t^S\\
dv_t = \kappa (\theta - v_t) dt + \xi \sqrt{v_t} d W_t^v
" />

where 

<img src="https://tex.s2cms.ru/svg/E%5BdW_t%5ES%2CdW_t%5Ev%5D%3D%5Crho%20dt" alt="E[dW_t^S,dW_t^v]=\rho dt" />

In risk neutral measure, <img src="https://tex.s2cms.ru/svg/%5Cmu%20%3D%20r%20-%20q" alt="\mu = r - q" />. 

## 3 Methodology

### 3.1 Simbol Definition
Option price always requires to work backward. That is because the option price is known exactly at the maturity. As with other schemes, we determine the option price step by step from the maturity <img src="https://tex.s2cms.ru/svg/t%3DK%5Ctau%3DT" alt="t=K\tau=T" /> to the present time <img src="https://tex.s2cms.ru/svg/t%3D0" alt="t=0" />. The unit of time being <img src="https://tex.s2cms.ru/svg/%5Ctau" alt="\tau" />, for example, one day. We simulate <img src="https://tex.s2cms.ru/svg/N" alt="N" /> trajectories. In trajectory <img src="https://tex.s2cms.ru/svg/i" alt="i" />, the price of the underlying asset at time <img src="https://tex.s2cms.ru/svg/k%5Ctau" alt="k\tau" /> is denoted as <img src="https://tex.s2cms.ru/svg/S_k%5E%7B(i)%7D" alt="S_k^{(i)}" />. The price of the derivative at time <img src="https://tex.s2cms.ru/svg/k%5Ctau" alt="k\tau" /> is denoted as <img src="https://tex.s2cms.ru/svg/C_k" alt="C_k" />, and the hedge function is <img src="https://tex.s2cms.ru/svg/H_k" alt="H_k" />. We define an optimal hedged portfolio as

<img src="https://tex.s2cms.ru/svg/W_k%5E%7B(i)%7D%20%3D%20C_k(S_k%5E%7B(i)%7D)%20%2B%20H_k(S_k%5E%7B(i)%7D)S_k%5E%7B(i)%7D" alt="W_k^{(i)} = C_k(S_k^{(i)}) + H_k(S_k^{(i)})S_k^{(i)}" />

The one-step change of our portfolio is

<img src="https://tex.s2cms.ru/svg/%5CDelta%20W_k%5E%7B(i)%7D%3D%20df(k%2Ck%2B1)%20C_%7Bk%2B1%7D(S_%7Bk%2B1%7D%5E%7B(i)%7D)%20-%20C_k(S_k%5E%7B(i)%7D)%20%2B%20H_k(S_%7Bk%7D%5E%7B(i)%7D)%20(df(k%2Ck%2B1)%20S_%7Bk%2B1%7D%5E%7B(i)%7D%20-%20S_%7Bk%7D%5E%7B(i)%7D)" alt="\Delta W_k^{(i)}= df(k,k+1) C_{k+1}(S_{k+1}^{(i)}) - C_k(S_k^{(i)}) + H_k(S_{k}^{(i)}) (df(k,k+1) S_{k+1}^{(i)} - S_{k}^{(i)})" />

Where <img src="https://tex.s2cms.ru/svg/df(k%2Ck%2B1)" alt="df(k,k+1)" /> is the discounted factor from time <img src="https://tex.s2cms.ru/svg/k%5Ctau" alt="k\tau" /> to <img src="https://tex.s2cms.ru/svg/(k%2B1)%20%5Ctau" alt="(k+1) \tau" />, <img src="https://tex.s2cms.ru/svg/df2(k%2Ck%2B1)" alt="df2(k,k+1)" /> is the discounted factor considering dividend <img src="https://tex.s2cms.ru/svg/e%5E%7B-(r-q)(t_%7Bk%2B1%7D-t_k)%7D" alt="e^{-(r-q)(t_{k+1}-t_k)}" />

### 3.2 Objective
The optimal hedged algorithm can be interpreted as the following optimal problem

<img src="https://tex.s2cms.ru/svg/%0A%5Cbegin%7Balign%7D%0A%5Cbegin%7Bsplit%7D%0A%5Cmbox%7Bminimize%7D%5Cquad%20%26%20%5Cquad%20Var%5B%5CDelta%20W_k%5D%5C%5C%0A%5Cmbox%7Bsubject%20to%7D%5Cquad%20%26%20%5Cquad%20E%5B%5CDelta%20W_k%5D%3D0%0A%5Cend%7Bsplit%7D%5Cnonumber%0A%5Cend%7Balign%7D%0A" alt="
\begin{align}
\begin{split}
\mbox{minimize}\quad &amp; \quad Var[\Delta W_k]\\
\mbox{subject to}\quad &amp; \quad E[\Delta W_k]=0
\end{split}\nonumber
\end{align}
" />

It means we should try to minimize the realized volatility of hedged portfolio while maintaining the expected value of portfolio unchanged.

### 3.3 Basis Functions
The original optimization is very difficult to solve. Thus we assume a set of basis function and solved it in such subspace. We use $$N_C$$ and $$N_H$$ to denote the number of basis functions for price and hedge.

<img src="https://tex.s2cms.ru/svg/%0A%5Cbegin%7Balign%7D%0A%5Cbegin%7Bsplit%7D%0AC_k(%5Ccdot)%20%26%3D%20%5Csum_%7Bi%3D0%7D%5E%7BN_C%7D%20a_%7Bk%2Ci%7D%20A_i(%5Ccdot)%5C%5C%0AH_k(%5Ccdot)%20%26%3D%20%5Csum_%7Bi%3D0%7D%5E%7BN_H%7D%20b_%7Bk%2Ci%7D%20B_i(%5Ccdot)%0A%5Cend%7Bsplit%7D%5Cnonumber%0A%5Cend%7Balign%7D%0A" alt="
\begin{align}
\begin{split}
C_k(\cdot) &amp;= \sum_{i=0}^{N_C} a_{k,i} A_i(\cdot)\\
H_k(\cdot) &amp;= \sum_{i=0}^{N_H} b_{k,i} B_i(\cdot)
\end{split}\nonumber
\end{align}
" />

The basis functions <img src="https://tex.s2cms.ru/svg/A_i" alt="A_i" /> and <img src="https://tex.s2cms.ru/svg/B_i" alt="B_i" /> are priori determined and need not to be identical. The coefficients <img src="https://tex.s2cms.ru/svg/a_i" alt="a_i" /> and <img src="https://tex.s2cms.ru/svg/b_i" alt="b_i" /> can be calibrated by solving the optimal problem.

### 3.4 Numerical Solution
<img src="https://tex.s2cms.ru/svg/%0A%5Cbegin%7Balign%7D%0A%5Cbegin%7Bsplit%7D%0A%5Cmbox%7Bminimize%7D%5Cquad%20%26%20%5Cquad%20%5Cfrac%7B1%7D%7BN%7D%20%5Csum_%7Bi%3D1%7D%5EN%20%5CDelta%20W_k%5E%7B(i)2%7D%5C%5C%0A%5Cmbox%7Bsubject%20to%7D%5Cquad%20%26%20%5Cquad%20%5Cfrac%7B1%7D%7BN%7D%20%5Csum_%7Bi%3D1%7D%5EN%20%5CDelta%20W_k%5E%7B(i)%7D%3D0%0A%5Cend%7Bsplit%7D%5Cnonumber%0A%5Cend%7Balign%7D%0A" alt="
\begin{align}
\begin{split}
\mbox{minimize}\quad &amp; \quad \frac{1}{N} \sum_{i=1}^N \Delta W_k^{(i)2}\\
\mbox{subject to}\quad &amp; \quad \frac{1}{N} \sum_{i=1}^N \Delta W_k^{(i)}=0
\end{split}\nonumber
\end{align}
" />

Denote the discounted forward underlying price change at time <img src="https://tex.s2cms.ru/svg/k%5Ctau" alt="k\tau" /> as

<img src="https://tex.s2cms.ru/svg/%5CDelta%20S_k%20%3D%20df2(k%2Ck%2B1)%20S_%7Bk%2B1%7D%20-%20S_k" alt="\Delta S_k = df2(k,k+1) S_{k+1} - S_k" />

Define
<img src="https://tex.s2cms.ru/svg/%0A%5Cbegin%7Balign%7D%0A%5Cbegin%7Bsplit%7D%0AQ_k%20%26%3D%20%5Cbegin%7Bbmatrix%7D%0A%20%20%20%20-A_%7Bk%2C1%7D(S_k%5E%7B(1)%7D)%20%26%20%5Ccdots%20%26%20-A_%7Bk%2CN_C%7D(S_k%5E%7B(1)%7D)%20%26%20B_%7Bk%2C1%7D(S_k%5E%7B(1)%7D)%5CDelta%20S_k%5E%7B(1)%7D%26%20%5Ccdots%20%20%26%20B_%7Bk%2CN_H%7D(S_k%5E%7B(1)%7D)%5CDelta%20S_k%5E%7B(1)%7D%0A%5C%5C%0A%20%20%20%20-A_%7Bk%2C1%7D(S_k%5E%7B(2)%7D)%20%26%20%5Ccdots%20%26%20-A_%7Bk%2CN_C%7D(S_k%5E%7B(2)%7D)%20%26%20B_%7Bk%2C1%7D(S_k%5E%7B(2)%7D)%5CDelta%20S_k%5E%7B(2)%7D%26%20%5Ccdots%20%20%26%20B_%7Bk%2CN_H%7D(S_k%5E%7B(1)%7D)%5CDelta%20S_k%5E%7B(2)%7D%20%5C%5C%0A%20%20%20%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%5C%5C%0A%20%20%20%20-A_%7Bk%2C1%7D(S_k%5E%7B(N)%7D)%20%26%20%5Ccdots%20%26%20-A_%7Bk%2CN_C%7D(S_k%5E%7B(N)%7D)%20%26%20B_%7Bk%2C1%7D(S_k%5E%7B(N)%7D)%5CDelta%20S_k%5E%7B(N)%7D%26%20%5Ccdots%20%20%26%20B_%7Bk%2CN_H%7D(S_k%5E%7B(N)%7D)%5CDelta%20S_k%5E%7B(N)%7D%0A%20%20%20%20%5Cend%7Bbmatrix%7D%5C%5C%5C%5C%0Ac_k%20%26%3D%20(a_%7Bk%2C1%7D%2C%20%5Ccdots%20a_%7Bk%2CN_C%7D%2C%20b_%7Bk%2C1%7D%2C%20%5Ccdots%2C%20b_%7Bk%2CN_H%7D)%5ET%5C%5C%5C%5C%0Av_%7Bk%7D%20%26%3D%20df(k%2Ck%2B1)%20C_%7Bk%2B1%7D(S_%7Bk%2B1%7D%5E%7B%7D)%0A%5Cend%7Bsplit%7D%5Cnonumber%0A%5Cend%7Balign%7D%0A" alt="
\begin{align}
\begin{split}
Q_k &amp;= \begin{bmatrix}
    -A_{k,1}(S_k^{(1)}) &amp; \cdots &amp; -A_{k,N_C}(S_k^{(1)}) &amp; B_{k,1}(S_k^{(1)})\Delta S_k^{(1)}&amp; \cdots  &amp; B_{k,N_H}(S_k^{(1)})\Delta S_k^{(1)}
\\
    -A_{k,1}(S_k^{(2)}) &amp; \cdots &amp; -A_{k,N_C}(S_k^{(2)}) &amp; B_{k,1}(S_k^{(2)})\Delta S_k^{(2)}&amp; \cdots  &amp; B_{k,N_H}(S_k^{(1)})\Delta S_k^{(2)} \\
    \vdots &amp; \vdots &amp; \vdots &amp; \vdots &amp; \vdots &amp; \vdots\\
    -A_{k,1}(S_k^{(N)}) &amp; \cdots &amp; -A_{k,N_C}(S_k^{(N)}) &amp; B_{k,1}(S_k^{(N)})\Delta S_k^{(N)}&amp; \cdots  &amp; B_{k,N_H}(S_k^{(N)})\Delta S_k^{(N)}
    \end{bmatrix}\\\\
c_k &amp;= (a_{k,1}, \cdots a_{k,N_C}, b_{k,1}, \cdots, b_{k,N_H})^T\\\\
v_{k} &amp;= df(k,k+1) C_{k+1}(S_{k+1}^{})
\end{split}\nonumber
\end{align}
" />

As for <img src="https://tex.s2cms.ru/svg/v_k" alt="v_k" />, note that we know the exact value at maturity, which means there is no need to approximate price in terms of basis functions, that is

<img src="https://tex.s2cms.ru/svg/%0A%5Cbegin%7Balign%7D%0Av_k%20%3D%20%5Cbegin%7Bcases%7D%0Adf(N-1%2CN)%5C%20payoff(S_N)%2C%5Cquad%20%26%20k%3DN-1%5C%5C%0Adf(k%2Ck%2B1)%5C%20%5Csum_%7Bi%3D1%7D%5E%7BN_C%7D%20a_%7Bk%2B1%2Ci%7D%20A_i(S_%7Bk%2B1%7D)%2C%20%5Cquad%20%26%200%3Ck%3CN-1%5C%5C%0Adf(0%2C1)%5C%20C_1(S_1)%2C%20%5Cquad%20%26%20k%3D0%0A%5Cend%7Bcases%7D%5Cnonumber%0A%5Cend%7Balign%7D%0A" alt="
\begin{align}
v_k = \begin{cases}
df(N-1,N)\ payoff(S_N),\quad &amp; k=N-1\\
df(k,k+1)\ \sum_{i=1}^{N_C} a_{k+1,i} A_i(S_{k+1}), \quad &amp; 0&lt;k&lt;N-1\\
df(0,1)\ C_1(S_1), \quad &amp; k=0
\end{cases}\nonumber
\end{align}
" />

Then, the optimization problem can be expressed as

<img src="https://tex.s2cms.ru/svg/%0A%5Cbegin%7Balign%7D%0A%5Cbegin%7Bsplit%7D%0A%5Carg%5Cmin_%7Bc_k%7D%5Cquad%20%26%20%5Cquad%20(v_%7Bk%7D%20%2B%20Q_k%20c_k)%5ET%20(v_%7Bk%7D%20%2B%20Q_k%20c_k)%5C%5C%0A%5Cmbox%7Bsubject%20to%7D%5Cquad%20%26%20%5Cquad%201_%7B%5BN%5Ctimes1%5D%7D%5ET%20(v_%7Bk%7D%20%20%2B%20Q_k%20c_k)%3D0%0A%5Cend%7Bsplit%7D%5Cnonumber%0A%5Cend%7Balign%7D%0A" alt="
\begin{align}
\begin{split}
\arg\min_{c_k}\quad &amp; \quad (v_{k} + Q_k c_k)^T (v_{k} + Q_k c_k)\\
\mbox{subject to}\quad &amp; \quad 1_{[N\times1]}^T (v_{k}  + Q_k c_k)=0
\end{split}\nonumber
\end{align}
" />

In step k, since we already know the information ($v_{k}$) in step k+1. By canceling the constant term, the optimal problem can be simplified as the following 

<img src="https://tex.s2cms.ru/svg/%0A%5Cbegin%7Balign%7D%0A%5Cbegin%7Bsplit%7D%0A%5Carg%5Cmin_%7Bc_k%7D%5Cquad%20%26%20%5Cquad%202%20v_%7Bk%7D%5ET%20Q_k%20c_k%20%2B%20c_k%5ET%20Q_k%5ET%20Q_k%20c_k%5C%5C%0A%5Cmbox%7Bsubject%20to%7D%5Cquad%20%26%20%5Cquad%201_%7B%5BN%5Ctimes1%5D%7D%5ET%20v_%7Bk%7D%20%20%2B%201_%7B%5BN%5Ctimes1%5D%7D%5ET%20Q_k%20c_k%3D0%0A%5Cend%7Bsplit%7D%5Cnonumber%0A%5Cend%7Balign%7D%0A" alt="
\begin{align}
\begin{split}
\arg\min_{c_k}\quad &amp; \quad 2 v_{k}^T Q_k c_k + c_k^T Q_k^T Q_k c_k\\
\mbox{subject to}\quad &amp; \quad 1_{[N\times1]}^T v_{k}  + 1_{[N\times1]}^T Q_k c_k=0
\end{split}\nonumber
\end{align}
" />

### 3.5 Convex Optimization Problem

Let us first review the standard form of linear constrained quadratic programming problem:

<img src="https://tex.s2cms.ru/svg/%5Cmin_%7Bx%7D%20%5Cquad%20%20%5Cfrac%7B1%7D%7B2%7D%20x%5ET%20P%20x%20%2B%20q%5ET%20x" alt="\min_{x} \quad  \frac{1}{2} x^T P x + q^T x" />

<img src="https://tex.s2cms.ru/svg/%5Cmbox%7Bsubject%20to%7D%20%5Cquad%20G%20x%20%5Cpreceq%20h" alt="\mbox{subject to} \quad G x \preceq h" />

<img src="https://tex.s2cms.ru/svg/A%20x%20%3D%20b" alt="A x = b" />

Note that <img src="https://tex.s2cms.ru/svg/x%5ET" alt="x^T" /> means the transpose of vector x, and <img src="https://tex.s2cms.ru/svg/G%20x%20%5Cpreceq%20h" alt="G x \preceq h" /> denotes the inequality is taken element-wise over the vectors <img src="https://tex.s2cms.ru/svg/G%20x" alt="G x" /> and <img src="https://tex.s2cms.ru/svg/h" alt="h" />. The objective function is convex if and only if the matrix <img src="https://tex.s2cms.ru/svg/P" alt="P" /> is positive-semidefinite(Hermitian matrix all of whose eigenvalues are nonnegative), which is the realm we concern with.

Recall that the constrained optimization problem:

<img src="https://tex.s2cms.ru/svg/%5Carg%5Cmin%7Bc_k%7D%5Cquad%20%20%5Cquad%20%20v%7Bk%7D%5ET%20Q_k%20c_k%20%2B%20%5Cfrac%7B1%7D%7B2%7Dc_k%5ET%20Q_k%5ET%20Q_k%20c_k" alt="\arg\min{c_k}\quad  \quad  v{k}^T Q_k c_k + \frac{1}{2}c_k^T Q_k^T Q_k c_k" />

<img src="https://tex.s2cms.ru/svg/%5Cmbox%7Bsubject%20to%7D%5Cquad%20%20%5Cquad%201%7B%5BN%5Ctimes1%5D%7D%5ET%20v%7Bk%7D%20%20%2B%201_%7B%5BN%5Ctimes1%5D%7D%5ET%20Q_k%20c_k%3D0" alt="\mbox{subject to}\quad  \quad 1{[N\times1]}^T v{k}  + 1_{[N\times1]}^T Q_k c_k=0" />

Correspondingly, we make the connection by letting

<img src="https://tex.s2cms.ru/svg/x%20%3D%20c_k" alt="x = c_k" />

<img src="https://tex.s2cms.ru/svg/P%20%3D%20Q_k%5ET%20Q_k" alt="P = Q_k^T Q_k" />

<img src="https://tex.s2cms.ru/svg/q%20%3D%20Q_k%5ET%20v_k" alt="q = Q_k^T v_k" />

<img src="https://tex.s2cms.ru/svg/A%20%3D%201_%7B%5BN%5Ctimes1%5D%7D%5ET%20Q_k" alt="A = 1_{[N\times1]}^T Q_k" />

<img src="https://tex.s2cms.ru/svg/b%20%3D%20-1%7B%5BN%5Ctimes1%5D%7D%5ET%20v%7Bk%7D" alt="b = -1{[N\times1]}^T v{k}" />

The hard work is almost over right now. As you would always find, formulating the problem is usually the hard step. Invoking a solver is straightforward.

Note that when $k=0$, the degree of freedom of the quadratic problem decreases to 2. Because here the only concerns are price and hedge at time zero (we don't need to project them into a high dimension space). Let <img src="https://tex.s2cms.ru/svg/x%3D%5BC_0%2C%20H_0%5D%5ET" alt="x=[C_0, H_0]^T" />

<img src="https://tex.s2cms.ru/svg/Q_0%20%3D%20%5Cbegin%7Bbmatrix%7D%20-1%20%26%20%5CDelta%20S_0%5E%7B(1)%7D%5C%5C%0A%20%20%20%20%5Cvdots%20%26%20%5Cvdots%5C%5C%0A%20%20%20%20-1%20%26%20%5CDelta%20S_0%5E%7B(N)%7D%0A%20%20%20%20%5Cend%7Bbmatrix%7D" alt="Q_0 = \begin{bmatrix} -1 &amp; \Delta S_0^{(1)}\\
    \vdots &amp; \vdots\\
    -1 &amp; \Delta S_0^{(N)}
    \end{bmatrix}" />

<img src="https://tex.s2cms.ru/svg/P%20%3D%20Q_0%5ET%20Q_0" alt="P = Q_0^T Q_0" />

<img src="https://tex.s2cms.ru/svg/q%20%3D%20Q_0%5ET%20v_0" alt="q = Q_0^T v_0" />

<img src="https://tex.s2cms.ru/svg/A%20%3D%201_%7B%5BN%20%5Ctimes%201%5D%7D%5ET%20Q_0" alt="A = 1_{[N \times 1]}^T Q_0" />

<img src="https://tex.s2cms.ru/svg/b%20%3D%20-1_%7B%5BN%20%5Ctimes%201%5D%7D%5ET%20v_0" alt="b = -1_{[N \times 1]}^T v_0" />

## 4 Variance reduction and other methods
The rate of convergence of the Monte Carlo simulation is <img src="https://tex.s2cms.ru/svg/O%5Cleft(%5Cmax%20%5Cleft(%20%5CDelta%20t%2C%20%5Cfrac%7B1%7D%7BN_x%7D%20%5Cright)%5Cright)" alt="O\left(\max \left( \Delta t, \frac{1}{N_x} \right)\right)" />. The variance reduction techniques are used to reduce the constant factor corresponding to the Monte Carlo approximation <img src="https://tex.s2cms.ru/svg/O%20%5Cleft(%5Cfrac%7B1%7D%7BN_x%7D%5Cright)" alt="O \left(\frac{1}{N_x}\right)" />. Some of the most used variance reduction techniques are:

* Control Variates
* Antithetic Variates
* Moment Matching

In this part we selected antithetic variates and delta-based control variates methods as a supplement to optimal hedged monte carlo simulation.

### 4.1 Antithetic variates
The main idea of this technique is to look at the asset equation that you aretrying to simulate:
<img src="https://tex.s2cms.ru/svg/d%20S_t%5E%7B(1)%7D%20%3D%20r%20S_t%5E%7B(1)%7D%20dt%20%2B%20%5Csigma%20S_t%5E%7B(1)%7D%20d%20W_t" alt="d S_t^{(1)} = r S_t^{(1)} dt + \sigma S_t^{(1)} d W_t" />
and recognize that sinceztis a standard Brownian motion so will be−ztandthey will have the same exact distribution.  This means that the equation:
<img src="https://tex.s2cms.ru/svg/d%20S_t%5E%7B(2)%7D%20%3D%20r%20S_t%5E%7B(2)%7D%20dt%20-%20%5Csigma%20S_t%5E%7B(2)%7D%20d%20W_t" alt="d S_t^{(2)} = r S_t^{(2)} dt - \sigma S_t^{(2)} d W_t" />
will also generate paths of the same asset.
The variance depends on the sign of the covariance of <img src="https://tex.s2cms.ru/svg/%5Ctextit%7Bpayoff%7D(S_t%5E%7B(1)%7D)" alt="\textit{payoff}(S_t^{(1)})" /> and <img src="https://tex.s2cms.ru/svg/%5Ctextit%7Bpayoff%7D(S_t%5E%7B(2)%7D)" alt="\textit{payoff}(S_t^{(2)})" />. It can increase the eventual variance or decrease it, both case do arise. One sufficient condition to insure variance reduction is the monotony of the payoff function. Then, when using both in the calculation of the final Monte Carlo value the variance of the estimate will be reduced.

### 4.2 Delta-based control variates
Delta hedging can be summarized succinctly in the following way:  Suppose that at time <img src="https://tex.s2cms.ru/svg/t%3D%200" alt="t= 0" />, we receive <img src="https://tex.s2cms.ru/svg/C_0" alt="C_0" /> the price of an option that pays <img src="https://tex.s2cms.ru/svg/C_T" alt="C_T" /> at time T.  The price of this option at any time <img src="https://tex.s2cms.ru/svg/t" alt="t" /> is a function <img src="https://tex.s2cms.ru/svg/C(t%2CS)" alt="C(t,S)" />. Then, if we hold at any moment in time <img src="https://tex.s2cms.ru/svg/%5Cfrac%7B%5Cpartial%20C%7D%7B%5Cpartial%20S%7D(t%2CS)%20%3D%20%5Cfrac%7B%5Cpartial%20C_t%7D%7B%5Cpartial%20S%7D" alt="\frac{\partial C}{\partial S}(t,S) = \frac{\partial C_t}{\partial S}" /> units of stock, then we will be able to replicate the payout of this option <img src="https://tex.s2cms.ru/svg/C_T" alt="C_T" /> at time T. This is in theory since of course we cannot trade continuously. So in practice we perform a partial hedge where we only rebalance at some discrete moments in time say <img src="https://tex.s2cms.ru/svg/t_1%2Ct_2%2C%5Ccdots%2Ct_N" alt="t_1,t_2,\cdots,t_N" />. The replicating strategy can be expressed as follow:
<img src="https://tex.s2cms.ru/svg/W(t_i%2CS_i)%20%3D%20C(t_0%2CS_0)%20e%5E%7Br(t_i%20-%20t_0)%7D%20%2B%20%5Csum_%7Bj%3D0%7D%5E%7Bi%7D%20%5CDelta(t_j%2CS_j)%20(%20S_%7Bj%2B1%7D%20e%5E%7B-r(t_%7Bj%2B1%7D%20-%20t_j%20)%7D%20-%20S_%7Bj%7D)e%5E%7Br(t_i%20-%20t_j)%7D%20%3D%20C(t_i%2CS_i)" alt="W(t_i,S_i) = C(t_0,S_0) e^{r(t_i - t_0)} + \sum_{j=0}^{i} \Delta(t_j,S_j) ( S_{j+1} e^{-r(t_{j+1} - t_j )} - S_{j})e^{r(t_i - t_j)} = C(t_i,S_i)" />
which is similar to the strategy in the optimal hedged Monte Carlo simulation where the only difference is that in OHMC, we use option and delta hedging to replicate the cash flow and here we do the opposite operation. But when implementing the delta-based control variates, we should move the hedging term to the right hand side which make it identical to the OHMC strategy. Note that here we are assumed to know the delta hedging function. It explains a lot why OHMC can reduce the variance.

### 4.3 Optimal hedged Monte Carlo simulation
**In conclusion, OHMC is just a control variates method with an optimization on top and it is more practical because we do not have an analytical formula for the hedge sensitivity (i.e. delta, gamma, etc.)**



## 5 Add Hedging portfolio with the Least Square Monte Carlo (LSM)

In order to price American type options, we need to consider the problem of optimal exercise. LSM is a well-defined method to tackle this problem. In contrast, here we only utilize the information of exercise points along each simulation path using cross-sectional regression. Different from the original LSM, here we equipe basis functions to approximate price and hedge at each step similar to OHMC. And discuss independently at the inception.

This combination create a magic reaction. Now we can not only price the American options but also hedge it! Moreover, it's model independent, model parameters or construction, dimension doesn't matter at all! We use Black-Scholes and Heston model as examples. What only matters is the underlying price trials. With it, we can calculate the following stuffs.

* American options price
* American options Greeks
* American options optimal exercise boundary

Here, Bouchard and Warin concluded two main dynamic strategy in American options pricing, A1 and A2. Besides, I equiped them with a hedging strategy:

### 5.1 A1 strategy with optimal exercise time estimate

* Initialization: <img src="https://tex.s2cms.ru/svg/%5Ctau(t_J)%20%3D%20T" alt="\tau(t_J) = T" />
* Backward induction: <img src="https://tex.s2cms.ru/svg/%5Ctau(t_j)%20%3D%20t_j%20%5Cmathbf%7B1%7D_%7B%5C%7Bg(t_j)%5Cgeq%20C(t_j)%5C%7D%7D%20%2B%20%5Ctau(t_%7Bj%2B1%7D)%5Cmathbf%7B1%7D_%7B%5C%7BZ(t_j)%3CC(t_j)%5C%7D%7D" alt="\tau(t_j) = t_j \mathbf{1}_{\{g(t_j)\geq C(t_j)\}} + \tau(t_{j+1})\mathbf{1}_{\{Z(t_j)&lt;C(t_j)\}}" />
* Price estimator at 0: <img src="https://tex.s2cms.ru/svg/P_0%20%3D%20E%5Bg(%5Ctau(t_0)%2CX_%7B%5Ctau(t_0)%7D)%5D" alt="P_0 = E[g(\tau(t_0),X_{\tau(t_0)})]" />


### 5.2 A2 strategy with American values estimate

* Initialization: <img src="https://tex.s2cms.ru/svg/P_T%20%3D%20g(T%2CX_T)" alt="P_T = g(T,X_T)" />
* Backward induction: <img src="https://tex.s2cms.ru/svg/P_%7Bt_j%7D%20%3D%20max%5C%7Bg(t_j%2CX_%7Bt_j%7D)%2CE%5BP_%7Bt_%7Bj%2B1%7D%7D%5D%5C%7D" alt="P_{t_j} = max\{g(t_j,X_{t_j}),E[P_{t_{j+1}}]\}" />
* Price estimator at 0: <img src="https://tex.s2cms.ru/svg/P_0" alt="P_0" />

### 5.3 A2b strategy with optimal exercise time estimate and American values estimate

* Initialization: <img src="https://tex.s2cms.ru/svg/%5Ctau(t_J)%20%3D%20T" alt="\tau(t_J) = T" />
* Backward induction: 
    * <img src="https://tex.s2cms.ru/svg/%5Ctau(t_j)%20%3D%20t_j%20%5Cmathbf%7B1%7D_%7B%5C%7Bg(t_j)%5Cgeq%20C(t_j)%5C%7D%7D%20%2B%20%5Ctau(t_%7Bj%2B1%7D)%5Cmathbf%7B1%7D_%7B%5C%7BZ(t_j)%3CC(t_j)%5C%7D%7D" alt="\tau(t_j) = t_j \mathbf{1}_{\{g(t_j)\geq C(t_j)\}} + \tau(t_{j+1})\mathbf{1}_{\{Z(t_j)&lt;C(t_j)\}}" />
    * Price estimator at j: <img src="https://tex.s2cms.ru/svg/P_j%20%3D%20E%5Bg(%5Ctau(t_j)%2CX_%7B%5Ctau(t_j)%7D)%5D%24%20for%20%24j%3DJ%2CJ-1%2C%5Ccdots%2C1" alt="P_j = E[g(\tau(t_j),X_{\tau(t_j)})]$ for $j=J,J-1,\cdots,1" />
* Price estimator at 0 (one-step hedged MC): <img src="https://tex.s2cms.ru/svg/%7Barg%5C%2Cmin%7D_%7BP_0%2CH_0%7D%20E%5B(%5CDelta%20W_0)%5E2%5D" alt="{arg\,min}_{P_0,H_0} E[(\Delta W_0)^2]" />


### 5.5 Performance Test:


Black-Scholes model: HLSM-BlackScholes-American.ipynb
Heston model: HLSM-Heston-American.ipynb

In this document, we take Black-Scholes model as an example

#### Parameters:

    risk_free_rate = 0.06
    dividend = 0.0
    time_to_maturity = 1
    volatility = 0.3
    strike = 1.1
    stock_price = 1
    n_trials = 4000
    n_steps = 20
    func_list = [lambda x: x**0, lambda x: x] # basis for OHMC part
    option_type = 'p'
    
#### Results:

**American Options**

|Algorithm | Price  | Delta  |
|---|---|---|
| A1  | 0.1499 | N/A  | 
|  A2 |  0.1590 | 0.585  | 
|  A2b |  0.1500 | 0.491  | 


**European Options**

* BS Formula: 0.1401
* BS Binomial Tree: 0.1410
* Regular MC: 0.1453
* OHMC: 0.1426
