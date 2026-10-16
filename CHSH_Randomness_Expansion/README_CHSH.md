# Local Finite-Rate Semi-Device-Independent Randomness Expansion from Quantum Contextuality

**Authors:** Dhananjoy Ghosh, Nandana J, and Jaskaran Singh

**Affiliation:** Centre for Quantum Science and Technology (CQST), IIT Mandi 

**arXiv:** https://doi.org/10.48550/arXiv.2610.21765  

---


# Overview

This README documents the MATLAB implementation for the even-cycle contextuality scenario, including the CHSH case `N = 4`, associated with the paper **“Local finite rate semi device independent randomness expansion from quantum contextuality”**.

The code evaluates a finite-statistics randomness-expansion (RE) rate using an entropy semidefinite program (SDP) and the entropy accumulation theorem (EAT). In this implementation, the CHSH-type measurements are treated as measurements on one localized four-dimensional system, as described in the paper's semi-device-independent contextuality setting.

## What the code does

The calculation has two main stages:

1. **Entropy SDP:** `hae_chsh` builds the even-cycle measurement projectors and solves an SDP using YALMIP and MOSEK. It returns an entropy bound and the dual variables associated with the normalization and observed-cycle-value constraints.
2. **Finite-rate calculation:** `RE_rate` uses these dual variables to construct the min-tradeoff function, apply the EAT correction terms, estimate the smooth min-entropy, subtract the input-randomness cost, and return the net RE rate in bits per round.

The additional function `cycle_projectors_even` constructs the four-dimensional measurement projectors used by the even-cycle implementation.

## Requirements

- MATLAB
- [YALMIP](https://yalmip.github.io/) for modelling and solving the semidefinite program
- [MOSEK](https://www.mosek.com/) as the SDP solver, installed and configured for use through YALMIP

The quadrature and binary-entropy calculations use standard MATLAB operations.

## Functions and inputs

| Function | Purpose | Main inputs / outputs |
|---|---|---|
| `hae_chsh(betaQ, delta, m, N)` | Solves the entropy SDP. | Returns `H`, `lambda0`, and `lambda1`. `betaQ` is the observed normalized even-cycle inequality value, `delta` is the confidence-interval offset, `m` is the Gauss–Radau quadrature order, and `N` is the number of measurements in the cycle. The code uses a four-dimensional Hilbert space. |
| `RE_rate(betaQ, delta, q, N, n, lambda, m, epsEA, epss, gamma)` | Calculates the net randomness expansion rate using the EAT bound. | Returns `r` in bits per round. `q` is the spot-checking probability, `n` is the number of rounds, `lambda = [lambda0, lambda1]` contains the dual variables as a **row vector**, `epsEA` and `epss` are security parameters, and `gamma` is the EAT-order parameter. |
| `cycle_projectors_even(N)` | Constructs the even-cycle projectors in dimension four. | `N` is the number of measurements; this code is intended for even-cycle scenarios, particularly `N = 4`. |
| `gauss_radau(m)` | Calculates the Gauss–Radau quadrature weights and nodes. | Returns `w` and `t`. |
| `binary_entropy(q)` | Calculates the binary entropy \(h(q)\). | Used in the input-randomness cost. |


## Example: CHSH case (`N = 4`)

For the CHSH/even-cycle scenario, the paper gives the noncontextual bound as \(3/4\) and the maximum quantum value as \((2+\sqrt{2})/4 \approx 0.853553\). The values below use an illustrative observed value close to the one used in the paper; the numerical rate depends on all protocol parameters and on the solver output.

```matlab
% CHSH / even-cycle parameters
N      = 4;
betaQ  = 0.8526;
delta  = 1e-3;
m      = 8;

% Stage 1: solve the entropy SDP and obtain the dual variables
[H, lambda0, lambda1] = hae_chsh(betaQ, delta, m, N);

% IMPORTANT: RE_rate uses lambda * [1; betaEff], so lambda must be a row vector.
lambda = [lambda0, lambda1];

% Stage 2: calculate the finite-rate randomness expansion rate
n      = 1e12;
q      = 1e-2;
epsEA  = 1e-6;
epss   = 1e-6;
gamma  = 1/sqrt(n);

r = RE_rate(betaQ, delta, q, N, n, lambda, m, ...
            epsEA, epss, gamma);

fprintf('Entropy bound H = %.10f\n', H);
fprintf('Net RE rate    = %.10f bits/round\n', r);
```

These parameters are an example to show how the functions are connected; they do not guarantee a particular numerical result. The result depends on the observed cycle value, `delta`, `q`, `n`, the EAT parameters, the SDP solver, and numerical precision.

## Notes on the implementation

- The effective inequality value used in the code is `betaEff = betaQ - delta`.
- The input-randomness cost is calculated as
  \[
  \ell_{\mathrm{in}} = n\bigl[h(q) + q\log_2 N\bigr] + 6,
  \]
  where \(h(q)\) is the binary entropy.
- The variance bound in `RE_rate` includes division by `q`; use `0 < q < 1` for the intended spot-checking calculation.
- `lambda` must be a row vector, `[lambda0, lambda1]`, for the matrix multiplication in the provided `RE_rate` implementation.
- The routine calculates the rate expression implemented in the supplied code. It does not include the fixed extractor entropy-loss term in the plotted rate expression discussed in the paper.
- Results are numerical and may vary slightly with MATLAB, YALMIP, MOSEK, or solver settings.

## Reference

Dhananjoy Ghosh, Nandana J, and Jaskaran Singh, **“Local finite rate semi device independent randomness expansion from quantum contextuality.”** 

