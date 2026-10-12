# Local Finite-Rate Semi-Device-Independent Randomness Expansion from Quantum Contextuality

**Authors:** Dhananjoy Ghosh, Nandana J, and Jaskaran Singh

**Affiliation:** Centre for Quantum Science and Technology (CQST), IIT Mandi 

**arXiv:** [arxiv link]  

---


# Overview

This code accompanies the paper **“Local finite rate semi device independent randomness expansion from quantum contextuality”**.

The paper studies semi-device-independent randomness expansion (SDI-RE) using contextual correlations and finite-statistics security analysis based on the entropy accumulation theorem (EAT). The MATLAB routines supplied here focus on the entropy SDP and randomness-expansion-rate calculation for odd-cycle contextuality scenarios, including the KCBS case \(N=5\).

## What the code does

The implementation has two main computational stages:

1. **Entropy SDP:** constructs the qutrit projectors for an odd \(N\)-cycle and solves an SDP using YALMIP and MOSEK. It returns an entropy bound and dual variables associated with the normalization and observed-violation constraints.
2. **Finite-rate calculation:** uses the dual variables to construct the min-tradeoff function, evaluate the EAT correction terms, estimate the smooth min-entropy, subtract the input-randomness cost, and return the net randomness expansion rate in bits per round.

The rate calculated by the supplied `RE_rate` routine follows the paper's reported rate definition, which omits the fixed randomness-extractor entropy loss from the plotted rate.

## Requirements

- MATLAB
- [YALMIP](https://yalmip.github.io/) for modelling and solving the semidefinite program
- [MOSEK](https://www.mosek.com/) as the SDP solver, installed and configured for use through YALMIP

The code uses MATLAB's built-in linear algebra routines for the Gauss–Radau quadrature calculation.

## Functions and inputs

| Function | Purpose | Main inputs / outputs |
|---|---|---|
| `hae_kcbs(betaQ, delta, m, N)` | Solves the entropy SDP. | Returns `H`, `lambda0`, and `lambda1`. Here `betaQ` is the observed normalized inequality value, `delta` is the confidence-interval offset used by the code, `m` is the Gauss–Radau quadrature order, and `N` is the number of measurements in the odd cycle. |
| `RE_rate(betaQ, delta, q, N, n, lambda, m, epsEA, epss, gamma)` | Calculates the net randomness expansion rate using the EAT bound. | `q` is the spot-checking probability; `n` is the number of rounds; `lambda = [lambda0; lambda1]` contains the dual variables from `hae_kcbs`; `epsEA` and `epss` are security parameters; `gamma` is the EAT-order parameter. Returns `r`, in bits per round. |
| `kcbs_projectors(N)` | Builds the qutrit rank-one projectors for an odd-cycle realization. | `N` must be odd. |
| `gauss_radau(m)` | Calculates Gauss–Radau quadrature weights and nodes. | Returns `w` and `t`. |
| `binary_entropy(q)` | Calculates the binary entropy \(h(q)\). | Used to calculate the input-randomness cost. |


## Example: KCBS case (N = 5)

After placing the callable functions in MATLAB files as described above, run the following example from the MATLAB Command Window:

```matlab
% KCBS / odd-cycle parameters
N      = 5;
betaQ  = 0.4460;
delta  = 1e-3;
m      = 8;

% Stage 1: solve the entropy SDP and obtain its dual variables
[H, lambda0, lambda1] = hae_kcbs(betaQ, delta, m, N);
lambda = [lambda0; lambda1];

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

These are illustrative parameter choices for an \(N=5\) calculation, not a guarantee of a particular numerical rate. The result depends on the observed violation, confidence parameter, spot-checking probability, number of rounds, EAT parameters, and SDP solver precision.

## Notes on parameters

- `betaEff = betaQ - delta` is the effective inequality value used in the SDP and rate calculation.
- The code's input-randomness cost is
  \[
  \ell_{\mathrm{in}} = n\bigl[h(q) + q\log_2 N\bigr] + 6,
  \]
  where \(h(q)\) is the binary entropy.
- The code implements an odd-cycle qutrit projector construction. `kcbs_projectors(N)` explicitly rejects even values of `N`; therefore this file should not be used directly to reproduce the paper's separate \(N=4\) CHSH plots.
- Use `0 < q < 1` for the rate calculation, since the variance bound in the code contains division by `q`.
- The output is numerical. Small differences may occur with different MATLAB, YALMIP, or MOSEK versions and solver tolerances.

## Reference

Dhananjoy Ghosh, Nandana J, and Jaskaran Singh, **“Local finite rate semi device independent randomness expansion from quantum contextuality.”** 

