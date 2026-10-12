# Semi-Device-Independent Quantum Randomness Expansion (SDI-QRE)

This repository contains MATLAB code for calculating finite-statistics randomness-expansion (RE) rates in semi-device-independent (SDI) protocols based on quantum contextuality. The analysis uses semidefinite programming (SDP) to bound the conditional entropy and the entropy accumulation theorem (EAT) to estimate the net randomness expansion rate.

## Repository structure

```text
SDI-QRE/
├── KCBS_Randomness_Expansion/
│   └── MATLAB code and README for the KCBS (odd-cycle) scenario
├── CHSH_Randomness_Expansion/
│   └── MATLAB code and README for the CHSH (N = 4 even-cycle) scenario
└── README.md
```

### `KCBS_Randomness_Expansion/`
Contains the implementation for the KCBS scenario, corresponding to an odd-cycle noncontextuality inequality. The code constructs the KCBS projectors, solves the entropy SDP, obtains the dual variables used in the min-tradeoff function, and calculates the finite-size randomness expansion rate.

### `CHSH_Randomness_Expansion/`
Contains the implementation for the CHSH scenario, treated in the paper as an even-cycle contextuality scenario implemented by a single localized system. The code constructs the even-cycle projectors and uses a real-valued SDP formulation to handle the complex-valued quantum operators. It then calculates the finite-size randomness expansion rate using EAT.

## General workflow

For either scenario, the calculation follows these main steps:

1. **Build the measurements:** construct the projectors for the chosen cycle scenario.
2. **Bound the entropy:** solve an SDP for the conditional von Neumann entropy bound and obtain the associated dual variables.
3. **Calculate the RE rate:** use the entropy bound, spot-checking probability, number of rounds, confidence parameters, and EAT correction terms to compute the net randomness expansion rate in bits per round.

## Requirements

- MATLAB
- [YALMIP](https://yalmip.github.io/) for formulating the SDP
- [MOSEK](https://www.mosek.com/) as the SDP solver

Make sure YALMIP and MOSEK are installed and configured in MATLAB before running the codes.

## How to use

Open the relevant scenario folder and follow the instructions in its `README.md` for the required filenames, input parameters, and example calls. The KCBS and CHSH implementations use different measurement constructions and should be run using the code provided in their respective folders.

The calculated rate depends on the chosen physical and protocol parameters, including the observed inequality value, confidence interval, spot-checking probability, total number of rounds, quadrature order, and EAT security parameters.

## Reference

The code accompanies the work:

**Dhananjoy Ghosh, Nandana J, and Jaskaran Singh, “Local finite rate semi device independent randomness expansion from quantum contextuality.”**

Please refer to the paper for the protocol, assumptions, security analysis, and definitions of the quantities implemented in the code.

## Citation

If you use this code in your research, please cite the associated paper.


## License

This project is licensed under the MIT License. See the LICENSE file for details.









