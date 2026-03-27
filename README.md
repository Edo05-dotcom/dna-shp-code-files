# DNA Supercoiling Code Files

This repository contains the code used for my Senior Honours Project on the **stochastic model of supercoiling-dependent transcription**.

The project studies how DNA supercoiling can drive transcription from a **Poissonian regime** to a **supercoiling-regulated regime**, and uses quantities such as:

- Shannon entropy
- Permutation entropy
- Transcription rate / mean-field fitting

## What is in this repository?

The repo includes code used to:

- run the **Fortran 77 stochastic simulation** (also computes Shannon entropy over 1 simulation)
- generate heat map plots and Gaussian Process fitting with RBF and Matern kernels
- compute permutation entropy
- fit analytical and semi-analytical mean-field approximations to transcription rate

## Main files

Depending on the folder / file names, this repository includes:

- **Fortran code** for the travelling-polymerase simulation
- **Python scripts** for post-processing, fitting, and plotting
## How to use

A typical workflow is:

1. Run the **Fortran simulation** with the desired parameters and gene configuration.
2. Save the output data files.
3. Use the **Python** scripts to:
   - compute permutation entropy
   - normalise transcription rates
   - fit mean-field curves
   - plot heat maps and fit a Gaussian Process with different kernels

## Notes

This repository was created mainly to store the code used in the report, so it is not packaged as a polished software project.  
The scripts are intended to reproduce the figures and analysis in the dissertation rather than serve as a general-purpose toolkit.

## Project context

This code supports the project:

**Stochastic Model of Supercoiling-Dependent Transcription**

which investigates how supercoiling-regulated feedback affects transcription dynamics and whether the crossover to the supercoiling-regulated regime is smooth.

## Author

**Edoardo John Harris**  
School of Physics and Astronomy  
The University of Edinburgh
