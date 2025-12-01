# Reproducibility Materials for Generalized Entropy Calibration (GEC)

This repository contains reproducibility materials for the manuscript titled  
**"Generalized Entropy Calibration for the Analysis of Voluntary Survey Data"**  
submitted to *Biometrics*.

## Folder Structure

- **code/**
  - `Simulation.R`: Code for running the simulation study.
  - `RealData.R`: Code for analyzing the real data example. The real dataset is not publicly available due to confidentiality restrictions.
  - `Simulation_extra.R`: Code for running the additional simulation in the supplementary material.  
  - `Archive/`: Additional scripts and legacy code.

- **data/**
  - `nhis.Rdata`: Required dataset for running `Simulation.R`.

- **output/**
  - `CI.png`: Output figure showing confidence intervals from the real data analysis.

## Instructions

To reproduce the results:
1. Install and load the required packages listed in the beginning of `Simulation.R`.
2. Run `Simulation.R` to generate simulation results.

For questions, please contact yhkwon@kma.ac.kr.
