# Reproducibility Materials for Generalized Entropy Calibration (GEC)

This repository contains reproducibility materials for the manuscript titled  
**"Generalized Entropy Calibration for the Analysis of Voluntary Survey Data"**  
submitted to *Biometrics*.

## Folder Structure

- **code/**
  - `RealData.R`: Code for analyzing the real data example.
  - `Simulation.R`: Code for running the simulation study.
  - `Archive/`: Additional scripts and legacy code.

- **data/**
  - `nhis.Rdata`: Required dataset for running `Simulation.R`.

- **output/**
  - `CI.png`: Output figure showing confidence intervals from the simulation study.

## Instructions

To reproduce the results:
1. Load the required packages listed in the beginning of each script.
2. Run `Simulation.R` to generate simulation results and figures.
3. Run `RealData.R` to reproduce the real data analysis.

For questions, please contact yhkwon@kma.ac.kr.
