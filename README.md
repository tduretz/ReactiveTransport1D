# Reactive Transport 1D

## Introduction

The goal of this project is to model the infiltration of ionic fluids into rocks by coupling the speciation calculation, reactivity with minerals in the infiltrated rock and rock mechanics to study the evolution of permeability and to compare with traditional models without complex fluids.
At present, the project is composed of test for brucite dissolution with comparison with Eugster and Baumgartner 1987 and computations from EQ3 using the DEW model at 1 GPa and 400 °C.

## Installation

1. Clone/Download the repository
2. Run Julia from within the folder 
3. In Julia's REPL switch to package mode: type `]`
4. Activate the environement: type `activate .`
5. Install all necessary dependencies: type `instantiate`

## Test against previous Brucite dissolution results

The results from this code are in good agreement with both Brucite dissolutions computations from Eugster and Baumgartner 1987 and EQ3 computations at 1 GPa and 400 °C. In both case the main species for Mg is the same and the molalities are within a few percents to the results. It is worth noting that EQ3 uses activity coefficients for aqueous species (hin the test the B-dot formulation was used), while this code assumes activity coefficients are equal to 1. This will be changed in the near future.

## How to modify the computation

To test this speciation calculator, you can change the matrix in the data folder and change the log Ks, line 59 (b variable). You can also change the chlorinity (Cl<sup>tot</sup> variable) line 50.
