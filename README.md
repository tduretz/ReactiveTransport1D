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

The code here follows the method described in Eugster and Baumgartner 1987, which solves a set of non linear equations made of the charge balance, the tot chlorine dissolved into the fluid, the dissolution of brucite reaction and the dissociation reactions for all aqueous species.
For the comparison with EQ3, reactions are written slightly differently to the E&B 1987 paper.

### Reactions:
1) Brucite dissolution:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Brucite + 2 H<sup>+</sup> = Mg<sup>2+</sup> + H<sub>2</sub>O
2) MgCl<sup>+</sup> dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MgCl<sup>+</sup> = Mg<sup>2+</sup> + Cl<sup>-</sup>
3) HCl dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;HCl = H<sup>+</sup> + Cl<sup>-</sup>
4) H<sub>2</sub>O dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;H<sub>2</sub>O = H<sup>+</sup> + OH<sup>-</sup>
5) Mg(OH)<sub>2</sub> dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Mg(OH)<sub>2</sub> + 2 H<sup>+</sup> = 2 H<sub>2</sub>O + Mg<sup>2+</sup>
6) Mg(OH)<sup>+</sup> dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Mg(OH)<sup>+</sup> + H<sup>+</sup> = H<sub>2</sub>O + Mg<sup>2+</sup>

### Equations:
Charge balance:&emsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;m<sub>MgCl<sup>+</sup></sub> + 2 m<sub>Mg<sup>2+</sup></sub> - m<sub>Cl<sup>-</sup></sub> + m<sub>H<sup>+</sup></sub> - m<sub>OH<sup>-</sup></sub> + m<sub>MgOH<sup>+</sup></sub> = 0<br/>
Chlorinity:&emsp;&emsp;&emsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;m<sub>MgCl<sup>+</sup></sub> + m<sub>Cl<sup>-</sup></sub> + m<sub>HCl</sub> = 0.01<br/>
Br dissolution:&emsp;&emsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;log m<sub>Mg<sup>2+</sup></sub> - 2 log m<sub>H<sup>+</sup></sub> - log K<sub>Br</sub> = 0<br/>
MgCl<sup>+</sup> dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;log m<sub>Mg<sup>2+</sup></sub> + log m<sub>Cl<sup>-</sup></sub> - log m<sub>MgCl<sup>+</sup></sub> - log K<sub>MgCl<sup>+</sup></sub> = 0<br/>
HCl dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;log m<sub>H<sup>+</sup></sub> + log m<sub>Cl<sup>-</sup></sub> - log m<sub>HCl</sub> - log K<sub>HCl</sub> = 0<br/>
H<sub>2</sub>O dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;log m<sub>H<sup>+</sup></sub> + log m<sub>OH<sup>-</sup></sub> - log K<sub>w</sub> = 0<br/>
Mg(OH)<sub>2</sub> dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;log m<sub>Mg<sup>2+</sup></sub> - log m<sub>Mg(OH)<sub>2</sub></sub> - 2 log m<sub>H<sup>+</sup></sub> - log K<sub>Mg(OH)<sub>2</sub></sub> = 0<br/>
Mg(OH)<sup>+</sup> dissociation:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;log m<sub>Mg<sup>2+</sup></sub> - log m<sub>MgOH<sup>+</sup></sub> + log m<sub>H<sup>+</sup></sub> - log K<sub>MgOH<sup>+</sup></sub> = 0<br/>

## How to modify the computation

To test this speciation calculator, you can change the matrix in the data folder and change the log Ks, line 59 (b variable). You can also change the chlorinity (Cl<sup>tot</sup> variable) line 50.

## Authors
Thibault Duretz
Guillaume Siron
