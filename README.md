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

The results from this code are in good agreement with both Brucite dissolutions computations from Eugster and Baumgartner 1987 and EQ3 computations at 1 GPa and 400 °C. In both case the main species for Mg is the same and the molalities are within a few percents to the results.

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

### Mathematical formulation

The set of 8 non linear equations above can be rewritten in the form of the Jacobian matrix where:

$$
\begin{align}
\begin{bmatrix}
      0 & 1 & 2 & -1 & 1 & -1 & 0 & 1 \\
      1 & 1 & 0 & 1 & 0 & 0 & 0 & 0 \\
      0 & 0 & 1 & 0 & -2 & 0 & 0 & 0 \\
      0 & -1 & 1 & 1 & 0 & 0 & 0 & 0 \\
      -1 & 0 & 0 & 1 & 1 & 0 & 0 & 0 \\
      0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 \\
      0 & 0 & 1 & 0 & -2 & 0 & -1 & 0 \\
      0 & 0 & 1 & 0 & -1 & 0 & 0 & -1
\end{bmatrix} \cdot
\begin{bmatrix}
       m_{HCl}\\
       m_{MgCl^+}\\
       m_{Mg^{2+}}\\
       m_{Cl^-}\\
       m_{H^+}\\
       m_{OH^-}\\
       m_{Mg(OH)2}\\
       m_{MgOH^+}
\end{bmatrix} =
\begin{bmatrix}
       0\\
       0.01\\
       6.8466\\
       -1.0841\\
       -0.6078\\
       -8.1764\\
       6.6296\\
       4.9398
\end{bmatrix}
\end{align}
$$

## BackCalc algorithm
Another way to compute the speciation of a fluid was presented in Galvez et al. (2015). The basic idea involve switching the dissolution reactions that conencts the fluid to the rock with a hydrolysis reaction with the general form:

$$
\begin{align}
M_xO_y + 2yH^+ = xM^{2y/x+} + yH_2O
\end{align}
$$

With M being a cation. This allow to link the activity of a basis species (a the activity of H+) with the chemical potential of the related oxide and the Gibbs free energy of this same basis species using the following general equation:

$$
\begin{align}
log(a(M^{2y/x+})) - 2y/x·log(a(H^+)) = \dfrac{1}{RT}(µ(M_xO_y) - y·µ(H_2O)/x - G^0(M^{2y/x+}))
\end{align}
$$

Taking the exemple of brucite dissolution above, the equation above replaces reaction (1) and all other reactions/constraints are kept. In this case, the equation takes the form:


$$
\begin{align}
log(a(Mg^{2+})) - 2·log(a(H^+)) = \dfrac{1}{RT}(µ(MgO) - 2·µ(H_2O) - G^0(Mg^{2+}))
\end{align}
$$

In the case of a more complex system, for exemple a antigorite-brucite serpentinite at 400 °C and 5 kbar in the system MgO-SiO<sub>2</sub>-H<sub>2</sub>O-HCl. The additional constrain takes the form:

$$
\begin{align}
log(a(SiO_{2}(aq))) = \dfrac{1}{RT}(µ(SiO_{2}) - G^0(SiO_{2}(aq)))
\end{align}
$$

## How to modify the computation

To test this speciation calculator, you can change the matrix in the data folder and change the log Ks, line 98 (b variable). You can also change the chlorinity (Cl<sup>tot</sup> variable) line 85.

## Authors
Thibault Duretz,
Guillaume Siron
