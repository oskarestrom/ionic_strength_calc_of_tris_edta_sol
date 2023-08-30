# ionic_strength_calc_of_tris_edta_sol
Calculating the Ionic Strength of a multivalent buffer solution such as Tris-EDTA is not straight forward. All the ions of the buffer affect each other and form an equilibrium. With Tris-EDTA, you end up with about a system of 7-9 equations to solve. 

This package calculates the Ionic Strength and pH of a Tris-EDTA solution (with or without betamercaptoethanol or BME) for you and extracts the calculated data to a .csv file.

## Quick Start Guide

### Dependencies

First install [Sympy](https://www.sympy.org/en/index.html) (version 1.10.1) by writing this in your anaconda prompt:
```shell
pip install sympy==1.10.1
```
You also need numpy and pandas:
```shell
pip install numpy >= 1.21.5
pip install pandas >= 1.4.4
```


### Installing the package
Currently the package can be found on [test.pypi.org](https://test.pypi.org/project/ionic-strength-calc-of-tris-edta-sol/0.0.2/). Install it on your anaconda prompt by:
```shell 
pip install -i https://test.pypi.org/simple/ ionic-strength-calc-of-tris-edta-sol==0.0.2
```
### Using the package
This is how a script using the package could look like:
```python
import ionic_strength_calc_of_tris_edta_sol.calc as ion
import numpy as np

#Set the concentration of the TE ions in mM in the array  concTE (here 1 mM and 5 mM):
concTE = np.array([1, 5])

#Dictionary with settings for the calculation:

settings = {
    'with_BME': False, #True if BME is present
    'c_BME_percent': 0, #BME concentration in percent
    'save_to_file': True, #True if the results should be saved to a file
    'showAllSol': False, #True if all solutions should be shown
    'showEq': False, #True if the equilibrium equations should be shown
    'showConc': False, #True if the concentrations should be shown
    'showCoeffs': False, #True if the coefficients should be shown
    'basepath' : r'', #Path to the folder where the results should be saved
}

ion.calc_ionic_strength_TE(concTE, settings)
```
### Output
See the files "df_TE_0_BME.csv" and "df_all_TE_0_BME.csv" in this repository for the output by the code above. I also post the content here:
#### Output - "df_TE_0_BME.csv":
This file contains a list of the samples used. Here, two concentrations of Tris-EDTA were used (1 mM and 5 mM). The resulting ionic strength and pH are presented.
|FIELD1|conc_TE|conc_BME|I                 |pH               |
|------|-------|--------|------------------|-----------------|
|0     |1      |0.0     |6.102551731739553 |8.422372930122124|
|1     |5      |0.0     |30.992848261160407|8.417746234353766|

#### Output - "df_all_TE_0_BME.csv":
This file contains **all** information about all ionic species in the equilibrium. This can be useful if you want to know a certain concentration of e.g. a reactive species.
|FIELD1|species_f|species|val               |pK_a             |K_a                  |pK_b|K_b                   |g                  |c                    |N_TE|I                 |pH               |c_mM                 |cont_mM              |c_BME_percent|
|------|---------|-------|------------------|-----------------|---------------------|----|----------------------|-------------------|---------------------|----|------------------|-----------------|---------------------|---------------------|-------------|
|0     |$[TH^{+}]$|[TH^{+}]|1                 |0.0              |                     |5.94|1.1481536214968817e-06|0.8480252333759359 |0.01524850479767449  |5   |30.992848261160407|8.417746234353766|15.24850479767449    |7.624252398837245    |0            |
|1     |$[T]$    |[T]    |0                 |0.0              |                     |    |                      |1.0                |0.03475149520232551  |5   |30.992848261160407|8.417746234353766|34.75149520232551    |0.0                  |0            |
|2     |$[HE^{3-}]$|[HE^{3-}]|-3                |6.16             |6.918309709189363e-07|    |                      |0.22681876200946888|0.004735671049012061 |5   |30.992848261160407|8.417746234353766|4.735671049012061    |21.310519720554275   |0            |
|3     |$[H_2 E^{2-}]$|[H_2 E^{2-}]|-2                |2.67             |0.0021379620895022326|    |                      |0.5171721148332381 |9.729442230604571e-06|5   |30.992848261160407|8.417746234353766|0.009729442230604572 |0.019458884461209143 |0            |
|4     |$[E^{4-}]$|[E^{4-}]|-4                |10.26            |5.495408738576248e-11|    |                      |0.0715385941424201 |0.000254599508757334 |5   |30.992848261160407|8.417746234353766|0.254599508757334    |2.036796070058672    |0            |
|5     |$[H^{+}]$|[H^{+}]|1                 |14.0             |1e-14                |    |                      |0.8480252333759359 |3.821675125839302e-09|5   |30.992848261160407|8.417746234353766|3.821675125839302e-06|1.910837562919651e-06|0            |
|6     |$[OH^{-}]$|[OH^{-}]|-1                |0.0              |                     |    |                      |0.8480252333759359 |3.638552822886325e-06|5   |30.992848261160407|8.417746234353766|0.0036385528228863252|0.0018192764114431626|0            |
|7     |$[TH^{+}]$|[TH^{+}]|1                 |0.0              |                     |5.94|1.1481536214968817e-06|0.8480252333759359 |0.01524850479767449  |5   |30.992848261160407|8.417746234353766|15.24850479767449    |7.624252398837245    |0            |
|8     |$[T]$    |[T]    |0                 |0.0              |                     |    |                      |1.0                |0.03475149520232551  |5   |30.992848261160407|8.417746234353766|34.75149520232551    |0.0                  |0            |
|9     |$[HE^{3-}]$|[HE^{3-}]|-3                |6.16             |6.918309709189363e-07|    |                      |0.22681876200946888|0.004735671049012061 |5   |30.992848261160407|8.417746234353766|4.735671049012061    |21.310519720554275   |0            |
|10    |$[H_2 E^{2-}]$|[H_2 E^{2-}]|-2                |2.67             |0.0021379620895022326|    |                      |0.5171721148332381 |9.729442230604571e-06|5   |30.992848261160407|8.417746234353766|0.009729442230604572 |0.019458884461209143 |0            |
|11    |$[E^{4-}]$|[E^{4-}]|-4                |10.26            |5.495408738576248e-11|    |                      |0.0715385941424201 |0.000254599508757334 |5   |30.992848261160407|8.417746234353766|0.254599508757334    |2.036796070058672    |0            |
|12    |$[H^{+}]$|[H^{+}]|1                 |14.0             |1e-14                |    |                      |0.8480252333759359 |3.821675125839302e-09|5   |30.992848261160407|8.417746234353766|3.821675125839302e-06|1.910837562919651e-06|0            |
|13    |$[OH^{-}]$|[OH^{-}]|-1                |0.0              |                     |    |                      |0.8480252333759359 |3.638552822886325e-06|5   |30.992848261160407|8.417746234353766|0.0036385528228863252|0.0018192764114431626|0            |


## Contributions
You can contribute to this project! Please send me an email if you want to collaborate on this. A major contribution would be making a general ionic-strength-calculator for all kinds of solutions. 

## How the script works
The script calculates the ionic strength based on the detailed description by [Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4) ( See the bottom of the page for a reference list). For an even more detailed description of the chemistry, look at the two papers from Persat et al.. We basically get a series of equations that needs to be solved. In these equations, the disassociation constants together with the total species concentrations are known but the concentrations of the ionic species are not. The script works iteratively, calculating the ionic strength and pH stepwise (see all steps below). 

### SymPy
This script relies on the equation solving functionality of the python library SymPy (v. 1.10.1).

See the following links for some information of equation solving in SymPy:
- https://apmonitor.com/che263/index.php/Main/PythonSolveEquations
- https://docs.sympy.org/latest/tutorial/solvers.html
- https://docs.sympy.org/latest/modules/solvers/solveset.html


Below I give brief chemistry explanations of some of the terms used. See [Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4) for a more detailed description.


### Ionic strength
The ionic strength, $I$, is calculated by summing the product of the concentration and squared charge of all ions in the solution:

$$I = \dfrac{1}{2} \sum_{i=1}^{n} c_i z_i^2$$

where the one half is added to include both anions and cations, $c_i$ is the molar concentration of the ion i (in M) and $z_i$ is the charge of the ion i. 

The buffer consists of (1x TE with BME):
- 10 mM Tris
- 1 mM EDTA
- 3% (v/v) Betamercaptoethanol (BME)

However, even though we are using TE, I first calculated the ionic strength for TBE buffer as we have good literature values to validate our calculation results. The percent error (100 × [calculated value – literature value] / literature value)
is below 0.6 % based on the values from Hsieh et al.

Note that [Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4) writes that because of most of the EDTA is triply ionized at their pH 8.5 (our pH is pH 8) they ignore the negliable concentrations of neutral and singly ionized EDTA.

We make the following assumptions:
- Valency of Tris: 1, EDTA (multivalent): 1-4, BME: 1.
- The temperature, T = 25$`^{\circ}`$ C for determining the pKs
- That the following values are correct:
    - Dissociation constant for water $pK_w$ = 14.0 at T = 25$`^{\circ}`$ C.
    - Tris base dissociation constant, $pK_{b_{Tris}}$ = 5.94 for T = 25$`^{\circ}`$ C from [Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4).
    - EDTA acid dissociation constants, $pK_{a_{EDTA}}$ = 1.99,2.67,6.16,10.26 for T = 25$`^{\circ}`$ C from [Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4).
    - BME acid dissociation constant $pK_{a_{BME}}$ = 9.6 for T = 25$`^{\circ}`$ C from [Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4).
    - The constant for the Davies equation, A = 0.51 from [Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4).
    - Molar mass for BME, $M_{BME}$ = 78.13 g/mol, Source: Sigma-Aldrich
    - Density for BME, $\rho_{BME}$ = 1.114 g/mL for T = 25$`^{\circ}`$, Source: Thomas Scientific


### Calculation of the ionic strength for TE with BME

We have the following product from Sigma-Aldrich:  - Tris-EDTA buffer solution, 100 x, for molecular biology
(Code: T9285-100ML)
The molar concentration of a species X is given by [X] and its activity coefficient given by $\gamma_X$

At N x TE, the solution contains N x 0.01 M Tris-HCl and N x 0.001 M EDTA. The pH should be approximately 8. The molar concentration of BME depends on the volumetric concentration, see calculations below.

We have the following system of equations for TE. Equations 1 to 6 are the equilibrium equations:

$$\dfrac{[TH^+]\gamma_{TH^+}[OH^-]\gamma_{OH^-}}{[T]\gamma_{T}} = 10^{-5.94} \ \  (1)$$ 

$$\dfrac{[HE^{3-}]\gamma_{HE^{3-}}[H^+]\gamma_{H^+}}{[H_2 E^{2-}]\gamma_{H_2 E^{2-}}} = 10^{-6.16} \  \ (2)$$ 

$$\dfrac{[E^{4-}]\gamma_{E^{4-}}[H^+]\gamma_{H^+}}{[HE^{3-}]\gamma_{HE^{3-}}} = 10^{-10.26} \  \ (3)$$ 

$$\dfrac{[\beta^{-}]\gamma_{\beta^{-}}[H^+]\gamma_{H^+}}{[H\beta]\gamma_{H\beta}} = 10^{-9.6} \  \ (4)$$ 

$$[H^+]\gamma_{H^+}[OH^-]\gamma_{OH^-} = 10^{-14.0} \  \ (5)$$ 

Equations 6 to 8 is to ensure the conservation of the concentration of the substances:

$$[T] + [TH^+] = C[T]  \  \ (6) $$

$$[H_2E^{2-}] + [HE^{3-}] + [E^{4-}] = C[E]  \  \ (7) $$

$$[H\beta]  + [\beta^{-}] = C[\beta]  \  \ (8) $$

Equation 9 is to ensure charge neutrality:

$$[H^+] + [TH^+] = [OH^-] + 2[H_2E^{2-}] + 3[HE^{3-}] + 4[E^{4-}] + [\beta^{-}] \  \ (9) $$

where $E$ stands for EDTA, $T$ for Tris and $\beta$ for BME.

### Henderson-Hasselbalch equation

The Henderson-Hasselbalch equation is given by
$$pH = pK_a + log \left( \dfrac{[A^-]}{[HA]} \right)$$

where $pK_a$ is the acid disassociation constant, $[A^-]$ is the concentration of the buffer's base form and $[HA]$ is the concentration of the buffer's acid form.

### The Davies equation
The Davies equation [Davies, 1938] gives the mean molal efficiency coefficient $\gamma_X$ of a n electrolyte at 25°C that dissociates into ions having valence (charge) $z_X$ as a function of ionic strength, $I$:

$$log_{10}(\gamma_X) = -0.51z_X^2 \left( \dfrac{\sqrt{I}}{1+\sqrt{I}} -0.3I \right)$$

Note: it could be worth calculating the ionic strength using Debye-Hückel model instead of the Davies equation, see [https://www.aqion.de/site/101](https://www.aqion.de/site/101)

Then, the ionic strength is given by the following for TE:

$$I = 1/2\left( [TH^+] + [H^+] + [OH^-] + 4[H_2E^{2-}] + 9[HE^{3-}] + 16[E^{4-}] + [\beta^{-}] \right)$$

### Steps to calculate the ionic strength of a multivalent ion solution:
[Iarko et al.](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.92.062701?casa_token=XRW2tXi736wAAAAA%3AKD0YkiBiHr__Hf6wHgsKtXdIQTb6tmdhWEhxoqcUC6J4nm0WNqYeHUvNyV1-pWcVZvrY2hMzQmA4)describes in the appendix C of their article how to calculate the ionic strength of a solution:
1. Guess the activity coefficients, $\gamma_X$=1 for all species. 
1. Solve the equations 1-9 and get a value of the ionic strength.
1. Get new values of the activity coefficients $\gamma_X$ based on the new value of the ionic strength.
1. Repeat the two previous steps till the ionic strength converges.

### References
Iarko, V., Werner, E., Nyberg, L. K., Müller, V., Fritzsche, J., Ambjörnsson, T., ... & Mehlig, B. (2015). Extension of nanoconfined DNA: Quantitative comparison between experiment and theory. Physical review E, 92(6), 062701.

Persat, A., Chambers, R. D., & Santiago, J. G. (2009). Basic principles of electrolyte chemistry for microfluidic electrokinetics. Part I: acid–base equilibria and pH buffers. Lab on a Chip, 9(17), 2437-2453.

Persat, A., Suss, M. E., & Santiago, J. G. (2009). Basic principles of electrolyte chemistry for microfluidic electrokinetics. Part II: Coupling between ion mobility, electrolysis, and acid–base equilibria. Lab on a Chip, 9(17), 2454-2469.

Hsieh, C. C., Balducci, A., & Doyle, P. S. (2008). Ionic effects on the equilibrium dynamics of DNA confined in nanoslits. Nano letters, 8(6), 1683-1688.
