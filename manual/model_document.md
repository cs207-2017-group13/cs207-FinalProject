# Chemical Kinetics (`chemical_kinetics`) Library User Manual

## 1. Introduction

`Chemical kinetics` is the study of chemical reactions with respect to reaction rates, effect of factors like  temperature, activation energy, stoichiometric coeffecients of reactants and products, etc. on these reactions.

The library `chemical_kinetics` is a simple library for handling such chemical reaction systems. A system of chemical reactions are input in a standard XML format and the library computes the reaction rates and progress rate of these reactions, based on the parameters passed in the XML.

In a system of $N$ species undergoing $M$ **elementary** reactions, the reaction rate of species $i$ is computed by - 

$$
\begin{align}
  f_{i} &= \sum_{j=1}^{M}{\nu_{ij}\omega_{j}}, \qquad i = 1, \ldots, N
\end{align}
$$
The progress rate for each reaction is computed by 
$$
\begin{align}
  \omega_{j} &= k_{j}^{\left(f\right)}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime}}} - k_{j}^{\left(b\right)}\prod_{i=1}^{N}{x_{i}^{\nu_{ij}^{\prime\prime}}}, \qquad j = 1,\ldots, M
\end{align}
$$
where,
$$
\begin{align}
  k_{j}^{\left(b\right)} &= \frac{k_{j}^{\left(f\right)}}{k_{j}^{e}}, \qquad j =1, \ldots, M\\
  k_{j}^{e} &= \left(\frac{p_{0}}{RT}\right)^{\gamma_{j}}\exp\left(\frac{\Delta S_{j}}{R} - \frac{\Delta H_{j}}{RT}\right), \qquad j =1, \ldots, M\\
  \gamma_{j} &= \sum_{i=1}^{N}{\nu_{ij}}\\
  \Delta S_{j} &= \sum_{i=1}^{N}{\nu_{ij}S_{i}} \quad \textrm{and} \quad \Delta H_{j} = \sum_{i=1}^{N}{\nu_{ij}H_{i}}, , \qquad j =1, \ldots, M\\
  H_{i} &= \int_{T_{0}}^{T}{C_{p,i}\left(T\right) \ \mathrm{d}T}, \qquad i = 1, \ldots, N \\
  S_{i} &= \int_{T_{0}}^{T}{\frac{C_{p,i}\left(T\right)}{T} \ \mathrm{d}T}, \qquad i = 1, \ldots, N\\
  C_{p,i} &= \left(\sum_{k=1}^{5}{a_{ik}T^{k-1}}\right)R, \qquad i = 1, \ldots, N
\end{align}
$$

The $7$th order NASA polynomials are given by 
$$\frac{C_{p,i}}{R} = a_{i1} + a_{i2}T + a_{i3}T^{2} + a_{i4}T^{3} + a_{i5}T^{4}$$
$$\frac{H_{i}}{RT} = a_{i1} + \frac{1}{2}a_{i2}T + \frac{1}{3}a_{i3}T^{2} + \frac{1}{4}a_{i4}T^{3} + \frac{1}{5}a_{i5}T^{4} + \frac{a_{i6}}{T}$$
$$\frac{S_{i}}{R} = a_{i1}\ln\left(T\right) + a_{i2}T + \frac{1}{2}a_{i3}T^{2} + \frac{1}{3}a_{i4}T^{3} + \frac{1}{4}a_{i5}T^{4} + a_{i7}$$
for $i = 1,\dots, N$.

Notation:

$\nu_{ij}^{\prime}$ : Stoichiometric coefficients of reactants,

$\nu_{ij}^{\prime\prime}$ : Stoichiometric coefficients of products,

$\omega_{j}$ : Progress rate of reaction $j$,

$x_{i}$ : Concentration of specie $i$,

$k_{j}^{\left(f\right)}$: forward reaction rate coefficient for reaction $j$,

$k_{j}^{\left(b\right)}$: backward reaction rate coefficient for reaction $j$,

$k_{j}^{e}$: *equilibrium coefficient* for reaction $j$,

$p_{0}$: pressure of the reactor (usually $10^{5}$ Pa),

$\Delta S_{j}$: the entropy change of reaction $j$,

$\Delta H_{j}$: the enthalpy change of reaction $j$,

$C_{p,i}$: specific heat at constant pressure (given by the NASA polynomial)

The clients could call the chemkin package and obtained the right-hand-side of an ODE. They can then use it as the righ-hand-side of the ODE, or in a neural net code to learn new reaction pathways.
 

## 2. Installation
You can obtain the `chemical_kinetics` Library [here](https://github.com/cs207-2017-group13/cs207-FinalProject).

### Installation instructions
Obtain the latest version from github, change to the directory, and install using pip.

    git clone https://github.com/cs207-2017-group13/cs207-FinalProject.git
	cd cs207-FinalProject

	# Either
	pip3 install ./
	# or alternatively
	pip3 install -e ./

The latter install command invokes "editable" mode. Package files will not be copied to your Python package directory, and source files may be edited from the installation directory.

### Testing
Users can run the test suite by calling pytest from the main directory, i.e. 

    pytest --cov=src


## 3. Basic Usage and Examples
### `XMLReader` class: Read and parse XML input file

The user should have an XML input file containing all the chemical reactions. `XMLReader` will read in the file and parse the file with the `xml.etree` library. It will output a list of dictionaries containing all the elements needed to calculate reaction rate coefficients, progress rates and reaction rates. It will create `ElementaryReaction` objects inside the `get_reaction_systems` function, and then put all `ElementaryReaction` objects in a reaction system into a list. It will then create `ReactionSystem` objects.

This class contains one method: `get_reaction_systems()`.

Example:
```python
reader = XMLReader("tests/rxns.xml")
reaction_systems = reader.get_reaction_systems()
```
`reaction_systems` is a list containing multiple `ReactionSystem` instances. The length of `reaction_systems` is the number of reaction systems, and the length of each list element is the number of reactions in a reaction system.


### `ElementaryReaction` class: Class for each elementary reaction

Takes a dictionary of properties from the XMLReader class for each elementary reaction. Calculates the rate coefficient for each elementary reaction and passes it to the ReactionSystem class. It also returns a dictionary of recatants and products to the ReactionSystem class.

This class has thee public methods, two private methods a special method:
 - `__repr__()`: Returns a string containing basic information about the elementary reaction
 - `get_reactants()`: Returns a dictionary with reactants as key and the stoichiometric coeff of reactants as value for each elementary reaction.
 - `get_products()`:  Returns a dictionary with products as key and the stoichiometric coeff of products as value for each elementary reaction.
 - `calculate_rate_coefficient(T)`: Calculates and returns the rate coeffiecient of the reaction based on the type of rate coefficient required.
 - `_constant_rate(k=1.0)`: Returns a constant reaction rate coefficient with default value as 1.0
 - `_k_arrhenius(A, E, T, b=0.0, R=8.314)`: Returns a reaction rate coefficient according to the Arrhenius equation. The Arrhenius equation relates the rate constant, k, of a chemical reaction to parameters A (pre-exponential factor), E (activation energy), T (absolute temperature), and b (exponential indicating temperature-dependence of pre-exponential factor)::
      
A nonzero value for b gives the modified Arrhenius equation.

$$
\begin{align}
  &k_{\textrm{mod arr}} = A T^{b} \exp\left(-\frac{E}{RT}\right) \tag{Modified Arrhenius}
\end{align}
$$ 

When b = 0, the above formula corresponds to the Arrhenius equation.

$$ 
\begin{align}
  &k_{\textrm{arr}}     = A \exp\left(-\frac{E}{RT}\right) \tag{Arrhenius}
\end{align}
$$

Example:

```python
elementary_reaction = ElementaryReaction(reaction_properties)
reactants = elementary_reaction.get_reactants()
products = elementary_reaction.get_products()
rate_coeff = elementary_reaction.calculate_rate_coefficient(1000)
repr(elementary_reaction)
 ```


### `ReactionSystem` class: Class for a system of reactions

Takes a list of ElementaryReaction instances and a list of species. Builds stoichiometric coefficient matrices for the reactants and products and calculates the corresponding progress rates and reaction rates.

This class has five methods, and two special methods:
 - `__repr__`: Returns a string containing basic information for the reaction system
 - `__len__`: Returns the number of species in the reaction system
 - `calculate_progress_rate(concs, temperature)`: Returns the progress rate of a system of elementary reactions
 - `calculate_reaction_rate(concs, temperature)`: Returns the reaction rate of a system of elementary reactions
 - `get_rate_coefficients()`: Calculate reaction rate coefficients
 - `build_reactant_coefficient_matrix()`: Build a reactant coefficients matrix for the reaction system
 - `build_product_coefficient_matrix()`: Build a product coefficients matrix for the reaction system
 - `check_reversible`: Check if each elementary reaction is reversible

Example:
```python
len(reaction_system[0])
concs = [1., 2., 1., 3., 1.]
reaction_system[0].build_reactant_coefficient_matrix()
reaction_system[0].build_product_coefficient_matrix()
print(reaction_system[0])
reaction_system[0].calculate_progress_rate(concs, 300)
reaction_system[0].calculate_reaction_rate(concs, 300)
```


### `Thermochem` class: Class for calculating the backward reaction rate

Construct with Rxnset class object, the default pressure of the reactor and the default ideal gas constant. The values of pressure of the reactor and ideal gas constant can be changed. It calculates backward reaction rate using the temperature passed into the functions and the corresponding NASA polynomial coefficients.

This class has four methods:
 - `Cp_over_R(T)`: Returns specific heat of each species given by the NASA polynomials
 - `H_over_RT(T)`: Returns the enthalpy of each species given by the NASA polynomials
 - `S_over_R(T)`: Returns the entropy of each species given by the NASA polynomials
 - `backward_coeffs(kf, T)`: Returns the backward reaction rate coefficient for reach reaction

Example:
```python
reader = chemkin.XMLReader("tests/rxns_reversible.xml")
reaction_system = reader.get_reaction_systems()[0]
nu = reaction_system.product_coefficients - reaction_system.reactant_coefficients
rxnset = thermodynamics.Rxnset(reaction_system.species, nu)
thermo =  thermodynamics.Thermochem(rxnset)
thermo.Cp_over_R(800)
thermo.H_over_RT(800)
thermo.S_over_R(800)
kf = reaction_system.get_rate_coefficients(800)
thermo.backward_coeffs(kf, 800)
```


### `Rxnset` class: Read and store NASA polynomial coefficients 

This class reads the NASA polynomial coefficients for all species in the reaction system from the SQL database which contains coefficients for all species. It stores the coefficients and temperature ranges in a dictionary of dictionaries where the name of species is the key. For each reaction system, the class just need to read from the database once, and check the range the given temperature is in every time the temperature of the reaction system changes afterwards. If the `get_nasa_coefficients` function is called twice for the same reaction system and the same temperature, the cached value is returned.

This class has two methods:
 - `get_nasa_coefficients(T)`: Returns the corresponding NASA polynomial coefficients for all species at the given temperature
 - `read_nasa_coefficients()`: Return NASA polynomial coefficients for all species involved in the reaction system

Example:
```python
reader = chemkin.XMLReader("tests/rxns_reversible.xml")
reaction_system = reader.get_reaction_systems()[0]
nu = reaction_system.product_coefficients - reaction_system.reactant_coefficients
rxnset = thermodynamics.Rxnset(reaction_system.species, nu)
rxnset.get_nasa_coefficients(800)
```


### `memoized` class: Caches a function's return value each time it is called

Decorator class. Caches a function's return value each time it is called. If called later with the same arguments, the cached value is returned (not reevaluated).

This class has three special methods:
 - `__call__(*args)`: Cached a function's return value
 - `__repr__`: Return the function's docstring
 - `__get__(obj, objtype)`: Support instance methods


