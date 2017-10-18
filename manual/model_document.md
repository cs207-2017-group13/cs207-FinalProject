
# Chemical Kinetics (`Chemkin`) Library User Manual

## 1. Introduction

`chemical_kinetics` is a simple library for handling chemical reaction systems. Chemical reactions are input using a standard XML format, and reaction rates can be computed.

The clients could call the chemkin package and obtained the right-hand-side of an ODE. They can then use it as the righ-hand-side of the ODE, or in a neural net code to learn new reaction pathways.
 

## 2. Installation
You can obtain the `Chemkin` Library [here](https://github.com/cs207-2017-group13/cs207-FinalProject).
Users can run the test suite by calling pytest from the main directory, i.e. 
```python
!pytest --cov=src
```


## 3. Basic Usage and Examples
#### `XMLReader` class: Read and parse XML input file

The user should have an XML input file containing all the chemical reactions. `XMLReader` will read in the file and parse the file with the `xml.etree` library. It will output a list of dictionaries containing all the elements needed to calculate reaction rate coefficients, progress rates and reaction rates. It will create `ElementaryReaction` objects inside the `get_reaction_systems` function, and then put all `ElementaryReaction` objects in a reaction system into a list. It will then create `ReactionSystem` objects.

This class contains one method: `get_reaction_systems()`.

Example:
```python
reader = XMLReader("tests/rxns.xml")
reaction_systems = reader.get_reaction_systems()
```
`reaction_systems` is a list containing multiple `ReactionSystem` instances. The length of `reaction_systems` is the number of reaction systems, and the length of each list element is the number of reactions in a reaction system.

#### `ElementaryReaction` class: Class for each elementary reaction
Takes a dictionary of properties from the XMLReader class for each elementary reaction. Calculates the rate coefficient for each elementary reaction and passes it to the ReactionSystem class. It also returns a dictionary of recatants and products to the ReactionSystem class.

This class has thee public methods, two private methods a special method:
 - `__repr__()`: Returns a string containing basic information about the elementary reaction
 - `get_reactants()`: Returns a dictionary with reactants as key and the stoichiometric coeff of reactants as value for each elementary reaction.
 - `get_products()`:  Returns a dictionary with products as key and the stoichiometric coeff of products as value for each elementary reaction.
 - `calculate_rate_coefficient(T)`: Calculates and returns the rate coeffiecient of the reaction based on the type of rate coefficient required.
 - `_constant_rate(k=1.0)`: Returns a constant reaction rate coefficient with default value as 1.0
 - `_k_arrhenius(A, E, T, b=0.0, R=8.314)`: Returns a reaction rate coefficient according to the Arrhenius equation. The Arrhenius equation relates the rate constant, k, of a chemical reaction to parameters A (pre-exponential factor), E (activation energy), T (absolute temperature), and b (exponential indicating temperature-dependence of pre-exponential factor)::
      
    A nonzero value for b gives the modified Arrhenius equation.



$$  &k_{\textrm{mod arr}} = A T^{b} \exp\left(-\frac{E}{RT}\right) \tag{Modified Arrhenius} $$
 
    When b = 0, the above formula corresponds to the Arrhenius equation.
 
$$  &k_{\textrm{arr}}     = A \exp\left(-\frac{E}{RT}\right) \tag{Arrhenius} $$
 
 
Example:
 ```python
elementary_reaction = ElementaryReaction(reaction_properties)
reactants = elementary_reaction.get_reactants()
products = elementary_reaction.get_products()
rate_coeff = elementary_reaction.calculate_rate_coefficient(1000)
repr(elementary_reaction)
 ```
 

 



#### `ReactionSystem` class: Class for a system of reactions

Takes a list of ElementaryReaction instances and a list of species. Builds stoichiometric coefficient matrices for the reactants and products and calculates the corresponding progress rates and reaction rates.

This class has five methods, and two special methods:
 - `__repr__`: Returns a string containing basic information for the reaction system
 - `__len__`: Returns the number of species in the reaction system
 - `calculate_progress_rate(concs, temperature)`: Returns the progress rate of a system of irreversible, elementary reactions
 - `calculate_reaction_rate(concs, temperature)`: Returns the reaction rate of a system of irreversible, elementary reactions
 - `get_rate_coefficients()`: Calculate reaction rate coefficients
 - `build_reactant_coefficient_matrix()`: Build a reactant coefficients matrix for the reaction system
 - `build_product_coefficient_matrix()`: Build a product coefficients matrix for the reaction system

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
