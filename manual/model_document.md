
# Chemical Kinetics (`Chemkin`) Library User Manual

## 1. Introduction

`chemical_kinetics` is a simple library for handling chemical reaction systems. Chemical reactions are input using a standard XML format, and reaction rates can be computed.

The clients could call the chemkin package and obtained the right-hand-side of an ODE. They can then use it as the righ-hand-side of the ODE, or in a neural net code to learn new reaction pathways.
 

## 2. Installation
You can obtain the `Chemkin` Library [here](https://github.com/cs207-2017-group13/cs207-FinalProject)


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
 Take a dictionaries of properties from the XMLReader class. Calculate the rate coefficient for each elementary reaction and pass it to the ReactionSystem class. Returns a dictionary of recatants and products to the ReactionSystem class.
 
 Example:
 ```python
 
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
