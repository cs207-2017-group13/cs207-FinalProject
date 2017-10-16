# Chemical Kinetics (`Chemkin`) Library User Manual
=======================================

## 1. Introduction
------------------ 

This is a chemical kinetics library. The clients could call the chemkin package and obtained the right-hand-side of an ODE. They can then use it as the righ-hand-side of the ODE, or in a neural net code to learn new reaction pathways.
 
## 2. Installation
------------------

## 3. Basic Usage and Examples
------------------ 
- `XMLReader` class: Read and parse XML input file
The user should have an XML input file containing all the chemical reactions. `XMLReader` will read in the file and parse the file with the `xml.etree` library. It will output a list of dictionaries containing all the elements needed to calculate reaction rate coefficients, progress rates and reaction rates. It will create `ElementaryReaction` objects inside the `get_reaction_systems` function, and then put all `ElementaryReaction` objects in a reaction system into a list. It will then create `ReactionSystem` objects.
Example:
```python
reader = XMLReader("tests/rxns.xml")
reaction_systems = reader.get_reaction_systems()
```
`reaction_systems` is a list containing multiple `ReactionSystem` instances. The length of `reaction_systems` is the number of reaction systems, and the length of each list element is the number of reactions in a reaction system. 
- `ElementaryReaction` class: Class for each elementary reaction
 Take a dictionaries of properties from the XMLReader class. Calculate the rate coefficient for each elementary reaction and pass it to the ReactionSystem class. Returns a dictionary of recatants and products to the ReactionSystem class.
 Example:
 ```python

 ```

- `ReactionSystem` class: Class for a system of reactions
