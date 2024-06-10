# Manipulating Atoms

## getter
* `get_atomic_numbers()`
* `get_initial_charges()`
* `get_charges()`
* `get_chemical_symbols()`
* `get_initial_magnetic_moments()`
* `get_magnetic_moments()`
* `get_masses()`
* `get_momenta()`
* `get_forces()`
* `get_positions()`
* `get_potential_energies()`
* `get_scaled_positions()`
* `get_stresses()`
* `get_tags()`
* `get_velocities()`

## setter
* `set_atomic_numbers()`
* `set_initial_charges()`
* `set_chemical_symbols()`
* `set_initial_magnetic_moments()`
* `set_masses()`
* `set_momenta()`
* `set_positions()`
* `set_scaled_positions()`
* `set_tags()`
* `set_velocities()`

## attribute
* `symbol`: element name of Atom

## Example
* Replacing elements
```python{cmd}
# Here we will replace C atom to Si
from ase import Atoms

mol  = Atoms("C2H6")
elem = ["Si" if atom.symbol == "C" else atom.symbol for atom in mol]
mol.set_chemical_symbols(elem)
```