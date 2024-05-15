# Energy calculation
* In this page, we see the basic usage of ASE. This starts with calculating the potential energy, which is a fudamental property of atoms/molecules/bulks/surfaces.
* Here, we assume the readers already have rough knowledge on Python. We use the following issues in Python, so if you are not sure about them please refer some basic textbook/webpage for Python.
    + importing library
    + list and tuple
    + function and its argument (positional argument and keyword argument)
    + class (use only)
    + for loop

# Calculating energy of molecules
* In the following script, we calculate the potential energy of slab with EMT potential.
* The main object in ASE is `Atoms` object. Atoms object is the list of `Atom` object, and Atom object has attributes like an element symbol, positions (in Cartesian coordinate).
* So first thing to do with ASE is define the Atoms object, and assign the elements and positions of consisting Atom object.

```python{cmd}
from ase import Atoms
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.visualize import view

mol = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.7]])
mol.calc = EMT()
e_mol = mol.get_potential_energy()

print(f"Potential energy is {e_mol:.3f} eV")

view(mol)
```
* `Atoms("H2", positions=...)` speficies the H2 (hydrogen) molecules as Atoms. Positions of the two H atoms are given by the list or tuple with three numbers, corresponding to x-, y-, and z-coordinates. So in this case, the positions are "list of lists".
* This script uses following functions:
    + `ase.calculators.emt()`: calculator to use EMT potential
    + `ase.get_potential_energy()`: calculate potential energy using the defined calculator
    + `ase.visualize.view()`: show the molecular structure with viewer

## Exercise
1. Try to calculate the potential energy with different bond length.
2. Try to calculate the potential energy of CO molecule.

# Calculating the energy of bulk material
* In this section, we calculate the potential energy of bulk mateial, which is periodic in three-dimensions.
* The following script treats the copper (Cu) bulk material.
* The simulation is done for **unit cell**, which is an elemental block of the bulk material. The periodic condition is applied to the unit cell, so it is repeated in x-, y-, and z-directions to represent the bulk.
```python{cmd}
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.visualize import view

bulk = bulk("Cu", crystalstructure="fcc", a=3.6)
bulk.calc = EMT()
e_bulk = bulk.get_potential_energy()

print(f"Potential energy is {e_bulk:.3f} eV")

view(bulk)
```
* `ase.build.bulk` makes the bulk structure of specified element. `crystalstructure` should be like `sc`, `fcc`, `bcc` etc. which means simple cubic, face-centered-cubic, and body-centered-cubic structures. `a` is the lattice constant in Angstrom (10^-10 meter) unit, which is the size of the simulation box.
* You can use square simulation box by setting `orthorhombic=True` in `bulk`.
* Lager unit cell can be used repeating the bulk: this can be made with `bulk = bulk*[2, 2, 2]` or `bulk *= [2, 2, 2]`. This extended unit cell is called **supercell**.

## Exercise
1. Compare the bulk structures of `orthorhombic=False` (default) and `True`.
2. Compare the potential energy of bulk material with different lattice constant (a) value.
3. Try to use `for` loop to do the above calculation, and make the plot using `matplotlib`.

# Calculating the energy of surfaces
* In this section, we calculate the potential energy of a surface or slab, in a very similar way to the above.
* Slab is the model for the surfaces, which is a little different from the bulk.
* Surface is periodic in two-diemsion, so we don't have any repeated strucuture in some one dimension. Usugally this direction is set to z-direction.
* Surfaces play important role in material science, because it can act as the place for the chemcial reaction; external molecules cannot enter into the bulk material but they can attach to the surface of bulk material. This often induces chemical reactions.
* Surface is made by "cutting" the bulk material. Thus, we need to define the direction of cutting, like horizontally cutting or vertically cutting the butter. This direction is specified with **Miller index**, which has three integer numbers (e.g. "111", "100").
* Since a bulk material has its crystal structure (like fcc, bcc), the speficiation of the surface should be like "fcc-111", "fcc-100", "bcc-100".

```python{cmd}
from ase import Atoms
from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.visualize import view

slab = fcc111("Cu", size=(3, 3, 3), a=3.6, vacuum=10.0)
slab.calc = EMT()
e_slab = slab.get_potential_energy()

print(f"Potential energy is {e_slab:.3f} eV")

view(slab)
```
* `ase.build.fcc111`: gives `Atoms` object for slab.
    + `size`: size of the slab, speficied by tuple
    + `a`: lattice constant in Angstrom
    + `vacuum`: length of the vacuum region in z-direction

# Calculating adsorption energy
* Here, we calculate the **adsorption energy** of a molecule on a surface.
* An adsorption energy is the energy required or released upon the adsorption.
* According to the definition of the potential energy,
    + positive adsorption energy: energy is required to adsorb
    + negative adsorption energy: energy is released upon adsorption
* Therefore more negative adsorption energy means stronger binding of molecule on the surface.
* Adsorbing atoms or molecules are called **adsorbates**.
* Adsorption energy can be calculation from three terms:
    + $E_{mol}$: potential energy of adsorbate
    + $E_{slab}$: potential energy of slab
    + $E_{mol+slab}$: potential energy of adsorbate plus slab
    + $E_{ads} = E_{mol+slab} - (E_{mol} + E_{slab})$

```python{cmd}
from ase import Atoms
from ase.build import add_adsorbate, fcc111
from ase.calculators.emt import EMT
from ase.visualize import view

height = 1.9

slab = fcc111("Cu", size=(4, 4, 2), vacuum=10.0)

slab.calc = EMT()
e_slab = slab.get_potential_energy()

molecule = Atoms("N2", positions=[[0.0, 0.0, 0.0], [0.0, 0.0 ,1.2]])
molecule.calc = EMT()
e_N2 = molecule.get_potential_energy()

add_adsorbate(slab=slab, adsorbate=molecule, height=height, position="ontop")

e_adslab = slab.get_potential_energy()

print(f"Adsorption energy is {e_adslab - (e_N2 + e_slab):.3} eV")

view(slab)
```
* `ase.build.add_adsorbate` has following arguments:
    + `slab`: Atoms for slab.
    + `adsorbate`: Atoms for adsorbate.
    + `height`: Height above the surface.
    + `position`: The x- and y-position of the adsorbate. Two numbers (list or tuple) or keyword like "ontop", "fcc", "hcp" can be used.
    + `offset`: Offsets the adsorbate position by a number of unit cell (not supercell). `offset=(1,1)` means adsorbate position is shifted to the next unit cell. Fractional numbers can be used for this argument, so `offset=(0.5, 0.5)` places the adsorbate in the middle of the cell.

## Exercise
1. Try using `fcc100` instead of fcc111. Just use `ase.build.fcc100`.
2. Compare the potential energies of different height values.