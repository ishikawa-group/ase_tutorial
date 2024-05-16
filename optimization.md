# Geometry optimization
* In this section, we will see how to do the geometry optimization with ASE.
* Geometry optimization tries to locate the molecular/bulk/surface structure that has the minimum potential energy. By changing the position of consisting atoms, the potentital energy rise or fall, so geometry optimization reduces the potential energy step by step and finally leads to the energy minimum.
* There are several algorithms to perform the geometry optimization:
    + BFGS

* The following script does the geometry optimization of a molecule.

```python{cmd}
from ase import Atoms
from ase.build import molecule
from ase.optimize.bfgs import BFGS
from ase.calculators.emt import EMT
import os

mol = molecule("H2O")
mol.set_calculator(EMT())
opt = BFGS(mol, trajectory="test.traj")
opt.run(fmax=0.01, steps=100)

os.system("ase gui test.traj")
os.system("rm test.traj")
```

* The above script uses
    * `ase.optimize.bfgs`: Uses BFGS optimizer, and the result is seved on the `trajectory` file. The optimization starts with `run` method, and it has `fmax`(force threshold for convergence) and `steps`(the number of maximum steps).
    * The trajectory file can be viewed with `ase gui`. To call it, `os.system` is used (but now it is recommended to use `subprocess` instead). `ase.visualize.view` gives the last (optimized) structure but does not show the whole trajectory.

## Exercise
1. Try to do optimize the Cu surface.

# Optimization with constraint
* When performing the slab (surface) calculation, it is common to fix several part of the slab (usually lower part).
* This is because we would like mimic the bulk behavior in the calculation: in reality, the bulk material has a large number of atomic layers (in z-direction, for example) and only some atomic layers near the surface would change positions while remaining part does not move.
* In calculation, however, we cannot use such large number of atomic layers. Thus, to mimic the above situmation, we constrain (or fix) the atomic position of lower part of the slab.
* The geometry optimization with constraint is done in the following way.

```python{cmd}
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.calculators.emt import EMT
from ase.optimize import BFGS
import subprocess

surf = fcc111("Cu", size=[4,4,4], vacuum=10.0)
indices = [atom.index for atom in surf if atom.tag > 2]
constraint = FixAtoms(indices)
surf.set_constraint(constraint)
surf.set_calculator(EMT())

opt = BFGS(surf, trajectory="test.traj")
opt.run(fmax=0.01, steps=100)

subprocess.run("ase gui test.traj", shell=True)
```

* The constraint is set by following steps:
    1. set indices that has the atom index to fix
    2. make `FixAtoms` object that has the index
    3. set constraint for Atoms by using `set_constraint` method (of Atoms)
* When building the surface, the index is set automatically. The top layer has index=1, and the second layer has index=2, etc. So, by getting atom's index list whose tag is larger than 2 (3 and 4 in this case), we can make the index list to fix. The atom's index is stored in `atom.index` attribute.
* We used *inline comprehension* to set index. This is very popular Python technique.

## Exercise
1. Try to optimize H2-adsorbed Cu surface, by fixing the lower half part of the slab.
