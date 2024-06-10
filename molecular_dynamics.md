# Molecular dynamics
```python{cmd}
from fairchem.core.preprocessing import AtomsToGraphs
from fairchem.core.datasets import LmdbDataset
import ase.io
from ase.build import bulk
from ase.build import fcc100, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
# from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase import units
import matplotlib.pyplot as plt
import lmdb
import pickle
from tqdm import tqdm
import torch
import os

# Generate toy dataset: Relaxation of CO on Cu
adslab = fcc100("Cu", size=(2, 2, 3))
ads = molecule("CO")
add_adsorbate(adslab, ads, 3, offset=(1, 1))
cons = FixAtoms(indices=[atom.index for atom in adslab if (atom.tag == 3)])
adslab.set_constraint(cons)
adslab.center(vacuum=13.0, axis=2)
adslab.set_pbc(True)
adslab.set_calculator(EMT())

temperature_K = 300
steps = 100
traj_name = "CuCO_adslab.traj"

MaxwellBoltzmannDistribution(adslab, temperature_K=temperature_K)
# dyn = Langevin(adslab, trajectory=traj_name, timestep=1.0*units.fs, temperature_K=temperature_K, friction=0.10/units.fs)
dyn = VelocityVerlet(adslab, trajectory=traj_name, timestep=1.0*units.fs)
dyn.run(steps=steps)
```