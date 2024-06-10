```python{cmd}
cluster = Icosahedron("Ru", noshells=4, latticeconstant=3.5)
```

```python{cmd}
from ase.visualize import view
from ase.cluster import Decahedron
cluster = Decahedron("Ru", p=4, q=4, r=0, latticeconstant=3.5)
view(cluster)
```

```python{cmd}
cluster = FaceCenteredCubic("Ru", surfaces=[(1,0,0),(1,1,1)], layers=[3,3], latticeconstant=3.5)
```

```python{cmd}
cluster = Octahedron("Ru", length=6, cutoff=2, latticeconstant=3.5)
```

## Wulff construction
To set up a Wulff construction, the surface energies should be specified, in units of energy per area (not energy per atom). The actual unit used does not matter, as only the ratio between surface energies is important. In addition, the approximate size of the nanoparticle should be given. As the Wulff construction is build from whole layers, it is not possible to hit the desired particles size exactly:


```python{cmd}
cluster = wulff_construction("Ru", surfaces=[(1,1,1),(1,0,0)], energies=[1.0, 1.3], structure="fcc", size=150, latticeconstant=3.5)
```