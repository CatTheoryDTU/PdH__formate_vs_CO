Package contains all the data and routines for the collaboration with the group
of Drew Higgins on the selectivity to CO and formate on PdH.

The data can be used to reproduce figure 6 and the figure S14. The raw
structure are given in the `data` directory and scripts for analysis are given
in the `scripts` directory.
The actual data is contained in `results/parsed_most_stable_data.pckl` and
`results/parsed_all_data.pckl`

In order to make figure 6a do

```
python3 scripts/plot_FED.py
```

for figure 6b do

```
cd catmap
# Make sure catmap is loaded
bash README.md
```

The README above contains additional instructions.


WARNING: Sadly, in order to parse the trajectory files properly a slightly changed ASE version from trunk is needed, as vanilla ASE does not allow readinig and writing of the number of electrons and potentials in the trajectory format. The necessary changes to ASE are:

In ase/calculators/singlepoint.py change line 25 from
```
if property in ['energy', 'magmom', 'free_energy']:
```
to
```
if property in ['energy', 'magmom', 'free_energy','ne','electrode_potential']:
```


and in ase/calculators/calculator line 98 exchange:
```
all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy']
```
with
```
all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy','ne','electrode_potential']
```

We also note that GPAW version 21.1.1b1 was used for the calculations and more recent version changed the names of `ne` and `electrode_potential`. However, ASE can still read the results if above's changes are made.

