Molecular Dynamics Simulation, Modeling and Docking/Virtual Screening Projects/Scripts
---
**NB:** It is important to pay attention to which python version is required for each process

- Modeller scripts require python2
- PDB coordinate files retrieval script requires python3: Can be activated by running `. config.sh`

Using PyMOL command line
---
Launch command line only

```
pymol -c
```

Launch command line and get help with a command

```
pymol -c -d 'help load'
```

Execute commands

```
pymol -c -d 'load 6VYB.pdb; indicate c. A; save 6VYB_chainA.pdb'
```
This loads a pdb coordinate file containing mltiple chains, selects the 'A' chain, and saves it


