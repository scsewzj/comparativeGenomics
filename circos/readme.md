# Creating CIRCOS plots


## Dec 03, 2025 

**Minh Ngoc** generated input to plot synteny in circos:

The homology files are in `circos/input`. This was extracted from BLASTP results, filtering was done minimally to plot the synteny:
```bash
$ ls circos/input
homology_ob.txt     # use this if you plot 1 unique circos for O.branchy
homology_os.txt     # use this if you plot 1 circos for O.sativa
homology_osob.txt   # use this if for a combined circos and if you only want to show cross-species synteny
homology_all.txt    # use this if for a combined circos and if you want to show both within- and cross-species syntenies
```
