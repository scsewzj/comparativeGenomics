# Gene duplication analysis for Oryza sativa japonica

**Minh Ngoc VU**

*These are snippets of a more detailed log file, generated only to communicate instant results with my teammates. The full log file will be uploaded when all analyses are completed.*

## Dec 02, 2025 - MCScanX output

BLASTP results were filtered using 2 sets of thresholds to create 2 datasets: H (highly stringent, `osativa_H`) and L (low stringent `osativa_L`). These were then used as input for MCScanX to classify duplication type.

The classification results for the H and L datasets are in:
```bash
$ ls MCScanX_output 
osativa_H.gene_type
osativa_L.gene_type
```

Each file (tab-delimited) has 2 columns: transcript id and duplication type:
```
Os01t0100100-01	0
Os01t0100200-01	0
Os01t0100300-00	0
Os01t0100466-00	0
Os01t0100400-01	1
```
The duplication types are encoded from 0 to 4:

```
Type_of_dup	        Code
Singleton	        0
Dispersed	        1
Proximal	        2
Tandem	            3
WGD_or_segmental	4
```

Further details can be found in `log/15.mcscan.log`.
