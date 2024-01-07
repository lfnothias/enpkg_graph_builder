# enpkg_graph_builder
Build the Experimental Natural Products Knowledge Graph

⚙️ Workflow part of [enpkg_workflow](https://github.com/enpkg/enpkg_workflow).

The aim of this repository is to format as RDF the data produced with GNPS/SIRIUS/tima.

## 0. Clone repository and install environment

1. Clone this repository.
2. Create environment: 
```console 
conda env create -f environment.yml
```
3. Activate environment:  
```console 
conda activate graph_builder
```

## Running it

```console
python src/rdf_builder.py --folder input --ion neg --ion_sirius auto --ion --ion
```

### Parameters
```
'--folder', required=True, help='The path to the directory where samples folders to process are located'
```
```
'-ion', '--ionization_mode', required=True, choices=['pos', 'neg'], help='The ionization mode to perform spectral library matching'
```
```
'-c', '--cpus', Number of cpu to use. Default is 80% of available CPUs.'
```
```
'-r', '--recompute', Recompute even if the files are already present
```
