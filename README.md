# nnp_testground

## How to run these examples

To use the newest `openmm-ml` version it is necessary to install from github. The following should set up a sufficient conda environment and installs both `nnpops` and `mace`.

```
conda create --name test_env python=3.11
conda activate test_env
conda install openmm 
conda install openmm-ml
conda remove --force openmm-ml
pip install git+https://github.com/openmm/openmm-ml.git  
pip install mace-torch
mamba install nnpops -c conda-forge
```


## Examples


- examples/waterbox contains simulation setup for a ~1K atom pure water simulation with periodic boundary conditions. Simulations are performed either with ANI2x or MACE.

- examples/protein_ligand contains a solvated protein-ligand simulation with periodic boundary conditions in which the ligand is treated with a NNP and the protein&water is treated with a classical force field.

