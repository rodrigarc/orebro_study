# orebro_study
 

## Installation

- Install Conda and Bioconda. Use the
  [Bioconda instructions](https://bioconda.github.io/) if you
  don’t have Conda and/or Bioconda.
- Optionally, but highly recommended is to install [Mamba](https://github.com/mamba-org/mamba),
  which is a faster alternative to Conda:
```
      conda install mamba
```
  If you don't install mamba, write `conda` instead of `mamba` in the next step,
  but the installation will take longer.

- Create a new Conda environment and install the dependencies into it:
```
      mamba env create -n OrebroEnv -f environment.yml
```
- Activate the environment

      conda activate OrebroEnv
