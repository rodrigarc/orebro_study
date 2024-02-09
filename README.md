## Modulation of innate immune response to mRNA vaccination after SARS-CoV-2 infection or sequential vaccination in humans
Hellgren F, Rosdahl A, Arcoverde Cerveira R, Lenart K, Ols S, Gwon Y-D, Joas G, Kurt S, Delis A-M, Evander M, Normark J, Ahlm C, Forsell M, Cajander S, Loré K

Corresponding author: Karin Loré, Division of Immunology and Allergy, Department of Medicine Solna, Karolinska Institutet, Visionsgatan 4, BioClinicum J7:30, Karolinska University Hospital, 171 64 Stockholm, Sweden.

Github repo for the bioinformatic analysis performed in the paper Hellgren et.al. 2024. (link to paper DOI), 
The metadata record is available at figshare scilifelab: [figshare](https://doi.org/10.17044/scilifelab.24941913)

## Table of Contents
* [Abstract](#abstract)
* [Preprocessing dataset](#preprocessing-dataset)
* [Repository structure](#repository-structure)
* [Reproducibility](#reproducibility)
* [Plots and rendered results](#plots-and-rendered-results)
* [Contributors](contributers)
* [License](#license)

## Abstract

Abstract mRNA vaccines are likely to become widely used for the prevention of infectious diseases in the future. Nevertheless, a notable gap exists in mechanistic data, particularly concerning the potential effects of sequential mRNA immunization or pre-existing immunity on the early innate immune response triggered by vaccination. In this study healthy adults, with or without documented prior SARS-CoV-2 infection, were vaccinated with the BNT162b2/Comirnaty mRNA vaccine. Prior infection conferred significantly stronger induction of proinflammatory and type I interferon related gene signatures, serum cytokines, and monocyte expansion after the prime vaccination. The response to the second vaccination further increased the magnitude of the early innate response in both study groups. The third vaccination did not additionally augment vaccine-induced inflammation. In vitro stimulation of PBMCs with toll-like receptor ligands showed no difference in cytokine responses between groups or pre/post prime vaccination, indicating absence of a trained immunity effect. We observed that levels of pre-existing antigen-specific CD4 T cells, antibody and memory B cells correlated with elements of the early innate response to the first vaccination. Our data thereby indicate that pre-existing memory formed by infection can augment the innate immune activation induced by mRNA vaccines.

## Preprocessing dataset

## Repository structure
 - `src` folder: contains all the source code used for the analysis.
 - `orebro_study.Rproj` file: contains the R project to be open within RStudio for reproducibility purposes.
 - `results` folder contains generated plots
 - `renv` folder: contains the files generated when renv is initiated. 
 - `renv.lock` filer: contains the package versions and dependencies used to generate the plots.

## Reproducibility

### System Requirements

[R](https://www.r-project.org/) version 4.3.2.

[RStudio](https://posit.co/download/rstudio-desktop/) (2023.12.1+402)

[renv](https://rstudio.github.io/renv/index.html) (1.0.3)

[Pandoc](https://pandoc.org/) is required for knitting R Markdown documents. `rmarkdown` relies on Pandoc but does not include it within the package. As such, Pandoc is a system dependency outside the scope of `renv`. To ensure consistent document rendering, please ensure that you have Pandoc installed on your system. The version used for this project is Pandoc 3.1.1, which can be downloaded from the [Pandoc releases page](https://github.com/jgm/pandoc/releases/tag/3.1.1).

After installing Pandoc, you can verify the installation and version by running the following command in your R console:

`rmarkdown::pandoc_version()`

### Installation

1. Open the terminal and clone the repository
```
git clone git@github.com:rodrigarc/orebro_study.git
```

3. Navigate to the project directory. If you are using RStudio, open the project by clicking on the `orebro_study.Rproj` file.

4. To restore the project environment, run the following command in the R console:
```
renv::restore()
```

To just access the HTML file of the rendered R markdown run
```
bash src/render_Rmd.sh
```

### Plots and rendered results

All the plots will be generated under a `results` folder under a date folder (eg. `results/2023-02-10`) automatically created. The Rmarkdown file is rendered to html and contains the information regarding the type of analysis and its code used to generate all plots. A rendered version is already online and uploaded under `src/orebro_rnaseq-analysis.html` if you wish to just check the code used for each plot. The rendered `html` containing the code and analysis can be accessed [HERE](need to setup github pages).

##  Contributors
The following people contributed to this repository:

[Rodrigo Arcoverde](https://github.com/rodrigarc)<br />
[Gustav Joas](https://github.com/GustavDavid)<br />
[Fredrika Hellgren](https://github.com/fredrihel)

## License

All the code in this repo is under a GNU General Public License v3.0.
