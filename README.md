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
* [Contributors](contributors)
* [License](#license)

## Abstract

Abstract mRNA vaccines are likely to become widely used for the prevention of infectious diseases in the future. Nevertheless, a notable gap exists in mechanistic data, particularly concerning the potential effects of sequential mRNA immunization or pre-existing immunity on the early innate immune response triggered by vaccination. In this study healthy adults, with or without documented prior SARS-CoV-2 infection, were vaccinated with the BNT162b2/Comirnaty mRNA vaccine. Prior infection conferred significantly stronger induction of proinflammatory and type I interferon related gene signatures, serum cytokines, and monocyte expansion after the prime vaccination. The response to the second vaccination further increased the magnitude of the early innate response in both study groups. The third vaccination did not additionally augment vaccine-induced inflammation. In vitro stimulation of PBMCs with toll-like receptor ligands showed no difference in cytokine responses between groups or pre/post prime vaccination, indicating absence of a trained immunity effect. We observed that levels of pre-existing antigen-specific CD4 T cells, antibody and memory B cells correlated with elements of the early innate response to the first vaccination. Our data thereby indicate that pre-existing memory formed by infection can augment the innate immune activation induced by mRNA vaccines.

## Preprocessing dataset

For bulk RNASeq analysis, whole blood was collected into PAXgene® Blood RNA tubes (PreAnalytiX, BD Biosciences) according to manufacturer instructions and stored at -20°C until sequencing. For library preparation, tubes were thawed overnight at RT and RNA extracted using the PAXGene Blood RNA kit (PreAnalytiX, Qiagen) according to manufacturer instructions. Total RNA concentration was measured using a NanoDrop ND-1000 spectrophotometer (Thermo Fisher Scientific, MA, USA). RNA integrity was assessed using the Agilent RNA ScreenTape assay and Agilent 2200 TapeStation (Agilent Technologies, Santa Clara, CA, USA) according to manufacturer instructions. In preparation for Illumina sequencing, isolation of mRNA, cDNA synthesis, anchor ligation, amplification and library indexing were performed using the Illumina ® Stranded mRNA Prep kit according to manufacturer instructions. Library yield was quantified by Qubit fluorometer (Thermo Fisher Scientific, MA, USA) and quality assessed by Agilent TapeStation. Indexed DNA libraries were normalized, pooled and sequenced using a NovaSeq 6000 instrument, S4 flowcell, 2x150 base pairs, in paired-end mode.


## Repository structure
 - `src`: all the source code used for the analysis.
 - `orebro_study.Rproj`: the R project to open for reproducibility purposes.
 - `results`: generated plots.
 - `renv`: files generated when renv is initiated.
 - `renv.lock`: contains the package versions and dependencies used to generate the plots.

## Reproducibility

### System Requirements

[R](https://www.r-project.org/) version 4.3.2.

[RStudio](https://posit.co/download/rstudio-desktop/) (2023.12.1+402)

[renv](https://rstudio.github.io/renv/index.html) (1.0.3)

[Pandoc](https://pandoc.org/) is required for knitting R Markdown documents. `rmarkdown` relies on Pandoc but does not include it within the package. As such, Pandoc is a system dependency outside the scope of `renv`. To ensure consistent document rendering, please ensure that you have Pandoc installed on your system. The version used for this project is Pandoc 3.1.1, which can be downloaded from the [Pandoc releases page](https://github.com/jgm/pandoc/releases/tag/3.1.1).

After installing Pandoc, you can verify the installation and version by running the following command in your R console: `rmarkdown::pandoc_version()`

### Download and run analysis

Open the terminal and clone the repository:
```
git clone git@github.com:rodrigarc/orebro_study.git
```

Navigate to the project directory. Open the project in RStudio by clicking on the `orebro_study.Rproj` file.

To restore the project environment, run the following command in the R console:
```
renv::restore()
```

### Plots and rendered results

The plots generated are automatically stored in the `results` directory, organized by date (e.g., results/2023-02-10). The Rmarkdown file is rendered to HTML, detailing the types of analyses performed and the corresponding code for generating all plots. A pre-rendered version is available online which can be accessed [HERE]() or at `src/orebro_rnaseq-analysis.html`.

##  Contributors
The following people contributed to this repository:

[Rodrigo Arcoverde](https://github.com/rodrigarc)<br />
[Gustav Joas](https://github.com/GustavDavid)<br />
[Fredrika Hellgren](https://github.com/fredrihel)

## License

All the code in this repo is under a GNU General Public License v3.0.
