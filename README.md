# Modulation of innate immune response to mRNA vaccination after SARS-CoV-2 infection or sequential vaccination in humans

Github repo for the bioinformatic analysis performed in the paper.

## Prerequisites
Ensure you have installed the latest version of [R](https://www.r-project.org/).

## Additional Dependencies

### Pandoc

[Pandoc](https://pandoc.org/) is required for knitting R Markdown documents. `rmarkdown` relies on Pandoc but does not include it within the package. As such, Pandoc is a system dependency outside the scope of `renv`. To ensure consistent document rendering, please ensure that you have Pandoc installed on your system. The version used for this project is Pandoc 3.1.1, which can be downloaded from the [Pandoc releases page](https://github.com/jgm/pandoc/releases/tag/3.1.1).

After installing Pandoc, you can verify the installation and version by running the following command in your R console:

`rmarkdown::pandoc_version()`


## Installation

1. Open the terminal and clone the repository
`git clone git@github.com:rodrigarc/orebro_study.git`

2. Navigate to the project directory. If you are using RStudio, open the project by clicking on the orebro_study.Rproj file.

3. To restore the project environment, run the following command in the R console:
`renv::restore()`

To just access the HTML file of the rendered R markdown run
`bash src/render_Rmd.sh`

##  Contributors
Thanks to the following people who have contributed to this project:

[Rodrigo Arcoverde](https://github.com/rodrigarc)
[Gustav Joas](https://github.com/GustavDavid)

## Contact
If you want to contact me you can reach me at rodrigo.arcoverde@ki.se.
