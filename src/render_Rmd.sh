#!/bin/bash

# Rendering Rmarkdown script to HTML file

Rscript -e "rmarkdown::render('orebro_rnaseq-analysis.Rmd')"

# For more info on rendering options please see:
# https://pkgs.rstudio.com/rmarkdown/reference/render.html
