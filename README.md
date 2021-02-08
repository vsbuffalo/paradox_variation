# Code and Data for *Why do species get a thin slice of π?*

This repository contains the code and data for Buffalo (2021), *Why do species
get a thin slice of π?  Revisiting Lewontin’s Paradox of Variation*.

## Overview

 - `data`: contains the data from other papers I have used, as well as all of
   the code to collate these data sources, estimate species ranges and
   population sizes, and create the final clean combined datasets. See the
   `README.md` file there for details; the `Makefile` in this repository
   generates all required datasets (note, for efficiency, there is a lot of
   caching that occurs).

- `notebooks/`: contains the main analysis Rmarkdown file
  (`main_analysis.Rmd`). This document contains the primary analyses used in
  the paper.

- `notebooks/figures/`: contains the R code for each figure. The `Makefile` in
  this repository creates these figures.


