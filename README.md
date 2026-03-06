# CBS_Sox2_manuscript
This repository contains the scripts and pipelines used in the manuscript "Effects of CTCF on the regulatory landscape of the mouse Sox2 locus" by Eder et al. 

## Introduction
In this study we investigated the location dependent effect of ectopic CTCF binding sites (CBS) with and without a Sox2 promotore reporter on the expression of the mouse Sox2 locus. We created high resolution functional maps of the locus by relocating either CBS-only or CBS-Sox2P reporter constructs to thousdands of alternative positions. 

## Repository guide
This repository contains information about
- Tagmentation: Tn5_tagmentation mapping (pipeline)
  - including all config files used for the data in the paper
    
- Info_files
  - file containing FACS sorting statistics for each experiment (needed for expression score calculation)
  - file containing FI from sorted experiments and population (needed for expression score calculation)

- General functions: includes the calculation of expression score
- Figures: R markdown files to create paper figures
- Pipline for genotyping 1xCBS_Fw_Sox2P clones
