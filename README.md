# Phylogenetic analysis of spatial and clonal dynamics in hepatocellular carcinoma

#### Maya Lewinsohn<sup>1,2</sup>, Trevor Bedford<sup>1</sup>

<sup>1</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA, <sup>2</sup>MSTP Program, University of Washington, Seattle, WA, USA

## Introduction

This project analyzes data published by Ling et al., "Extremely high genetic diversity in a single tumor points to prevalence of non-Darwinian cell evolution" (2015). In this paper, the authors extracted DNA from 286 tumor samples collected from a single hepatocellular carcinoma (hcc) cross section. WES was carried out on 23 of these samples and 35 unique SNPs were genotyped via PCR in the remaining samples.

## File descriptions

### Data

* [alleles.csv](data/alleles.csv) -- raw data
* [locations.csv](data/locations.csv) -- pixel locations of each sample measured from Figure 2A using Preview
* clones.csv -- clone labels of each sample *TB: I don't see this file*
* [alleles_fast_meta.txt](data/alleles_fasta_meta.txt) -- output of CancerProjectIData.ipynb for BEAST analysis
* [CancerProjectIData.ipynb](CancerProjectIData.ipynb) -- contains script to clean and parse raw data *TB: move this script to data/ directory*

### Analysis

* [CancerProjectI-BALTIC.ipynb](CancerProjectI-BALTIC.ipynb) -- analysis of trees draw from posterior distribution in BEAST output
