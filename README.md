##INTRODUCTION
This project analyzes data published by Ling et al., "Extremely high genetic diversity in a single tumor points to prevalence of non-Darwinian cell evolution" (2015). In this paper, the authors extracted DNA from 286 tumor samples collected from a single hepatocellular carcinoma (hcc) cross section. WES was carried out on 23 of these samples and 35 unique SNPs were genotyped via PCR in the remaining samples.
##FILE DESCRIPTIONS
alleles.csv -- raw data
locations.csv -- pixel locations of each sample measured from Figure 2A using Preview
clones.csv -- clone labels of each sample
CancerProjectIData.ipynb -- contains script to clean and parse raw data
alleles_fast_meta.txt -- output of CancerProjectIData.ipynb for BEAST analysis
CancerProjectI-BALTIC.ipynb -- analysis of trees draw from posterior distribution in BEAST output 
