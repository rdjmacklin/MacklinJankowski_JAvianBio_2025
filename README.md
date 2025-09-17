# Unifying species distributions, community science and the 'natural removal experiment' to explore species interactions at broad geographic scales.

Journal of Avian Biology

Rory D.J. Macklin <sup>1,2 \*</sup> and Jill E. Jankowski <sup>1,2</sup>

<sup>\*</sup> Corresponding author <macklin@zoology.ubc.ca>

<sup>1</sup> Biodiversity Research Centre, University of British Columbia, BC, Canada

<sup>2</sup> Department of Zoology, University of British Columbia, BC, Canada

# Contents

This repository contains the code needed to reproduce the analyses conducted in Macklin & Jankowski (2025). Due to licensing restrictions, users must obtain the data necessary to complete the analyses themselves, though we have made every effort to make this as simple a process as possible. We have included links to the data sources (permanent DOIs where possible) and built the folder structure of this repository such that it should be easy to drop data files into the correct location with ease.

The three necessary scripts for running this analysis are contained in the Scripts folder and are as follows:

1. eBird Data Filtering and Organization.R
2. Habitat Data Download and Organization.R
3. Species Distribution Modeling.R


## eBird Data Filtering and Organization.R

In Script 1: eBird Data Filtering and Organization.R the user will supply eBird data downloaded from <https://www.ebird.org/science/download-ebird-data-products>. To get the exact release we used in our analysis, users may have to reach out to eBird at <ebird@cornell.edu> requesting the Apr-2023 release. The user will need to download the eBird Basic Dataset Sampling Event Data and eBird Basic Dataset filtered to the Chestnut-backed Chickadee (*Poecile rufescens*; CBCH) and Black-capped Chickadee (*P. atricapillus*; BCCH) and to data from January, 2000 or later for the following 5 regions:

1. Alaska, USA
2. British Columbia, Canada,
3. Washington, USA
4. Oregon, USA
5. California, USA

Raw species data files should be stored in each species respective directory under Data/Raw/eBird-Basic-Dataset/CBCH or Data/Raw/eBird-Basic-Dataset/BCCH. The sampling event database should be stored under Data/Raw/eBird-Basic-Dataset/EBD-Sampling. The default names for files supplied by eBird upon download are retained in the R script, but so long as the user changes both the filename and the reference in the R script, users may customize filenames as they see fit.
