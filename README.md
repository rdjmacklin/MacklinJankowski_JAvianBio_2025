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

## Habitat Data Download and Organization.R

Script 2: Habitat Data Download and Organization.R imports data from a variety of sources, in some cases automated via R packages. Each source will be explained below.

* The R package 'rnaturalearth' acts to gather the necessary basemap data. See <https://github.com/ropensci/rnaturalearth> and <https://cran.r-project.org/web/packages/rnaturalearth/index.html>.
* Bird Conservation Region shapefiles are downloaded within the script from <https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip>. Should this link die, alternative sources of the data may be found at <https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions> or by contacting <birdmap@birdscanada.org>.
* MODIS MCD12Q1 v 6.1 data are downloaded via the R package 'MODIS' (see <https://github.com/fdetsch/MODIS>). At the time of publication, users must register with a NASA Earthdata Account to access this data at <https://urs.earthdata.nasa.gov/users/new>.

* Users must supply the derived elevation, slope, aspect northness, and aspect eastness data from Amatulli et al. (2018). At the time of publication, this data was available for download from <https://www.earthenv.org/topography>. The below screenshot gives an example of the settings used to download this data from the website. These files should retain their default naming convention and be placed in Data/Raw/Elevation.

  <img width="788" height="558" alt="Screen Shot 2025-09-17 at 4 58 33 PM" src="https://github.com/user-attachments/assets/f2e6ba5f-5ce0-4446-8f5c-0762ec4a886a" />

* Users must supply the mean (average) nighttime light intensity rasters from Elvidge et al. (2021). These data are downloadable from <https://eogdata.mines.edu/products/vnl/> under the header 'Annual VNL V2'. We used version 2.1. Users will need a NASA Earth Observation Group Login, obtainable by selecting "Register" after clicking the "Version 2.1" button on the website supplied above. From each annual directory, we downloaded the respective file ending in ".average.dat.tif.gz", unzipped once downloaded, and stored in the Data/Raw/VIIRS directory, retaining the original names.
