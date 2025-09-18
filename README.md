# Unifying species distributions, community science and the 'natural removal experiment' to explore species interactions at broad geographic scales.

Journal of Avian Biology

Rory D.J. Macklin <sup>1,2 \*</sup> and Jill E. Jankowski <sup>1,2</sup>

<sup>\*</sup> Corresponding author <macklin@zoology.ubc.ca>

<sup>1</sup> Biodiversity Research Centre, University of British Columbia, BC, Canada

<sup>2</sup> Department of Zoology, University of British Columbia, BC, Canada

# Contents

This repository contains the code needed to reproduce the analyses conducted in Macklin & Jankowski (2025). Due to licensing restrictions, users must obtain the data necessary to complete the analyses themselves, though we have made every effort to make this as simple a process as possible. We have included links to the data sources (permanent DOIs where possible) and built the folder structure of this repository such that it should be easy to drop data files into the correct location with ease.

The four necessary scripts for running this analysis are contained in the Scripts folder and are as follows:

1. 00_Setup.R
2. 01_eBird Data Filtering and Organization.R
3. 02_Habitat Data Download and Organization.R
4. 03_Species Distribution Modeling.R

## 00_Setup

This script installs our package management tool `librarian` (<https://cran.r-project.org/web/packages/librarian/>) and creates the necessary directories for the analysis.

## 01_eBird Data Filtering and Organization.R

In 01_eBird Data Filtering and Organization.R the user will supply eBird data downloaded from <https://www.ebird.org/science/download-ebird-data-products>. To get the exact release we used in our analysis, users may have to reach out to eBird at <ebird@cornell.edu> requesting the Apr-2023 release. The user will need to download the eBird Basic Dataset Sampling Event Data and eBird Basic Dataset filtered to the Chestnut-backed Chickadee (*Poecile rufescens*; CBCH) and Black-capped Chickadee (*P. atricapillus*; BCCH) and to data from January, 2000 or later for the following 5 regions:

1. Alaska, USA
2. British Columbia, Canada,
3. Washington, USA
4. Oregon, USA
5. California, USA

Raw species data files should be stored in each species respective directory under `Data/Raw/eBird-Basic-Dataset/CBCH` or `Data/Raw/eBird-Basic-Dataset/BCCH`. The sampling event database should be stored under `Data/Raw/eBird-Basic-Dataset/EBD-Sampling`. The default names for files supplied by eBird upon download are retained in the R script, but so long as the user changes both the filename and the reference in the R script, users may customize filenames as they see fit.

## 02_Habitat Data Download and Organization.R

02_Habitat Data Download and Organization.R imports data from a variety of sources, in some cases automated via R packages. Each source will be explained below.

* The R package 'rnaturalearth' acts to gather the necessary basemap data. See <https://github.com/ropensci/rnaturalearth> and <https://cran.r-project.org/web/packages/rnaturalearth/index.html>. Code is included to facilitate data download, but in case this code becomes defunct over time, we include the data modified and compiled into a geopackage under `Data/Raw/Natural-Earth` as the necessary data are released under the public domain (<https://www.naturalearthdata.com/about/terms-of-use/>).
* Bird Conservation Region shapefiles are downloaded within the script from <https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip>. Should this link die, alternative sources of the data may be found at <https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions> or by contacting <birdmap@birdscanada.org>. As we lack licensing information, we do not include the raw data files in this repository. These data should be stored in `Data/Raw/BCR_Terrestrial`.
* MODIS MCD12Q1 v 6.1 data are downloaded via the R package 'MODIS' (see <https://github.com/fdetsch/MODIS>). At the time of publication, users must register with a NASA Earthdata Account to access this data at <https://urs.earthdata.nasa.gov/users/new>. As we could not find licensing information releasing this data for redistribution, we do not include it here.
* Users must supply the derived elevation, slope, aspect northness, and aspect eastness data from Amatulli et al. (2018). At the time of publication, this data was available for download from <https://www.earthenv.org/topography> and archived permanently at <http://doi.org/10.17616/R31NJMVX>. The below screenshot gives an example of the settings used to download this data from the website. These files retain their default naming convention and are stored in `Data/Raw/Elevation`.

  <img width="788" height="558" alt="Screen Shot 2025-09-17 at 4 58 33 PM" src="https://github.com/user-attachments/assets/f2e6ba5f-5ce0-4446-8f5c-0762ec4a886a" />

* Users must supply the mean (average) nighttime light intensity rasters from Elvidge et al. (2021). These data are downloadable from <https://eogdata.mines.edu/products/vnl/> under the header 'Annual VNL V2'. We used version 2.1. If users wish to access the data, they will need a NASA Earth Observation Group Login, obtainable by selecting "Register" after clicking the "Version 2.1" button on the website supplied above. From each annual directory, we downloaded the respective file ending in `.average.dat.tif.gz`, unzipped once downloaded, and stored in the `Data/Raw/VIIRS` directory, retaining the original names.
* Users must supply the road vector data for North America from the GRIP Global Roads Database of Meijer et al. (2018). We accessed this data from <https://www.globio.info/download-grip-dataset>, and it is archived permanently at <https://zenodo.org/records/6420961>. These files should be stored in the `Data/Raw/GRIP4-Vector` directory.
* Global canopy height data from Potapov et al. (2021) were downloaded from <https://glad.geog.umd.edu/dataset/gedi/>. We downloaded the `Forest_height_2019_NAM.tif` file. In case this link dies, the associated publication can be accessed via the permanent identifier at <https://doi.org/10.1016/j.rse.2020.112165> and the corresponding author can be contacted at <potapov@umd.edu>. As there is no licensing information provided with the data, we do not include it here. These data should be stored in `Data/Raw/GEDI`.
* Users must supply the Daymet rasters derived from Thornton et al. (2022) in the `Data/Raw/Daymet` directory. Data can be accessed at <https://doi.org/10.3334/ORNLDAAC/2130> should users wish to access the data themselves.

With these datafiles, the user will be able to recreate both the checklist-level covariate datasets and the prediction grid across the study area needed to conduct the analysis.

# 03_Species Distribution Modeling.R

This file is used to conduct the bulk of the analysis and requires a few more data files that we have not already imported. We describe the process for obtaining this data here.
* We provide the eBird Status and Trend range polygons for CBCH and BCCH as the distribution of these data for the purposes of scientific publication is permitted (<https://science.ebird.org/en/status-and-trends/products-access-terms-of-use>). We accessed these from the eBird Status and Trends Data Version 2021 release (Fink et al. 2022). Should users wish to obtain the data themselves, they must contact <ebird@cornell.edu> to request the version needed to exactly replicate our analysis. The two `.gpkg` files should be placed within the `Data/Raw/eBird-Status-and-Trends` directory.
* Users will need to provide the Regional District Boundaries - Road centreline aligned polygon from Elections BC. We downloaded this data from <https://catalogue.data.gov.bc.ca/dataset/regional-district-boundaries-road-centreline-aligned>. Due to licensing restrictions, we do not include this data here. Users can submit a request for this data at the supplied link. We requested the data by selecting "Access/Download" under the "BC Geographic Warehouse Custom Download" header on the right on the webpage. We requested the file in GeoJSON format. The contents of the zipped data file should be placed in the `Data/Raw/BC_Regional-Districts` directory.

<img width="1824" height="505" alt="Screen Shot 2025-09-18 at 12 34 50 PM" src="https://github.com/user-attachments/assets/d513ac20-4a1a-424a-b09f-1f9ace6033be" />

* We provide the Natural Earth Counties and Coastlines polygons, as they are under public domain (<https://www.naturalearthdata.com/about/terms-of-use/>). These data are available directly from Natural Earth at <https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/> and <https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-2-counties/>. Users procuring the data themselves should place the contents of zipped datafiles in `Data/Raw/Natural-Earth/ne_10m_coastline` and `Data/Raw/Natural-Earth/ne_10m_admin_2_counties` respectively.

With these files, the user should have all they need to conduct the full analysis. Note that this script will take a long time (i.e. days to weeks depending on processing power), at times running in parallel. The script is set to use all but one of the available cores for parallel processing. Users may adjust this option by modifying the `cores` object within the script.
