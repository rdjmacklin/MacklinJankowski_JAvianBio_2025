################################################################################

# Script Title: Set-up script for Macklin & Jankowski (2025).

################################################################################

# Install librarian package if not already installed.

if(!requireNamespace("librarian"), quietly = TRUE) {
  install.packages("librarian")
}

# Create required directories not preserved by the repo file structure.

necessary.dirs <- c("./Data/Raw/eBird-Basic-Dataset", "./Data/Raw/eBird-Basic-Dataset/CBCH",
                    "./Data/Raw/eBird-Basic-Dataset/BCCH", "./Data/Raw/eBird-Basic-Dataset/EBD-Sampling",
                    "./Data/Raw/BCR_Terrestrial", "./Data/Raw/GEDI", "./Data/Raw/eBird-Status-and-Trends",
                    "./Data/Raw/BC_Regional-Districts", "./Data/Cleaned", "./Data/Cleaned/eBird-Basic-Dataset",
                    "./Data/Cleaned/eBird-Basic-Dataset/Filtered", "./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH",
                    "./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH", "./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered",
                    "./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered/BCCH", "./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered/CBCH",
                    "./Data/Cleaned/Environmental-Variables", "./Data/Cleaned/Prediction-Grids", "./Outputs",
                    "./Outputs/Plots", "./Outputs/Tables", "./Scripts/RData")

for(i in necessary.dirs) {
  
  if(!dir.exists(i)){
    
    dir.create(i)
    
  }
  
}