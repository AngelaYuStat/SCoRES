## Test environments
* Local R installation, Windows 11, R 4.4.0
* GitHub Actions: macOS-arm64 (R-devel), Linux (R-devel), Windows (R-devel)
* All checks completed successfully on these platforms.

## R CMD check results
There were **no ERRORs**, **WARNINGs** or **no NOTEs**.

There was one **INFO** message:

checking installed package size ... INFO
installed size is 5.0Mb
sub-directories of 1Mb or more:
data 4.0Mb

This size is expected and justified.  
The data directory contains essential **climate and geospatial datasets** used in the examples and vignettes.  
These datasets are necessary to demonstrate the package’s core functionality, including spatial modeling and reproducible analyses.  
All data have been appropriately reduced, subsetted, and compressed to the minimal size required for functionality testing and illustration.

## Downstream dependencies
There are currently no downstream dependencies.

## Additional notes
* The package passes R CMD check --as-cran on all major operating systems.
* This is the first submission of the package **SCoRES (version 0.1.0)** to CRAN.
