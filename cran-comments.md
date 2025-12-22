## Test environments
* local OS X 26.2, R 4.5.2
* macOS-latest (on GitHub-Actions), R 4.5.2
* windows-latest (on GitHub-Actions), R 4.5.2
* ubuntu-latest (on GitHub-Actions), R 4.5.2
* ubuntu-latest (on GitHub-Actions), R devel
* ubuntu-latest (on GitHub-Actions), R 4.4.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs. 

There was 1 NOTE:
> N  checking installed package size
     installed size is  5.5Mb
     sub-directories of 1Mb or more:
       data      1.8Mb
       extdata   3.1Mb
  
  This is a minor release for 'driveR', fixing the CRAN errors due to strong
  dependencies on packages from Bioconductor data annotation repository. The
  packages were moved to 'Suggests' and code was updated to execute if 
  installed, raising an error if not. Once again, the large files in the 
  indicated sub-directories are necessary for the vignette and the examples 
  to work. Therefore, I could not discard them.
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
