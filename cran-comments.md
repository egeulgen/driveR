## Test environments
* local OS X 13.5.1, R 4.3.1
* macOS-latest (on GitHub-Actions), R 4.3.1
* windows-latest (on GitHub-Actions), R 4.3.1
* ubuntu-latest (on GitHub-Actions), R 4.3.1
* ubuntu-latest (on GitHub-Actions), R devel
* ubuntu-latest (on GitHub-Actions), R 4.2.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs. 

There was 1 NOTE:
> N  checking installed package size
     installed size is  5.4Mb
     sub-directories of 1Mb or more:
       data      1.7Mb
       extdata   3.1Mb

  This is a patch release to driveR, fixing the bug in package documentation as
  pointed out by CRAN. Once again, the large files in the indicated 
  sub-directories are necessary for the vignette and the examples to work. 
  Therefore, I could not discard them.
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
