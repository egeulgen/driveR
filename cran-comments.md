## Test environments
* local OS X 12.5, R 4.2.1
* macOS-latest (on GitHub-Actions), R 4.2.1
* windows-latest (on GitHub-Actions), R 4.2.1
* ubuntu-latest (on GitHub-Actions), R 4.2.1
* ubuntu-latest (on GitHub-Actions), R devel
* ubuntu-latest (on GitHub-Actions), R 4.1.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs. 

There was 1 NOTE:
> checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      data      1.7Mb
      extdata   3.1Mb

  This is a minor update to driveR, adding support for GRCh38. The large files 
  in the indicated sub-directories are necessary for the vignette and the 
  examples to work. Therefore, I could not discard them.
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
