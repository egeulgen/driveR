## Test environments
* local OS X 12.0.1, R 4.1.2
* Ubuntu 16.04.7 LTS (on Travis-CI), R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs. 

There was 1 NOTE:
> checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      data      1.7Mb
      extdata   3.1Mb

  This is a minor update to driveR. The large files in the indicated 
  sub-directories are necessary for the vignette and the examples to work. 
  Therefore, I could not discard them.
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
