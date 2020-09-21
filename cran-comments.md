## Test environments
* local OS X 10.15.6, R 4.0.2
* Ubuntu 16.04.6 LTS (on Travis-CI), R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs. 

There were 2 NOTEs:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ege Ulgen <egeulgen@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Genomics (2:47)
  KEGG (20:75)
  driveR (15:58, 18:5)
  hotspot (19:44)
  metaprediction (18:31)

* checking for future file timestamps ... NOTE
unable to verify current time

  This is the first submission for driveR. The second NOTE occurs because the 
  resource R CMD check uses (http://worldclockapi.com/) is currently not 
  available. 
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
