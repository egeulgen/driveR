# driveR 0.4.0

## Major Changes
- Added GRCh38 support

***

# driveR 0.3.0

## Major Changes
- Updated the cancer-type-specific thresholds

## Minor changes and bug fixes
- Updated citation to the method's article where necessary

***

# driveR 0.2.1

## Major Changes

- Updated the MTL model after fixing the issue below
- Fixed issue in `MCR_table`. The coordinates were converted to hg19 (from hg18)

## Minor changes and bug fixes

- Made use of `caret::predict.train` explicit
- Removed 'Homo.sapiens' from Imports field (not used)

***

# driveR 0.2.0

## Major Changes

- Fixed an issue in `create_SCNA_score_df()` where the SCNA score was not calculated because the column name for "MCR_overlap_percent" was incorrectly assigned as "transcript_overlap_percent"
- Updated the MTL classification model `TCGA_MTL_fit` and `specific_thresholds` after fixing the issue above

## Minor changes and bug fixes

- Updated citation information, correcting author name
- Minor changes in utility functions
- Added the `na.string` argument to `create_noncoding_impact_score_df()`, `predict_coding_impact()` and `create_features_df()` as the string that was used to indicate when a score is not available during annotation with ANNOVAR (default = ".")
- Generalized `determine_hotspot_genes()` to able to use occurrence annotations from different versions of COSMIC 
- Updated all data

***

# driver 0.1.1

Initial release
