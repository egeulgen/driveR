# driveR 0.2.0

## Major Changes

- fixed an issue in `create_SCNA_score_df()` where the SCNA score was not calculated because the column name for "MCR_overlap_percent" was incorrectly assigned as "transcript_overlap_percent"
- updated the MTL classification model `TCGA_MTL_fit` and `specific_thresholds` after fixing the issue above

## Minor changes and bug fixes

- Updated citation information, correcting author name
- Minor changes in utility functions
- added the `na.string` argument to `create_noncoding_impact_score_df()`, `predict_coding_impact()` and `create_features_df()` as the string that was used to indicate when a score is not available during annotation with ANNOVAR (default = ".")
- generalized `determine_hotspot_genes()` to able to use occurrence annotations from different versions of COSMIC 
- updated all data

***

# driver 0.1.1

Initial release
