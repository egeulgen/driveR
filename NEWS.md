# driveR 1.0.1.9001
## to be released as 1.0.1

## Major Changes

## Minor changes and bug fixes
- Fixed issue in `create_metaprediction_score_df()` where gene symbols were turned into vectors
- To prevent warnings, in `create_SCNA_score_df()`, if there is no overlap between gene-level SCNAs and MCRs, an empty data frame is returned
