#' Create Coding Impact Meta-prediction Score Data Frame
#'
#' @param annovar_csv_path path to ANNOVAR csv output file
#'
#' @return data frame of meta-prediction scores containing 2 columns: \describe{
#'   \item{gene_symbol}{HGNC gene symbol}
#'   \item{metaprediction_score}{metapredictor impact score}
#' }
#'
#' @examples
#' path2annovar_csv <- system.file("extdata/imielinski.hg19_multianno.csv",
#'                                 package = "driveR")
#' metapred_df <- driveR:::create_metaprediction_score_df(path2annovar_csv)
create_metaprediction_score_df <- function(annovar_csv_path) {
    anno_bechmark_df <- utils::read.csv(annovar_csv_path)

    # filter out missing scores
    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$SIFT_score != ".", ]
    anno_bechmark_df$SIFT_score <- as.numeric(anno_bechmark_df$SIFT_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$Polyphen2_HDIV_score != ".", ]
    anno_bechmark_df$Polyphen2_HDIV_score <- as.numeric(anno_bechmark_df$Polyphen2_HDIV_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$LRT_score != ".", ]
    anno_bechmark_df$LRT_score <- as.numeric(anno_bechmark_df$LRT_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$MutationTaster_score != ".", ]
    anno_bechmark_df$MutationTaster_score <- as.numeric(anno_bechmark_df$MutationTaster_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$MutationAssessor_score != ".", ]
    anno_bechmark_df$MutationAssessor_score <- as.numeric(anno_bechmark_df$MutationAssessor_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$FATHMM_score != ".", ]
    anno_bechmark_df$FATHMM_score <- as.numeric(anno_bechmark_df$FATHMM_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$GERP.._RS != ".", ]
    anno_bechmark_df$GERP.._RS <- as.numeric(anno_bechmark_df$GERP.._RS)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$phyloP7way_vertebrate != ".", ]
    anno_bechmark_df$phyloP7way_vertebrate <- as.numeric(anno_bechmark_df$phyloP7way_vertebrate)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$CADD_phred != ".", ]
    anno_bechmark_df$CADD_phred <- as.numeric(anno_bechmark_df$CADD_phred)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$VEST3_score != ".", ]
    anno_bechmark_df$VEST3_score <- as.numeric(anno_bechmark_df$VEST3_score)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$SiPhy_29way_logOdds != ".", ]
    anno_bechmark_df$SiPhy_29way_logOdds <- as.numeric(anno_bechmark_df$SiPhy_29way_logOdds)

    anno_bechmark_df <- anno_bechmark_df[anno_bechmark_df$DANN_score != ".", ]
    anno_bechmark_df$DANN_score <- as.numeric(anno_bechmark_df$DANN_score)

    # Predict metapredictor probabilities
    pred_df <- stats::predict(metapredictor_model, newdata = anno_bechmark_df, type = "prob")
    anno_bechmark_df$metaprediction_score <- pred_df$non.neutral

    metapred_scores_df <- anno_bechmark_df[, c("Gene.refGene", "metaprediction_score")]
    colnames(metapred_scores_df) <- c("gene_symbol", "metaprediction_score")

    # keep first symbol if multiple symbols exist
    metapred_scores_df$gene_symbol[grepl(";", metapred_scores_df$gene_symbol)] <- sapply(metapred_scores_df$gene_symbol[grepl(";", metapred_scores_df$gene_symbol)], function(x) unlist(strsplit(x, ";"))[1])

    # keep highest score
    metapred_scores_df <- metapred_scores_df[order(metapred_scores_df$metaprediction_score, decreasing = TRUE), ]
    metapred_scores_df <- metapred_scores_df[!duplicated(metapred_scores_df$gene_symbol), ]
    return(metapred_scores_df)
}

#' Create SCNA Score Data Frame
#'
#' @param scna_df the SCNA segments data frame. Must contain: \describe{
#'   \item{chr}{chromosome the segment is located in}
#'   \item{start}{start position of the segment}
#'   \item{end}{end position of the segment}
#'   \item{log2ratio}{\ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} ratio of the segment}
#' }
#' @param log2_ratio_threshold the \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}}
#' ratio threshold for keeping high-confidence SCNA events (default = 0.25)
#' @param gene_overlap_threshold the percentage threshold for the overlap between
#' a segment and a transcript (default = 25). This means that if only a segment
#' overlaps a transcript more than this threshold, the transcript is assigned
#' the segment's SCNA event.
#' @param MCR_overlap_threshold the percentage threshold for the overlap between
#' a gene and an MCR region (default = 25). This means that if only a gene
#' overlaps an MCR region more than this threshold, the gene is assigned the
#' MCR's SCNA density
#'
#' @return data frame containing gene-level SCNA events (that are concordant with
#' MCRs) with proxy SCNA scores (`CNVdensity..CNV.Mb.`)

#' @return data frame of SCNA proxy scores containing 2 columns: \describe{
#'   \item{gene_symbol}{HGNC gene symbol}
#'   \item{SCNA_density}{SCNA proxy score. SCNA density (SCNA/Mb) of the minimal common region (MCR) in which the gene is located.}
#' }
#'
#' @details The function first aggregates SCNA \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} ratio
#' on gene-level (by keeping the ratio with the maximal \ifelse{html}{\out{|log<sub>2</sub>|}}{\eqn{|log_2|}}
#' ratio over all the SCNA segments overlapping a gene). Next, it identifies the
#' minimal common regions (MCRs) that the genes overlap and finally assigns the
#' SCNA density (SCNA/Mb) values as proxy SCNA scores.
#'
#' @examples
#' \dontrun{
#' SCNA_scores_df <- driveR:::create_SCNA_score_df(imielinski_scna_table)
#' }
create_SCNA_score_df <- function(scna_df,
                                 log2_ratio_threshold = 0.25,
                                 gene_overlap_threshold = 25,
                                 MCR_overlap_threshold = 25){
    ### format check
    nec_cols <- c("chr", "start", "end", "log2ratio")
    if (!all(nec_cols %in% colnames(scna_df))) {
        stop("`scna_df` must contain all of: ",
             paste(dQuote(nec_cols), collapse = ", "))
    }

    #### determine gene-level log2-ratios
    # GRanges objects
    scna_gr <- GenomicRanges::makeGRangesFromDataFrame(scna_df,
                                                       keep.extra.columns = TRUE)
    GenomeInfoDb::seqlevelsStyle(scna_gr) <- "UCSC"
    genes_gr <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene))

    # overlaps
    hits <- GenomicRanges::findOverlaps(scna_gr,
                                        genes_gr,
                                        type="any", select="all")

    # obtain overlap data frame between SCNA segments and gene transcripts
    ranges <- genes_gr[S4Vectors::subjectHits(hits)]
    S4Vectors::mcols(ranges) <- c(S4Vectors::mcols(ranges),
                                  S4Vectors::mcols(scna_gr[S4Vectors::queryHits(hits)]))

    genes_df <- data.frame(chr = as.vector(GenomeInfoDb::seqnames(ranges)),
                           as.data.frame(S4Vectors::mcols(ranges)))

    `%>%` <- dplyr::`%>%`
    genes_df <- as.data.frame(genes_df %>%
                                  dplyr::mutate(segment_start = as.integer(GenomicRanges::start(GenomicRanges::ranges(scna_gr[S4Vectors::queryHits(hits)])))) %>%
                                  dplyr::mutate(segment_end = as.integer(GenomicRanges::end(GenomicRanges::ranges(scna_gr[S4Vectors::queryHits(hits)])))) %>%
                                  dplyr::mutate(segment_length_Mb = round((as.numeric((segment_end - segment_start + 1) / 1e6)), digits = 4)) %>%
                                  dplyr::mutate(transcript_start = GenomicRanges::start(ranges)) %>%
                                  dplyr::mutate(transcript_end = GenomicRanges::end(ranges)) %>%
                                  dplyr::mutate(chrom = as.character(GenomeInfoDb::seqnames(ranges))) %>%
                                  dplyr::rowwise() %>%
                                  dplyr::mutate(transcript_overlap_percent =
                                                    round(as.numeric((min(transcript_end, segment_end) - max(segment_start, transcript_start)) / (transcript_end - transcript_start)) * 100, digits = 2)))

    # filter for `gene_overlap_threshold`
    genes_df <- genes_df[genes_df$transcript_overlap_percent >= gene_overlap_threshold, ]

    # add gene symbols
    all_syms_tbl <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)
    genes_df$symbol <- all_syms_tbl$symbol[match(genes_df$gene_id, all_syms_tbl$gene_id)]
    genes_df <- genes_df[!is.na(genes_df$symbol), ]

    # discard sex chromosomes
    genes_df <- genes_df[!genes_df$chr %in% c("chrX", "chrY"), ]

    # aggregate as max |log2-ratio| over all values per gene
    agg_ratios <- tapply(genes_df$log2ratio, genes_df$symbol, function(x) {
        max_val <- max(abs(x))
        if (max_val %in% x)
            return(max_val)
        return(-max_val)
    })

    gene_agg_df <- genes_df[, c("symbol", "chrom", "transcript_start", "transcript_end")]
    gene_agg_df <- gene_agg_df[!duplicated(gene_agg_df$symbol), ]
    gene_agg_df$agg_log2_ratio <- agg_ratios[match(gene_agg_df$symbol, names(agg_ratios))]

    # filter for high confidence
    gene_agg_df <- gene_agg_df[abs(gene_agg_df$agg_log2_ratio) > log2_ratio_threshold, ]

    ### assign MCR score to overlapping genes
    MCR_gr <- GenomicRanges::makeGRangesFromDataFrame(MCR_table, keep.extra.columns = TRUE)
    agg_gr <- GenomicRanges::makeGRangesFromDataFrame(gene_agg_df,
                                                      seqnames.field = "chrom",
                                                      start.field = "transcript_start",
                                                      end.field = "transcript_end",
                                                      keep.extra.columns = TRUE)
    # overlaps
    hits <- GenomicRanges::findOverlaps(agg_gr,
                                        MCR_gr,
                                        type="any", select="all")

    # obtain overlap data frame between MCR regions and gene-level SCNA
    ranges <- MCR_gr[S4Vectors::subjectHits(hits)]
    S4Vectors::mcols(ranges) <- c(S4Vectors::mcols(ranges), S4Vectors::mcols(agg_gr[S4Vectors::queryHits(hits)]))

    final_scna_df <- data.frame(chr = as.vector(GenomeInfoDb::seqnames(ranges)),
                                as.data.frame(S4Vectors::mcols(ranges)))

    final_scna_df <- as.data.frame(final_scna_df %>%
                                       dplyr::mutate(transcript_start = as.integer(GenomicRanges::start(GenomicRanges::ranges(agg_gr[S4Vectors::queryHits(hits)])))) %>%
                                       dplyr::mutate(transcript_end = as.integer(GenomicRanges::end(GenomicRanges::ranges(agg_gr[S4Vectors::queryHits(hits)])))) %>%
                                       dplyr::mutate(MCR_start = GenomicRanges::start(ranges)) %>%
                                       dplyr::mutate(MCR_end = GenomicRanges::end(ranges)) %>%
                                       dplyr::mutate(chrom = as.character(GenomeInfoDb::seqnames(ranges))) %>%
                                       dplyr::rowwise() %>%
                                       dplyr::mutate(MCR_overlap_percent =
                                                         round(as.numeric((min(MCR_end, transcript_end) - max(transcript_start, MCR_start)) / (MCR_end - MCR_start)) * 100, digits = 2)))
    # filter for `MCR_overlap_threshold`
    final_scna_df <- final_scna_df[final_scna_df$MCR_overlap_percent >= MCR_overlap_threshold, ]

    # keep only MCR-concordant SCNA events
    final_scna_df$scna_type <- ifelse(final_scna_df$agg_log2_ratio > 0, "Amp", "Del")
    final_scna_df$MCR_concordant <- final_scna_df$MCR_type == final_scna_df$scna_type
    final_scna_df <- final_scna_df[final_scna_df$MCR_concordant, ]

    final_scna_df <- final_scna_df[, c("symbol", "SCNA_density")]
    colnames(final_scna_df)[1] <- "gene_symbol"

    return(final_scna_df)
}
