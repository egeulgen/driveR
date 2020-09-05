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

#' Create Gene-level SCNA Data Frame
#'
#' @param scna_df the SCNA segments data frame. Must contain: \describe{
#'   \item{chr}{chromosome the segment is located in}
#'   \item{start}{start position of the segment}
#'   \item{end}{end position of the segment}
#'   \item{log2ratio}{\ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} ratio of the segment}
#' }
#' @param gene_overlap_threshold the percentage threshold for the overlap between
#' a segment and a transcript (default = 25). This means that if only a segment
#' overlaps a transcript more than this threshold, the transcript is assigned
#' the segment's SCNA event.
#'
#' @return data frame of gene-level SCNA events, i.e. table of genes overlapped
#' by SCNA segments.
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' gene_scna_df <- driveR:::create_gene_level_scna_df(imielinski_scna_table)
#' }
create_gene_level_scna_df <- function(scna_df, gene_overlap_threshold = 25) {
    ### argument checks
    nec_cols <- c("chr", "start", "end", "log2ratio")
    if (!all(nec_cols %in% colnames(scna_df))) {
        stop("`scna_df` should contain all of: ",
             paste(dQuote(nec_cols), collapse = ", "))
    }

    if (!is.numeric(gene_overlap_threshold))
        stop("`gene_overlap_threshold` should be numberic")

    if (gene_overlap_threshold < 0 | gene_overlap_threshold > 100)
        stop("`gene_overlap_threshold` should be between 0-100")

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
                                  dplyr::mutate(segment_length_Mb = round((as.numeric((.data$segment_end - .data$segment_start + 1) / 1e6)), digits = 4)) %>%
                                  dplyr::mutate(transcript_start = GenomicRanges::start(ranges)) %>%
                                  dplyr::mutate(transcript_end = GenomicRanges::end(ranges)) %>%
                                  dplyr::mutate(chrom = as.character(GenomeInfoDb::seqnames(ranges))) %>%
                                  dplyr::rowwise() %>%
                                  dplyr::mutate(transcript_overlap_percent =
                                                    round(as.numeric((min(.data$transcript_end, .data$segment_end) - max(.data$segment_start, .data$transcript_start)) / (.data$transcript_end - .data$transcript_start)) * 100, digits = 2)))

    # filter for `gene_overlap_threshold`
    genes_df <- genes_df[genes_df$transcript_overlap_percent >= gene_overlap_threshold, ]

    # add gene symbols
    all_syms_tbl <- as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)
    genes_df$symbol <- all_syms_tbl$symbol[match(genes_df$gene_id, all_syms_tbl$gene_id)]
    genes_df <- genes_df[!is.na(genes_df$symbol), ]

    return(genes_df)
}

#' Create SCNA Score Data Frame
#'
#' @inheritParams create_gene_level_scna_df
#' @param log2_ratio_threshold the \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}}
#' ratio threshold for keeping high-confidence SCNA events (default = 0.25)
#' @param MCR_overlap_threshold the percentage threshold for the overlap between
#' a gene and an MCR region (default = 25). This means that if only a gene
#' overlaps an MCR region more than this threshold, the gene is assigned the
#' MCR's SCNA density
#'
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
    ### argument checks
    if (!is.numeric(log2_ratio_threshold))
        stop("`log2_ratio_threshold` should be numberic")

    if (!is.numeric(MCR_overlap_threshold))
        stop("`MCR_overlap_threshold` should be numberic")

    if (MCR_overlap_threshold < 0 | MCR_overlap_threshold > 100)
        stop("`MCR_overlap_threshold` should be between 0-100")

    #### determine gene-level log2-ratios
    # gene-level SCNA df
    genes_df <- create_gene_level_scna_df(scna_df = scna_df,
                                          gene_overlap_threshold = gene_overlap_threshold)

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

    `%>%` <- dplyr::`%>%`
    final_scna_df <- as.data.frame(final_scna_df %>%
                                       dplyr::mutate(transcript_start = as.integer(GenomicRanges::start(GenomicRanges::ranges(agg_gr[S4Vectors::queryHits(hits)])))) %>%
                                       dplyr::mutate(transcript_end = as.integer(GenomicRanges::end(GenomicRanges::ranges(agg_gr[S4Vectors::queryHits(hits)])))) %>%
                                       dplyr::mutate(MCR_start = GenomicRanges::start(ranges)) %>%
                                       dplyr::mutate(MCR_end = GenomicRanges::end(ranges)) %>%
                                       dplyr::mutate(chrom = as.character(GenomeInfoDb::seqnames(ranges))) %>%
                                       dplyr::rowwise() %>%
                                       dplyr::mutate(MCR_overlap_percent =
                                                         round(as.numeric((min(.data$MCR_end, .data$transcript_end) - max(.data$transcript_start, .data$MCR_start)) / (.data$MCR_end - .data$MCR_start)) * 100, digits = 2)))
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

#' Determine Hotspot Containing Genes
#'
#' @inheritParams create_metaprediction_score_df
#' @param threshold (integer) threshold for the minimum number of cases with
#' a certain mutation in COSMIC (default = 5)
#'
#' @return vector of gene symbols of genes containing hotspot mutation(s)
#'
#' @examples
#' path2annovar_csv <- system.file("extdata/imielinski.hg19_multianno.csv",
#'                                 package = "driveR")
#' hotspot_genes <- driveR:::determine_hotspot_genes(path2annovar_csv)
determine_hotspot_genes <- function(annovar_csv_path, threshold = 5L) {
    # argument check
    if (!is.numeric(threshold))
        stop("`threshold` should be numeric")

    annovar_df <- utils::read.csv(annovar_csv_path)

    # parse COSMIC occurences
    annovar_df$coding_occurence <- vapply(annovar_df$cosmic91_coding,
                                          function(x) {
                                              tmp <- unlist(strsplit(x, ";"))[2]
                                              tmp <- sub("OCCURENCE=", "", tmp)
                                              tmp <- unlist(strsplit(tmp, ","))
                                              tmp <- as.numeric(sub("\\(.*\\)", "", tmp))
                                              return(sum(tmp))
                                          }, 1)

    annovar_df$noncoding_occurence <- vapply(annovar_df$cosmic91_noncoding,
                                             function(x) {
                                                 tmp <- unlist(strsplit(x, ";"))[2]
                                                 tmp <- sub("OCCURENCE=", "", tmp)
                                                 tmp <- unlist(strsplit(tmp, ","))
                                                 tmp <- as.numeric(sub("\\(.*\\)", "", tmp))
                                                 return(sum(tmp))
                                             }, 1)

    annovar_df$coding_occurence <- ifelse(is.na(annovar_df$coding_occurence), 0, annovar_df$coding_occurence)
    annovar_df$noncoding_occurence <- ifelse(is.na(annovar_df$noncoding_occurence), 0, annovar_df$noncoding_occurence)

    # keep first symbol if multiple symbols exist
    annovar_df$Gene.refGene[grepl(";", annovar_df$Gene.refGene)] <- vapply(annovar_df$Gene.refGene[grepl(";", annovar_df$Gene.refGene)], function(x) unlist(strsplit(x, ";"))[1], "char")

    # determine genes containing hotspot mutation
    cond <- (annovar_df$coding_occurence > threshold) | (annovar_df$noncoding_occurence > threshold)
    hotspot_genes <- unique(annovar_df$Gene.refGene[cond])

    return(hotspot_genes)
}

#' Determine Double-Hit Genes
#'
#' @inheritParams create_metaprediction_score_df
#' @inheritParams create_gene_level_scna_df
#' @param log2_threshold \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} threshold
#' for determining homozygous loss events (default = -1)
#'
#' @return vector of gene symbols that are subject to double-hit event(s), i.e.
#' non-synonymous mutation + homozygous CN loss
#'
#' @examples
#' path2annovar_csv <- system.file("extdata/imielinski.hg19_multianno.csv",
#'                                 package = "driveR")
#' \dontrun{
#' dhit_genes <- driveR:::determine_double_hit_genes(path2annovar_csv, imielinski_scna_table)
#' }
determine_double_hit_genes <- function(annovar_csv_path,
                                  scna_df,
                                  gene_overlap_threshold = 25,
                                  log2_threshold = -1) {
    ### argument checks
    if (!is.numeric(log2_threshold))
        stop("`log2_threshold` should be numberic")

    ### gene-level hom. loss df
    genes_df <- create_gene_level_scna_df(scna_df = scna_df,
                                          gene_overlap_threshold = gene_overlap_threshold)
    # keep only hom. loss
    loss_genes_df <- genes_df[genes_df$log2ratio < log2_threshold, ]

    ### non-synonymous mutations df
    annovar_df <- utils::read.csv(annovar_csv_path)
    # exclude synonymous muts
    non_syn_df <- annovar_df[annovar_df$ExonicFunc.refGene != "synonymous SNV", ]
    # keep first symbol if multiple symbols exist
    non_syn_df$Gene.refGene[grepl(";", non_syn_df$Gene.refGene)] <- vapply(non_syn_df$Gene.refGene[grepl(";", non_syn_df$Gene.refGene)], function(x) unlist(strsplit(x, ";"))[1], "char")

    ### determine double-hit genes
    tmp <- unique(loss_genes_df$symbol)
    tmp2 <- unique(non_syn_df$Gene.refGene)
    dhit_genes <- tmp[tmp %in% tmp2]

    return(dhit_genes)
}
