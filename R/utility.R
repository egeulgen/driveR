#' Create Non-coding Impact Score Data Frame
#'
#' @inheritParams  predict_coding_impact
#'
#' @return data frame of meta-prediction scores containing 2 columns: \describe{
#'   \item{gene_symbol}{HGNC gene symbol}
#'   \item{CADD_phred}{PHRED-scaled CADD score}
#' }
create_noncoding_impact_score_df <- function(annovar_csv_path, na.string = ".") {
    annovar_df <- utils::read.csv(annovar_csv_path)
    noncoding_df <- annovar_df[annovar_df$Func.refGene != "exonic", ]
    noncoding_df <- noncoding_df[noncoding_df$CADD_phred != na.string, ]
    noncoding_df$CADD_phred <- as.numeric(noncoding_df$CADD_phred)
    noncoding_scores_df <- noncoding_df[, c("Gene.refGene", "CADD_phred")]
    colnames(noncoding_scores_df) <- c("gene_symbol", "CADD_score")

    # keep first symbol if multiple symbols exist
    noncoding_scores_df$gene_symbol[grepl(";", noncoding_scores_df$gene_symbol)] <- vapply(noncoding_scores_df$gene_symbol[grepl(";", noncoding_scores_df$gene_symbol)], function(x) unlist(strsplit(x, ";"))[1], "char")

    # keep highest score
    noncoding_scores_df <- noncoding_scores_df[order(noncoding_scores_df$CADD_score, decreasing = TRUE), ]
    noncoding_scores_df <- noncoding_scores_df[!duplicated(noncoding_scores_df$gene_symbol), ]

    return(noncoding_scores_df)
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
    q_hits <- S4Vectors::queryHits(hits)
    S4Vectors::mcols(ranges) <- c(S4Vectors::mcols(ranges),
                                  S4Vectors::mcols(scna_gr[q_hits]))

    genes_df <- data.frame(chr = as.vector(GenomeInfoDb::seqnames(ranges)),
                           as.data.frame(S4Vectors::mcols(ranges)))
    genes_df$segment_start <- GenomicRanges::start(GenomicRanges::ranges(scna_gr[q_hits]))
    genes_df$segment_end <- GenomicRanges::end(GenomicRanges::ranges(scna_gr[q_hits]))
    genes_df$transcript_start <- GenomicRanges::start(ranges)
    genes_df$transcript_end <- GenomicRanges::end(ranges)

    tmp1 <- apply(genes_df[, c("transcript_end", "segment_end")], 1, min)
    tmp2 <- apply(genes_df[, c("segment_start", "transcript_start")], 1, max)
    genes_df$transcript_overlap_percent <- (tmp1 - tmp2) / (genes_df$transcript_end - genes_df$transcript_start) * 100

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
#' @param gene_SCNA_df data frame of gene-level SCNAs (output of \code{\link{create_gene_level_scna_df}})
#' @param log2_ratio_threshold the \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}}
#' ratio threshold for keeping high-confidence SCNA events (default = 0.25)
#' @param MCR_overlap_threshold the percentage threshold for the overlap between
#' a gene and an MCR region (default = 25). This means that if only a gene
#' overlaps an MCR region more than this threshold, the gene is assigned the
#' SCNA density of the MCR
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
create_SCNA_score_df <- function(gene_SCNA_df,
                                 log2_ratio_threshold = 0.25,
                                 MCR_overlap_threshold = 25){
    ### argument checks
    if (!is.numeric(log2_ratio_threshold))
        stop("`log2_ratio_threshold` should be numberic")

    if (!is.numeric(MCR_overlap_threshold))
        stop("`MCR_overlap_threshold` should be numberic")

    if (MCR_overlap_threshold < 0 | MCR_overlap_threshold > 100)
        stop("`MCR_overlap_threshold` should be between 0-100")

    #### determine gene-level log2-ratios
    # discard sex chromosomes
    gene_SCNA_df <- gene_SCNA_df[!gene_SCNA_df$chr %in% c("chrX", "chrY"), ]

    # aggregate as max |log2-ratio| over all values per gene
    agg_ratios <- tapply(gene_SCNA_df$log2ratio, gene_SCNA_df$symbol, function(x) {
        max_val <- max(abs(x))
        if (max_val %in% x)
            return(max_val)
        return(-max_val)
    })

    gene_agg_df <- gene_SCNA_df[, c("symbol", "chr", "transcript_start", "transcript_end")]
    gene_agg_df <- gene_agg_df[!duplicated(gene_agg_df$symbol), ]
    gene_agg_df$agg_log2_ratio <- agg_ratios[match(gene_agg_df$symbol, names(agg_ratios))]

    gene_agg_df <- gene_agg_df[!is.na(gene_agg_df$agg_log2_ratio), ]

    # filter for high confidence
    gene_agg_df <- gene_agg_df[abs(gene_agg_df$agg_log2_ratio) > log2_ratio_threshold, ]

    # return empty data frame if no gene passing log2 threshold
    if (nrow(gene_agg_df) == 0) {
        df <- data.frame(symbol = character(), SCNA_density = numeric())
        warning("No gene-level SCNA event passed the log2_ratio_threshold")
        return(df)
    }

    ### assign MCR score to overlapping genes
    MCR_gr <- GenomicRanges::makeGRangesFromDataFrame(MCR_table, keep.extra.columns = TRUE)
    agg_gr <- GenomicRanges::makeGRangesFromDataFrame(gene_agg_df,
                                                      seqnames.field = "chr",
                                                      start.field = "transcript_start",
                                                      end.field = "transcript_end",
                                                      keep.extra.columns = TRUE)
    # overlaps
    hits <- GenomicRanges::findOverlaps(agg_gr,
                                        MCR_gr,
                                        type="any", select="all")

    # return empty data frame if no overlap
    if (length(hits) == 0) {
        df <- data.frame(symbol = character(), SCNA_density = numeric())
        warning("No gene-level SCNA event overlapped any MCR region")
        return(df)
    }

    # obtain overlap data frame between MCR regions and gene-level SCNA
    ranges <- MCR_gr[S4Vectors::subjectHits(hits)]
    q_hits <- S4Vectors::queryHits(hits)
    S4Vectors::mcols(ranges) <- c(S4Vectors::mcols(ranges), S4Vectors::mcols(agg_gr[q_hits]))

    final_scna_df <- data.frame(chr = as.vector(GenomeInfoDb::seqnames(ranges)),
                                as.data.frame(S4Vectors::mcols(ranges)))

    final_scna_df$transcript_start <- GenomicRanges::start(GenomicRanges::ranges(agg_gr[q_hits]))
    final_scna_df$transcript_end <- GenomicRanges::end(GenomicRanges::ranges(agg_gr[q_hits]))
    final_scna_df$MCR_start <- GenomicRanges::start(ranges)
    final_scna_df$MCR_end <- GenomicRanges::end(ranges)

    tmp1 <- apply(final_scna_df[, c("MCR_end", "transcript_end")], 1, min)
    tmp2 <- apply(final_scna_df[, c("transcript_start", "MCR_start")], 1, max)
    final_scna_df$MCR_overlap_percent <- (tmp1 - tmp2) / (final_scna_df$MCR_end - final_scna_df$MCR_start) * 100

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
#' @inheritParams predict_coding_impact
#' @param hotspot_threshold to determine hotspot genes, the (integer) threshold
#' for the minimum number of cases with certain mutation in COSMIC (default = 5)
#'
#' @return vector of gene symbols of genes containing hotspot mutation(s)
determine_hotspot_genes <- function(annovar_csv_path, hotspot_threshold = 5L) {
    # argument check
    if (!is.numeric(hotspot_threshold))
        stop("`hotspot_threshold` should be numeric")

    annovar_df <- utils::read.csv(annovar_csv_path)

    # parse COSMIC occurences
    coding_column <- colnames(annovar_df)[grepl("cosmic\\d+_coding", colnames(annovar_df))]
    non_coding_column <- colnames(annovar_df)[grepl("cosmic\\d+_noncoding", colnames(annovar_df))]

    annovar_df$coding_occurence <- vapply(annovar_df[, coding_column],
                                          function(x) {
                                              tmp <- unlist(strsplit(x, ";"))[2]
                                              tmp <- sub("OCCURENCE=", "", tmp)
                                              tmp <- unlist(strsplit(tmp, ","))
                                              tmp <- as.numeric(sub("\\(.*\\)", "", tmp))
                                              return(sum(tmp))
                                          }, 1)

    annovar_df$noncoding_occurence <- vapply(annovar_df[, non_coding_column],
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
    cond <- (annovar_df$coding_occurence > hotspot_threshold) | (annovar_df$noncoding_occurence > hotspot_threshold)
    hotspot_genes <- unique(annovar_df$Gene.refGene[cond])

    return(hotspot_genes)
}

#' Determine Double-Hit Genes
#'
#' @inheritParams predict_coding_impact
#' @inheritParams create_SCNA_score_df
#' @param log2_hom_loss_threshold to determine double-hit events, the
#' \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} threshold for identifying
#' homozygous loss events (default = -1).
#' @param batch_analysis boolean to indicate whether to perform batch analysis
#' (\code{TRUE}, default) or personalized analysis (\code{FALSE}). If \code{TRUE},
#' a column named 'tumor_id' should be present in both the ANNOVAR csv and the SCNA
#' table.
#'
#' @return vector of gene symbols that are subject to double-hit event(s), i.e.
#' non-synonymous mutation + homozygous copy-number loss
determine_double_hit_genes <- function(annovar_csv_path,
                                       gene_SCNA_df,
                                       log2_hom_loss_threshold = -1,
                                       batch_analysis = FALSE) {
    ### argument checks
    if (!is.numeric(log2_hom_loss_threshold))
        stop("`log2_hom_loss_threshold` should be numberic")

    if (!is.logical(batch_analysis))
        stop("`batch_analysis` should be `TRUE` or `FALSE`")

    ### gene-level hom. loss df
    # keep only hom. loss
    loss_genes_df <- gene_SCNA_df[gene_SCNA_df$log2ratio < log2_hom_loss_threshold, ]
    # discard sex chromosomes
    loss_genes_df <- loss_genes_df[!loss_genes_df$chr %in% c("chrX", "chrY"), ]

    ### non-synonymous mutations df
    annovar_df <- utils::read.csv(annovar_csv_path)
    # exclude synonymous mutations
    non_syn_df <- annovar_df[annovar_df$ExonicFunc.refGene != "synonymous SNV", ]
    # keep first symbol if multiple symbols exist
    non_syn_df$Gene.refGene[grepl(";", non_syn_df$Gene.refGene)] <- vapply(non_syn_df$Gene.refGene[grepl(";", non_syn_df$Gene.refGene)], function(x) unlist(strsplit(x, ";"))[1], "char")

    ### determine double-hit genes
    if (batch_analysis) {
        if (!("tumor_id" %in% colnames(non_syn_df) & "tumor_id" %in% colnames(loss_genes_df)))
            stop("'tumor id' should be present in both ANNOVAR output and SCNA table if `batch_analysis == TRUE`")

        all_donors <- unique(non_syn_df$tumor_id)
        all_donors <- all_donors[all_donors %in% loss_genes_df$tumor_id]

        dhit_genes <- c()
        for (donor in all_donors) {
            tmp_loss <- loss_genes_df[loss_genes_df$tumor_id == donor, ]
            tmp_loss <- unique(tmp_loss$symbol)

            tmp_nonsyn <- non_syn_df$Gene.refGene[non_syn_df$tumor_id == donor]

            dhit <- intersect(tmp_loss, tmp_nonsyn)
            if (length(dhit) != 0)
                dhit_genes <- c(dhit_genes, dhit)
        }
        dhit_genes <- unique(dhit_genes)
    } else {
        tmp_loss <- unique(loss_genes_df$symbol)
        tmp_nonsyn <- unique(non_syn_df$Gene.refGene)
        dhit_genes <- intersect(tmp_loss, tmp_nonsyn)
    }
    return(dhit_genes)
}
