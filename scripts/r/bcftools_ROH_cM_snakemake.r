#!/usr/bin/env Rscript
# ROH lengths in cM from the averaged recombination maps, plus FROH and length classes.
# Usage: Rscript bcftools_ROH_cM_snakemake.r <roh_csv> <recomb_dir> <output_csv> [genome_len_bp]
# genome_len_bp defaults to EquCab3 chr1-31.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript scripts/compute_cM_from_roh.R <roh_csv> <recomb_dir> <output_csv>")
}
roh_csv     <- args[1]
recomb_dir  <- args[2]
output_csv  <- args[3]
genome_len_arg <- if (length(args) >= 4) as.numeric(args[4]) else NA

# EquCab3.0 chr1-31 lengths (EquCab3.0_chrom.fa.fai)
genome_len <- sum(c(
  188260577, 121350024, 121351753, 109462549, 96759418, 87230776, 100787686,
  97563019, 85793548, 85155674, 61676917, 36992759, 43784481, 94600235,
  92851403, 88962352, 80722430, 82641348, 62681739, 65343332, 58984458,
  50928189, 55556184, 48288683, 40282968, 43147642, 40254690, 47348498,
  34776120, 31395959, 26001039
))
if (!is.na(genome_len_arg)) genome_len <- genome_len_arg

read_recomb_map <- function(chr, dir_path) {
  file_name <- sprintf("Averaged_ECA%s_map.txt", chr)
  file_path <- file.path(dir_path, file_name)
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$chr <- as.integer(chr)
  df
}

compute_cM <- function(chrom, query_from, query_to, recom_all) {
  if (query_from > query_to) {
    tmp <- query_from
    query_from <- query_to
    query_to <- tmp
  }

  chr_recom <- recom_all %>% filter(chr == chrom)
  if (nrow(chr_recom) == 0) return(NA_real_)

  overlaps <- chr_recom %>%
    filter(from < query_to & to > query_from)

  if (nrow(overlaps) == 0) return(0)

  # (bp overlap) * (cM per Mb) / 1e6
  contrib <- ((pmin(overlaps$to, query_to) - pmax(overlaps$from, query_from)) * overlaps$`Rate.cM.Mb.`) / 1e6
  sum(contrib, na.rm = TRUE)
}

main <- function() {
  chr_combined_bcf <- read.csv(roh_csv, stringsAsFactors = FALSE)

  chr_combined_bcf <- chr_combined_bcf %>%
    dplyr::rename(
      from      = from,
      to        = to,
      lengthBps = lengthBps,
      ind       = Sample,
      chrom     = Chromosome
    ) %>%
    # maps are keyed by chromosome number; accept "chr31" as well as "31"
    mutate(chrom = as.integer(sub("^chr", "", as.character(chrom)))) %>%
    mutate(
      length_interval = dplyr::case_when(
        lengthBps < 1e6 ~ "< 1 Mb",
        lengthBps >= 1e6  & lengthBps < 2e6  ~ "1-2 Mb",
        lengthBps >= 2e6  & lengthBps < 4e6  ~ "2-4 Mb",
        lengthBps >= 4e6  & lengthBps < 8e6  ~ "4-8 Mb",
        lengthBps >= 8e6  & lengthBps < 1.6e7 ~ "8-16 Mb",
        lengthBps >= 1.6e7 ~ ">16 Mb"
      ),
      log_lengthBps = log(lengthBps),
      FROH = lengthBps / genome_len
    )

  chr_list_recomb_bcf <- 1:31
  recom_bcf <- bind_rows(lapply(chr_list_recomb_bcf, read_recomb_map, dir_path = recomb_dir))

  recom_all <- recom_bcf %>%
    arrange(chr, `Position.bp.`) %>%
    mutate(
      from = `Position.bp.`,
      to = ifelse(
        dplyr::lead(chr) == chr,
        dplyr::lead(`Position.bp.`),
        # last entry on a chromosome: estimate the end using cM increment / rate
        `Position.bp.` + ((`Map.cM.` - dplyr::lag(`Map.cM.`, default = 0)) / `Rate.cM.Mb.`) * 1e6
      )
    )

  chr_FROHints_bcf_cM <- chr_combined_bcf %>%
    rowwise() %>%
    mutate(cM = compute_cM(chrom, from, to, recom_all)) %>%
    ungroup() %>%
    mutate(
      length_interval_cM = case_when(
        cM >= 1  & cM < 2  ~ "1-2 cM",
        cM >= 2  & cM < 4  ~ "2-4 cM",
        cM >= 4  & cM < 8  ~ "4-8 cM",
        cM >= 8  & cM < 16 ~ "8-16 cM",
        cM >= 16          ~ ">16 cM",
        TRUE              ~ NA_character_
      )
    )

  write.csv(chr_FROHints_bcf_cM, file = output_csv, row.names = FALSE, quote = FALSE)
}

main()
