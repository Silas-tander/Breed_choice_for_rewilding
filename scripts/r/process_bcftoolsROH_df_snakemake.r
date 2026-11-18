# Combine the bcftools roh RG lines of chr1-31 into one CSV.
# Usage: Rscript process_bcftoolsROH_df_snakemake.r <roh_dir> <out.csv>

library(readr)
library(dplyr)
library(purrr)

chr_list <- paste0("chr", 1:31)
args <- commandArgs(trailingOnly = TRUE)
base_path <- args[1]
output_path <- args[2]

process_roh_file <- function(chr) {
  # {chr}.roh.txt (mustang_t2t) or {chr}.AF.roh.txt (EquCab3 pipelines)
  file_path <- file.path(base_path, paste0(chr, ".roh.txt"))
  if (!file.exists(file_path)) file_path <- file.path(base_path, paste0(chr, ".AF.roh.txt"))
  lines <- read_lines(file_path)
  rg_lines <- lines[grepl("^RG", lines)]
  if (length(rg_lines) == 0) return(NULL)
  # Sample and Chromosome as character: guessed types turn IDs like AB_17A and "chr31" into NA
  df <- read_table2(paste(rg_lines, collapse = "\n"),
                    col_names = c("Type", "Sample", "Chromosome", "Start", "End",
                                  "Length_bp", "Num_markers", "Quality"),
                    col_types = cols(Sample = col_character(), Type = col_character(),
                                     Chromosome = col_character(),
                                     .default = col_double()))
  return(df)
}

all_rg_df <- map_dfr(chr_list, process_roh_file)

all_rg_df <- all_rg_df %>%
  rename(from = Start, to = End, lengthBps = Length_bp)

write_csv(all_rg_df, output_path)