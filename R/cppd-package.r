#' cppd.
#'
#' @name cppd
#' @docType package
#' 
#' @import tidyr
#' @import dplyr
#' @import purrr
#' @import purrrlyr
#' @import stringr
#' @import ggplot2
#' @import misha
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom yaml yaml.load_file
#' @importFrom data.table fwrite
#' @importFrom readr write_csv
#' @importFrom readr write_tsv
#' @importFrom GetoptLong qq
#' @importFrom logging loginfo
#' @importFrom logging logwarn
#' @importFrom logging logerror
NULL

utils::suppressForeignCheck(c(".", ".GLIBDIR", "A", "C", "G", "TM_c", "TM_cm", "best_TM_diff", "center", "cg_num", "cgs", "chrom", "chrom1", "chunk", "dist", "e", "end", "end1", "end_reg", "gene", "geneSymbol", "genomic_seq", "homopol", "intervalID", "kmer", "kmer_id", "kmers", "l", "len", "mapping", "max_TM", "max_TM_diff", "min_TM", "min_TM_diff", "min_cgs_minus", "min_cgs_plus", "min_cgs_strand", "n_kmer", "n_reg", "n_reg_org", "new_id", "read", "reg_id", "retv", "revcomp", "revcomp_dup", "s", "seqc", "seqcm", "start", "start1", "start_reg", "strand", "type"))
