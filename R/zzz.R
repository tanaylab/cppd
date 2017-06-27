.onLoad <- function(libname, pkgname) {    
    cppd.set_parallel(parallel::detectCores() / 2)
    if(getRversion() >= "2.15.1"){
    	utils::globalVariables(c(".", ".GLIBDIR", "A", "C", "G", "TM_c", "TM_cm", "best_TM_diff", "center", "cg_num", "cgs", "chrom", "chrom1", "chunk", "dist", "e", "end", "end1", "end_reg", "gene", "geneSymbol", "genomic_seq", "homopol", "intervalID", "kmer", "kmer_id", "kmers", "l", "len", "mapping", "max_TM", "max_TM_diff", "min_TM", "min_TM_diff", "min_cgs_minus", "min_cgs_plus", "min_cgs_strand", "n_kmer", "n_reg", "n_reg_org", "new_id", "read", "reg_id", "retv", "revcomp", "revcomp_dup", "s", "seqc", "seqcm", "start", "start1", "start_reg", "strand", "type"))	
    }
	
}


