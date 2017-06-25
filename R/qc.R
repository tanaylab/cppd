check_sequence <- function(probes){
	pseq <- probes %>% 
		select(chrom, start, end, strand, cg_num, genomic_seq, seq) %>% 
		mutate(gseq = toupper(gseq.extract(.)), gseq_conv = gseq.extract_conv(.)) %>% 
		mutate(cgs = str_count(genomic_seq, 'CG'))
	browser()
	assert_probes(all(pseq$genomic_seq == pseq$gseq), 'genomic_seq', 'genomic_seq is corrupt')	
	assert_probes(all(pseq$cg_num > 0 | pseq$seq == pseq$gseq_conv), 'seq', 'seq is corrupt')	
	assert_probes(all(pseq$cg_num == pseq$cgs), 'cg_num', 'cg_num is corrupt')	
}

assert <- function (cond, error_msg='', ...){
  if (!cond) {
  	logerror(error_msg, ...)  
  	stop(error_msg)  
  }
}

assert_probes <- function(cond, test_msg, error_msg, ...){
	loginfo(sprintf('testing %s', test_msg))
	assert(cond, error_msg, ...)
	loginfo(sprintf('%s is OK', test_msg))
}