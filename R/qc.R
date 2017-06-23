check_sequence <- function(probes){
	pseq <- probes %>% select(chrom, start, end, strand, cg_num, genomic_seq, seq) %>% mutate(gseq = toupper(gseq.extract(.)), gseq_conv = gseq.extract_conv(.))	
	browser()
	assert_probes(all(pseq$genomic_seq == pseq$gseq), 'genomic_seq', 'genomic_seq is corrupt')
	assert_probes(all(pseq$cg_num > 0 | pseq$seq == pseq$gseq_conv), 'seq', 'seq is corrupt')
}

assert <- function (cond, ...){
  if (!cond) {
  	logerror(...)    
  }
}

assert_probes <- function(cond, test_msg, ...){
	loginfo(sprintf('testing %s', test_msg))
	assert(cond, ...)
}