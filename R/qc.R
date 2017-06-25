# check_probes <- function(probes, regions=NULL, max_dist=350, max_cg_num=2, TM_range=c(60, 72), optimal_TM_range=c(62,65), max_map=1, max_homopol=8, max_revcomp=15, max_kmers=40000){
check_probes <- function(probes, max_dist, max_cg_num, TM_range, optimal_TM_range, max_map, max_homopol, max_revcomp, max_kmers, regions=NULL){
	check_sequence(probes)

	probes_per_reg <- probes %>% count(chrom, start_reg, end_reg) %>% pull(n)
	assert_probes_warn(all(probes_per_reg == 2), 'probes per region', 'not all regions have 2 probes')

	assert_probes(all(probes$start >= probes$start_reg & probes$start <= probes$end_reg), 'probes start coordinates are within the region', 'start coordinates outside of region')
	assert_probes(all(probes$end >= probes$start_reg & probes$end <= probes$end_reg), 'probes end coordinates are within the region', 'end coordinates outside of region')

	# probe requirments
	assert_probes(all(probes$min_TM >= TM_range[1]), 'minimal TM', 'some probes have TM less than %d', TM_range[1])
	assert_probes(all(probes$max_TM <= TM_range[2]), 'maximal TM', 'some probes have TM more than %d', TM_range[2])

	assert_probes_warn(all(probes$min_TM >= optimal_TM_range[1]), 'optimal minimal TM', 'some probes have TM less than %d', optimal_TM_range[1])
	assert_probes_warn(all(probes$max_TM <= optimal_TM_range[2]), 'optimal maximal TM', 'some probes have TM more than %d', optimal_TM_range[2])

	assert_probes(all(probes$revcomp <= max_revcomp), 'maximal shared revcomp', 'some probes have shared revcomp more than %d', max_revcomp)
	assert_probes(all(probes$mapping <= max_map), 'maximal mappings', 'some probes are mapped more than %d times', max_map)
	assert_probes(all(probes$kmers <= max_kmers), 'maximal kmers', 'some probes have more than %d kmers in the genome', max_kmers)
	assert_probes(all(probes$homopol <= max_homopol), 'maximal homoplymers', 'some probes have more than %d homoplymers', max_homopol)

	assert_probes(all(probes$cg_num <= max_cg_num), 'maximal number of CpGs', 'some probes have more than %d cpgs', max_cg_num)

	if (!is.null(regions)){
		reg_ids <-  probes %>% separate_rows(reg_id) %>% pull(reg_id)
		assert_probes(all(reg_ids %in% regions$id), 'region ids appear in original regions', 'some region ids do not appear in original regions')		

		regs_nb <- regions %>% filter(id %in% reg_ids) %>% gintervals.neighbors1(probes %>% select(chrom, start, end, strand), maxneighbors=2) %>% mutate(dist = abs(dist))
		assert_probes_warn(all(regs_nb[['dist']] <= max_dist), 'probes distance', 'some probes are more than %d from the original region', max_dist)		
	}
}


check_sequence <- function(probes){
	pseq <- probes %>% 
		select(chrom, start, end, strand, cg_num, genomic_seq, seq) %>% 
		mutate(gseq = toupper(gseq.extract(.)), gseq_conv = gseq.extract_conv(.)) %>% 
		mutate(cgs = str_count(genomic_seq, 'CG'), seq_cgs = str_count(seq, 'CG'))
	
	assert_probes(all(pseq$genomic_seq == pseq$gseq), 'genomic_seq', 'genomic_seq is corrupt')	
	assert_probes(all(pseq$cg_num > 0 | pseq$seq == pseq$gseq_conv), 'seq', 'seq is corrupt')	
	assert_probes(all(pseq$cg_num == pseq$cgs), 'cg_num', 'cg_num is corrupt')	
	assert_probes(all(pseq$seq_cgs == 0), 'cgs in probe sequence', 'there are CpGs in probe sequence')	
}

assert <- function (cond, error_msg='', ...){
  if (!cond) {
  	logerror(error_msg, ...)  
  	stop(error_msg)  
  }
}

assert_probes <- function(cond, test_msg, error_msg, ...){	
	assert(cond, error_msg, ...)
	loginfo(sprintf('%s... OK', test_msg))
}

assert_probes_warn <- function(cond, test_msg, warning_msg, ...){	
	if (!cond) {
  		logwarn(warning_msg, ...)  
  		warning(warning_msg)  	
  	} else {
  		loginfo(sprintf('%s... OK', test_msg))	
  	}
}