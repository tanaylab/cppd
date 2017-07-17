#' QC test pbat capture probes
#'
#' @param probes data frame with probes (output of cppd.generate_probes)
#' @param conf_fn configuration file
#' @param defaults_fn default configuration file (optional)
#' @param log_fn log file (optional)
#'
#'
#' @export
cppd.check_probes <- function(probes, conf_fn, defaults_fn=NULL, log_fn=NULL){
	logging::basicConfig()

	if (!is.null(log_fn)){
        logging::addHandler(logging::writeToFile, file=log_fn)
    }
    conf <- read_yaml(conf_fn)
    if (!is.null(defaults_fn)){
        defaults <- read_yaml(defaults_fn)
        conf <- plyr::defaults(conf, defaults)
    }
    conf <- apply_genome_conf(conf)    

    cmd_args <- conf2args(conf, c(check_probes))    
    cmd_args[['probes']] <- NULL

    loginfo('testing with the following paramters:')
    walk2(names(cmd_args), cmd_args, function(x, y) loginfo("%s: %s", x, y))

    cmd_args[['probes']] <- probes    
    do.call(check_probes, cmd_args)
}

check_probes <- function(probes, max_dist=350, max_cg_num=2, TM_range=c(60, 72), optimal_TM_range=c(62,65), max_map=1, max_homopol=8, max_revcomp=15, max_kmers=40000, regions=NULL, misha_root=NULL){

    if (!is.null(misha_root)){
        gsetroot(misha_root)
    }

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

	seq_n <- count(probes, seq)
	assert_probes_warn(all(seq_n$n == 1), 'duplicates', 'some sequences appear more than once')	

	probes_per_region <- probes %>% count(chrom, start_reg, end_reg)
	assert_probes_warn(all(probes_per_region$n == 2), 'probes per region', 'some regions have more than 2 probes')	

	if (!is.null(regions)){
		if (is.character(regions)){
			regions <- fread(regions, sep=',') %>% as.tibble()
		}
		if ('reg_id' %in% colnames(probes)){
			reg_ids <-  probes %>% tidyr::separate_rows(reg_id) %>% pull(reg_id)
			assert_probes(all(reg_ids %in% regions$id), 'region ids appear in original regions', 'some region ids do not appear in original regions')		
			
			regs_nb <- regions %>% filter(id %in% reg_ids) %>% gintervals.neighbors1(probes %>% select(chrom, start, end, strand), maxneighbors=2) %>% mutate(dist = abs(dist))
			assert_probes_warn(all(regs_nb[['dist']] <= max_dist), 'probes distance', 'some probes are more than %dbp from the original region', max_dist)		

			regs_dist <- probes %>% select(chrom, start, end) %>% gintervals.neighbors1(regions) %>% pull(dist)		
			assert_probes_warn(all(regs_dist <= max_dist), 'probes distance from any region', 'some probes are more than %dbp from any original region', max_dist)		

		}		
	}

	check_sequence(probes)
}

check_sequence <- function(probes){
	opt <- getOption('gmax.data.size')
	on.exit(options(gmax.data.size=opt))

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