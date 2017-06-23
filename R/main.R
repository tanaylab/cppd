#' Generate pbat capture probes
#'
#' @param conf_fn configuration file
#' @param defaults_fn default configuration file (optional)
#' @param log_fn log file (optional)
#' @param return_probes return a data frame with the probes
#'
#' @return probes data frame if return_probes is TRUE, and NULL otherwise
#'
#' @export
cppd.generate_probes <- function(conf_fn, defaults_fn=NULL, log_fn=NULL, return_probes=FALSE){
    if (!is.null(log_fn)){
        logging::addHandler(logging::writeToFile, file=log_fn)
    }
    conf <- read_yaml(conf_fn)
    if (!is.null(defaults_fn)){
        defaults <- read_yaml(defaults_fn)
        conf <- plyr::defaults(conf, defaults)
    }
    conf <- apply_genome_conf(conf)

    cmd_args <- conf2args(conf, c(generate_probes, get_candidates, regions2seq))

    loginfo('run with the following paramters:')
    walk2(names(cmd_args), cmd_args, function(x, y) loginfo("%s: %s", x, y))

    probes <- do.call(generate_probes, cmd_args)

    if (return_probes){
        return(probes)
    }
    return(NULL)
}

#' Dump config file templates.
#'
#' Dump templates of config files required by the package.
#'
#' @param path directory to dump files to.
#' @param name name of yaml config file
#' @param probes_dir (name of the probes dir)
#' @param genome genome id
#' @param n_probes number of desired probes
#' @param downsample downsample to \code{n_probes} probes
#' @param insert_length insert length of the probes
#' @param annotations use regions annotations
#'
#' @export
cppd.dump_example_config <- function(path,
                                     name='example.yaml',
                                     probes_dir=NULL,
                                     genome=NULL,
                                     n_probes=20000,
                                     downsample=FALSE,
                                     insert_length=350,
                                     annotations=FALSE){
    config_files <- dir(system.file("config", package='cppd'), full.names=T)
    defaults_conf <- system.file('config/defaults.yaml', package='cppd')
    dir.create(path, recursive=T, showWarnings=FALSE)
    ret <- file.copy(defaults_conf, path, recursive=T, overwrite=TRUE)

    example <- yaml::yaml.load_file(system.file('config/example.yaml', package='cppd'))
    if (!is.null(probes_dir)){
        example[['probes_dir']] <- probes_dir
    }
    if (!is.null(genome)){
        example[['genome']] <- genome
    }
    example[['n_probes']] <- n_probes
    example[['downsample']] <- downsample
    example[['insert_length']] <- insert_length
    if (!annotations){
        example[['regions_annot']] <- NULL
    }
    write_lines(yaml::as.yaml(example), paste0(path, '/', name))

    if (!all(ret)) {
        logwarn("Couldn't dump config files to: %s", path)
    } else {
        loginfo("Dumped config files to: %s", path)
    }
}


#' @export
cppd.choose_from_cands <- function(conf_fn, defaults_fn=NULL, log_fn=NULL, return_probes=FALSE){
    if (!is.null(log_fn)){
        logging::addHandler(logging::writeToFile, file=log_fn)
    }
    conf <- read_yaml(conf_fn)
    if (!is.null(defaults_fn)){
        defaults <- read_yaml(defaults_fn)
        conf <- plyr::defaults(conf, defaults)
    }
    conf <- apply_genome_conf(conf)

    regions <- fread(conf$regions)
    cands <- fread(conf$candidates)
    exp_regions <- fread(conf$regions_expanded)
    if (is.null(conf$rm_revcomp)){
        conf$rm_revcomp <- TRUE
    }
    probes <- choose_probes(cands=cands, regions=regions, exp_regions=exp_regions, TM_range=conf$optimal_TM_range, probes_ofn=conf$probes, all_probes_ofn=conf$all_probes_ofn, regs_annots=conf$regions_annot, n_probes=conf$n_probes, kmer_len=conf$max_shared_revcomp, downsample=conf$downsample, rm_revcomp=conf$rm_revcomp)

    if (return_probes){
        return(probes)
    }
    return(NULL)
}


generate_probes <- function(regions, chunk_size, probes, TM_range, optimal_TM_range, regions_expanded=NULL, candidates=NULL, regions_annot=NULL, n_probes=NULL, max_shared_revcomp=15, downsample=FALSE, threads=1, misha_root=NULL, rm_revcomp=TRUE, workdir=tempdir(), use_sge_cands=NULL, use_sge_choose=NULL, use_sge=FALSE, ...){
    if (!is.null(misha_root)){
        gsetroot(misha_root)
    }

    use_sge_cands <- use_sge_cands %||% use_sge
    use_sge_choose <- use_sge_choose %||% use_sge

    if (is.character(regions)){
        regions <- fread(regions) %>% tbl_df
    }

    keep_field <- 'keep' %in% colnames(regions)
    if (!keep_field){
        regions[['keep']] <- FALSE
    }
    

    nchunks <- ceiling(nrow(regions) / chunk_size)
    loginfo(qq('number of chunks: @{nchunks}'))

    regs <- regions %>% mutate(chunk = ntile(1:n(), nchunks))
    run_chunk <- function(chunk_num, ...){   
        loginfo('chunk %d', chunk_num)     
        do.call_ellipsis(get_candidates, 
            list(regs=regs %>% filter(chunk == chunk_num) %>% select(-chunk), 
                TM_range=TM_range, 
                threads=threads, 
                cands_ofn=paste0(tempfile(tmpdir=workdir), '_chunk_', chunk_num, '_cands'), 
                regs_exp_ofn=paste0(tempfile(tmpdir=workdir), '_chunk_', chunk_num, '_regs_exp')), ...)
    }

    if (use_sge_cands){
        cmds <- paste0('run_chunk(', 1:nchunks, ', ...)')            
        res <- gcluster.run2(command_list=cmds, ...)
    } else {
        walk(1:nchunks, ~ run_chunk(.x, ...))
    }   
    
    loginfo('collecting chunks expanded regions')
    regs_exp <- map_df(list.files(workdir, pattern='.+_chunk_\\d+_regs_exp', full.names=TRUE), ~ fread(.)) %>% as.tibble()

    if (!is.null(regions_expanded)){
        fwrite(regs_exp, regions_expanded, sep=',')
    }

    loginfo('collecting chunks candidates')
    cands <- map_df(list.files(workdir, pattern='.+_chunk_\\d+_cands', full.names=TRUE), ~ fread(.)) %>% as.tibble()   
    if (!is.null(candidates)){
        fwrite(cands, candidates, sep=',')
    }
    
    cands <- cands %>% left_join(cands %>% distinct(chrom, start_reg, end_reg) %>% mutate(chunk = ntile(1:n(), nchunks)), by=c('chrom', 'start_reg', 'end_reg'))

    choose_chunk <- function(chunk_num){
        loginfo('chunk %d', chunk_num)     
        choose_probes_per_regions(cands=cands %>% filter(chunk == chunk_num) %>% select(-chunk), 
                exp_regions=regs_exp,
                TM_range=optimal_TM_range)       
    }

    loginfo('choosing candidates per chunk. # of chunks: %d', nchunks)
    if (use_sge_choose){
        cmds <- paste0('choose_chunk(', 1:nchunks, ')')            
        chosen_cands <- gcluster.run2(command_list=cmds, ...) %>% map_df('retv')
    } else {
        # chosen_cands <- map_df(1:nchunks, ~ choose_chunk(.x))
        chosen_cands <- plyr::adply(1:nchunks, 1, function(.x) choose_chunk(.x), .parallel = TRUE)
    }       

    loginfo('calculating multiple-regions statistics')
    probes <- choose_probes(probes=chosen_cands, regions=regions, exp_regions=regs_exp, probes_ofn=probes, regs_annots=regions_annot, n_probes=n_probes, kmer_len=max_shared_revcomp, downsample=downsample, rm_revcomp=rm_revcomp)

    return(probes)
}


############################################################################
apply_genome_conf <- function(conf){
    if ('genome' %in% names(conf) && 'genomes' %in% names(conf)){
        conf <- plyr::defaults(conf, conf$genomes[[conf$genome]])
    }
    return(conf)
}

func_args <- function(f){
    names(as.list(args(f)))
}

conf2args <- function(conf, funcs){
    arg_names <- funcs %>%
        map(func_args) %>%
        as_vector %>%
        unique %>%
        discard(~ .x == '' || .x == '...')
    f_args <- list()
    for (f in names(conf)){
        if (f %in% arg_names){
            f_args[[f]] <- conf[[f]]
        }
    }
    return(f_args)
}

