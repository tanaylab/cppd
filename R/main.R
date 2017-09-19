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

    cmd_args <- conf2args(conf, c(generate_probes, get_candidates, regions2seq, check_probes, gcluster.run2))

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
    readr::write_lines(yaml::as.yaml(example), paste0(path, '/', name))

    if (!all(ret)) {
        logwarn("Couldn't dump config files to: %s", path)
    } else {
        loginfo("Dumped config files to: %s", path)
    }
}

generate_probes <- function(regions, chunk_size, probes, TM_range, optimal_TM_range, regions_expanded=NULL, candidates=NULL, chosen_candidates=NULL, regions_annot=NULL, n_probes=NULL, max_shared_revcomp=15, downsample=FALSE, threads=1, misha_root=NULL, rm_revcomp=TRUE, workdir=tempdir(), use_sge=FALSE, verify_probes=TRUE, only_choose=FALSE, ...){

    opt <- getOption('gmax.data.size')
    on.exit(options(gmax.data.size=opt))

    if (!is.null(misha_root)){
        gsetroot(misha_root)
    }

    if (is.character(regions)){
        regions <- fread(regions) %>% tbl_df
    }

    keep_field <- 'keep' %in% colnames(regions)
    if (!keep_field){
        regions[['keep']] <- FALSE
    }

    if (!only_choose){
    
        nchunks <- ceiling(nrow(regions) / chunk_size)
        loginfo(qq('number of chunks: @{nchunks}'))

        regs <- regions %>% mutate(chunk = ntile(1:n(), nchunks))

        run_chunk <- function(chunk_num, temp_prefix, ...){   
            loginfo('chunk %d', chunk_num)   
            cands_ofn <- paste0(temp_prefix, '_chunk_', chunk_num, '_cands')        
            regs_exp_ofn <- paste0(temp_prefix, '_chunk_', chunk_num, '_regs_exp')
            
            do.call_ellipsis(get_candidates, 
                list(regs=regs %>% filter(chunk == chunk_num) %>% select(-chunk), 
                    TM_range=TM_range, 
                    threads=threads, 
                    cands_ofn=cands_ofn, 
                    regs_exp_ofn=regs_exp_ofn), ...)
            cands <- fread(cands_ofn, sep=',') %>% as.tibble()
            regs_exp <- fread(regs_exp_ofn, sep=',') %>% as.tibble()
            chosen_cands <- choose_probes_per_regions(cands=cands, 
                    exp_regions=regs_exp,
                    TM_range=optimal_TM_range)
            return(chosen_cands)
        }

        temp_prefix <- tempfile(tmpdir=workdir)
        on.exit(system(qq('rm -f @{temp_prefix}*')))
        
        if (use_sge){
            cmds <- paste0('run_chunk(', 1:nchunks, ', temp_prefix=temp_prefix, ...)')            
            res <- gcluster.run2(command_list=cmds, ...)
            chosen_cands <- map_df(res, 'retv')
        } else {            
            chosen_cands <- map_dfr(1:nchunks, run_chunk, temp_prefix=temp_prefix, ...)
        }   
        
        loginfo('collecting chunks - expanded regions')
        regs_exp <- map_df(paste0(temp_prefix, '_chunk_', 1:nchunks, '_regs_exp'), ~ fread(.)) %>% as.tibble()
        # regs_exp <- plyr::adply(paste0(temp_prefix, '_chunk_', 1:nchunks, '_regs_exp'), 1, fread, .parallel=TRUE) %>% select(-X1) %>% as.tibble()

        if (!is.null(regions_expanded)){
            fwrite(regs_exp, regions_expanded, sep=',')
        }

        loginfo('collecting chunks - candidates')
        cands <- map_df(paste0(temp_prefix, '_chunk_', 1:nchunks, '_cands'), ~ fread(.)) %>% as.tibble()
        # cands <- plyr::adply(paste0(temp_prefix, '_chunk_', 1:nchunks, '_cands'), 1, fread, .parallel=TRUE) %>% select(-X1) %>% as.tibble()
        
        if (!is.null(candidates)){
            fwrite(cands, candidates, sep=',')
        }

        if (!is.null(chosen_candidates)){
            fwrite(chosen_cands, chosen_candidates, sep=',')   
        }
    } else {
        loginfo('loading candidates and regions')
        chosen_cands <- fread(chosen_candidates) %>% as.tibble()
        regs_exp <- fread(regions_expanded) %>% as.tibble()
    }

    loginfo('calculating multiple-regions statistics')
    probes <- choose_probes(probes=chosen_cands, regions=regions, exp_regions=regs_exp, probes_ofn=probes, regs_annots=regions_annot, n_probes=n_probes, kmer_len=max_shared_revcomp, downsample=downsample, rm_revcomp=rm_revcomp)

    loginfo('checking probes')    
    if (verify_probes){
        do.call_ellipsis(check_probes, list(probes=probes, regions=regions, TM_range=TM_range, optimal_TM_range=optimal_TM_range), ...)
    }

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

