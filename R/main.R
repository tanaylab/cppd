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
#'
#' @export
cppd.dump_example_config <- function(path){
    config_files <- dir(system.file("config", package='cppd'), full.names=T)
    dir.create(path, recursive=T, showWarnings=FALSE)
    ret <- file.copy(config_files, path, recursive=T, overwrite=FALSE)
    if (!all(ret)) {
        warning("Couldn't dump config files to ", path, "\n  Perhaps they're already there? ")
    } else {
        message("Dumped config files to: ", path)
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


generate_probes <- function(regions, chunk_size, probes, TM_range, optimal_TM_range, regions_expanded=NULL, candidates=NULL, regions_annot=NULL, n_probes=NULL, max_shared_revcomp=15, downsample=FALSE, threads=1, misha_root=NULL, rm_revcomp=TRUE, ...){
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

    nchunks <- ceiling(nrow(regions) / chunk_size)
    loginfo(qq('number of chunks: @{nchunks}'))  

    cands <- regions %>% mutate(chunk = ntile(1:n(), nchunks)) %>% group_by(chunk) %>% by_slice(~ do.call_ellipsis(get_candidates, list(regs=., TM_range=TM_range, threads=threads), ...))
    
    regs_exp <- cands$.out %>% map_df('regs_exp')

    if (!is.null(regions_expanded)){
        write_csv(regs_exp, regions_expanded)
    }

    cands <- cands$.out %>% map_df('cands')
    if (!is.null(candidates)){
        write_csv(cands, candidates)
    }

    probes <- choose_probes(cands=cands, regions=regions, exp_regions=regs_exp, probes_ofn=probes, regs_annots=regions_annot, n_probes=n_probes, kmer_len=max_shared_revcomp, downsample=downsample, TM_range=optimal_TM_range, rm_revcomp=rm_revcomp)

    if (!keep_field){
        regions <- regions %>% select(-keep)
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

