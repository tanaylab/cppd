fread <- function(...) data.table::fread(..., data.table=FALSE)

.random_track_name <- function(len=12){
    stringi::stri_rand_strings(1, len, "[a-z]")
}

#' Finds neighbors between two sets of intervals (and does not return conflicting column names)
#'
#' @inheritParams misha::gintervals.neighbors
#'
#' @export
#'
gintervals.neighbors1 <- function(intervals1 = NULL,
                                  intervals2 = NULL,
                                  maxneighbors = 1,
                                  mindist = -1e+09,
                                  maxdist = 1e+09,
                                  na.if.notfound = TRUE) {

    res <-
        gintervals.neighbors(
            intervals1 = intervals1,
            intervals2 = intervals2,
            maxneighbors = maxneighbors,
            mindist = mindist,
            maxdist = maxdist,
            na.if.notfound = na.if.notfound
        ) %>%
        tibble::repair_names()

    return(res %>% as_tibble())
}


#' Filter intervals set by distance to another.
#' A wrapper around gintervals.neighbours and filter(dist <= max_distance)
#'
#' @param intervals1 intervals set
#' @param intervals2 intervals set by which to filter intervals1
#' @param max_distance maximal distance of every interval in intervals1 from intervals2 (defualt 0)
#' @param abs_dist take the absolute distance
#' @param bind_intervals2 cbind add intervals2 to result
#' @param ... additional parameters to gintervals.neighbours1
#'
#' @export
#'
gintervals.filter <- function(intervals1, intervals2, max_distance=0, abs_dist = TRUE, bind_intervals2 = FALSE, ...){
    intervals1_cols <- colnames(intervals1)
    res <- intervals1 %>% gintervals.neighbors1(intervals2, ...)
    if (abs_dist){
        res$dist <- abs(res$dist)
    }
    res <- res %>% filter(dist <= max_distance)
    if (!bind_intervals2){
        res <- res %>% select(one_of(intervals1_cols))
    }
    return(res)
}


gintervals.neighbors_y <- function(intervals1, intervals2, maxneighbors=1, mindist=-1e+09, maxdist=1e+09, na.if.notfound=FALSE)
{
    cols1 <- colnames(intervals1)
    cols2 <- colnames(intervals2)
    collisions <- cols2 %in% cols1
    cols2[collisions] <- paste0(cols2[collisions], '.y')
    cols12 <- c(cols1, cols2)

    if (nrow(intervals2) > 0) {
        neighbors <- gintervals.neighbors(intervals1, intervals2, maxneighbors=maxneighbors,
                                          mindist=mindist, maxdist=maxdist, na.if.notfound=na.if.notfound)
        colnames(neighbors)[1:length(cols12)] <- cols12
    }
    else {
        neighbors <- bind_cols(intervals1, data.frame(matrix(NA, nrow=nrow(intervals1), ncol=ncol(intervals2)+1)))
        colnames(neighbors) <- c(cols12, 'dist')
    }

    return(neighbors)
}


gintervals.left_join <- function(intervals1, intervals2, maxneighbors=1, mindist=-1e+09, maxdist=1e+09)
{
    joined <- gintervals.neighbors_y(intervals1, intervals2, maxneighbors=maxneighbors, mindist=mindist, maxdist=maxdist, na.if.notfound=TRUE)
    joined <- joined %>%
              group_by(chrom, start, end) %>%
              filter(is.na(dist) | (abs(dist) == min(abs(dist)))) %>%
              ungroup()

    return(joined)
}


gintervals.centers <- function(inv) {
    inv %>%
        mutate(center = floor((start + end) / 2), start=center, end=center+1) %>%
        select(-center)
}


gintervals.expand <- function(inv, expansion = 100) {
    if (is.character(inv)){
        inv <- gintervals.load(inv)
    }
    inv %>%
        mutate(start = start - expansion, end = end + expansion) %>%
        as.data.frame %>%
        gintervals.force_range()
}

gintervals.normalize <- function(inv, size) {
    centers <- gintervals.centers(inv) %>%
        mutate(end = end - 1)
    return(gintervals.expand(centers, floor(size/2)))
}


gintervals.distance <- function(start1,end1,start2,end2) {
    pmax(pmax(start1, start2) - pmin(end1, end2), 0)
}


gintervals.center <- function(intervals1, intervals2, max_dist=0, size=NULL){
    res <- intervals1 %>% 
        gintervals.neighbors1(intervals2, na.if.notfound=T) %>% 
        filter(abs(dist) <= max_dist) %>% 
        select(chrom=chrom1, start=start1, end=end1, chrom_orig=chrom, start_orig=start, end_orig=end)
    if (!is.null(size)){
        res <- res %>% gintervals.normalize(size) %>% as_tibble()
    }
    return(res)    
}

#' Creates a virtual track and runs gextract
#'
#' @param tracks tracks
#' @param intervals intervals
#' @param colnames colnames
#' @param iterator iterator
#' @param band band
#' @param file file
#' @param intervals.set.out intervals.set.out
#' @param func func
#' @param params params
#'
#' @export
gvextract <- function(tracks, intervals, colnames = NULL, iterator = NULL,
              band = NULL, file = NULL, intervals.set.out = NULL, func=NULL, params=NULL){
    vtracks_pref <- .random_track_name()
    vtracks <- paste0(vtracks_pref, '_', 1:length(tracks))

    walk2(vtracks, tracks, gvtrack.create, func=func, params=params)

    if (is.null(colnames)){
        colnames <- tracks
    }

    res <- gextract(vtracks,
                    intervals = intervals,
                    iterator = iterator,
                    colnames=colnames,
                    band=band,
                    file=file,
                    intervals.set.out=intervals.set.out)

    walk(vtracks, gvtrack.rm)

    return(res)
}



#' Returns the result of track expressions evaluation for each of the
#' iterator intervals, and cbinds the intervals (instead of intervalID)
#'
#' @inheritParams misha::gextract
#' @param suffix suffix for conflicting column names
#'
#' @export
#'
#' @seealso \link[misha]{gextract}
gextract.left_join <- function(expr, intervals = NULL, colnames = NULL, iterator = NULL, band = NULL, file = NULL, intervals.set.out = NULL, suffix='1'){
    if ('character' %in% class(intervals)){
        intervals <- gintervals.load(intervals)
    }
    d <- gextract(expr, intervals = intervals, colnames = colnames, iterator = iterator, band = band, file = file, intervals.set.out = intervals.set.out)
    conflict_names <- which(colnames(intervals) %in% colnames(d))
    colnames(intervals)[conflict_names] <- paste0(colnames(intervals)[conflict_names], suffix)
    intervals$intervalID <- 1:nrow(intervals)
    d <- d %>% arrange(intervalID) %>% left_join(intervals, by='intervalID') %>% select(-intervalID)
    return(d)
}

gseq.extract_conv <- function(..., methylated=TRUE) {
    res <- toupper(gseq.extract(...))
    if (methylated) {
        return(gsub('C(?!G)', 'T', res, perl=T))
    } else {
        return(gsub('C', 'T', res))
    }
}


gseq.rev_comp <- function(s) {
    chartr('acgtACGT', 'tgcaTGCA', s) %>% stringi::stri_reverse()
}

gintervals.mark_overlapping <- function(intervals, unify_touching_intervals=TRUE, var='overlap')
{
    canonic <- intervals %>%
               gintervals.force_range() %>%
               gintervals.canonic(unify_touching_intervals)
    return(intervals %>% mutate(overlap=attr(canonic, 'mapping')))
}


#' Runs R commands on a cluster that supports SGE
#'
#' @inheritParams misha::gcluster.run
#' @param command_list list of strings with R commands
#' @param packages packages to load for each command
#' @param jobs_title title of job names. names would have the title followed
#' by a serial number
#' @param job_names vector with the names of the jobs
#' @param collapse_results collapse return values of the jobs to a data frame.
#' if not possible - would return the usual list.
#'
#' @param queue queue to use
#' @param memory memory requirments (would be called using \code{memory_flag})
#' @param threads threads requirments (would be called using \code{threads_flag})
#' @param io_saturation io_saturation requirments (would be called using \code{io_saturation_flag})
#' @param queue_flag flag for queue (formatted as in qq)
#' @param memory_flag flag for memory requirment (formatted as in qq)
#' @param threads_flag flag for threads requirment (formatted as in qq)
#' @param io_saturation_flag flag for io_saturation requirment (formatted as in qq)
#' @param script sgjob.sh script to use
#'
#' @return if collapse_results is TRUE: data frame with the results of all jobs (rbinded).
#' if collapse_results is FALSE returns the same as: \link[misha]{gcluster.run}
#'
#' @export
#' @seealso \link[misha]{gcluster.run}
#'
gcluster.run2 <- function (...,
                           command_list = NULL,
                           opt.flags = "",
                           max.jobs = 400,
                           debug = FALSE,
                           R = paste0(R.home(component='bin'), '/R'),
                           packages = NULL,
                           jobs_title = NULL,
                           job_names = NULL,
                           collapse_results = FALSE,
                           queue = NULL,
                           memory = NULL,
                           threads = NULL,
                           io_saturation = NULL,
                           queue_flag = '-q @{queue}',
                           memory_flag = '-l mem_free=@{memory}G',
                           threads_flag = '-pe threads @{threads}',
                           io_saturation_flag = '-l io_saturation=@{io_saturation}',
                           script = system.file("cluster", "sgjob.sh", package="cppd")){
    
    if (!is.null(command_list)){ 
        commands <- purrr::map(command_list, function(x) parse(text=x))
    } else {
        commands <- as.list(substitute(list(...))[-1L])
    }

    if (!is.null(queue)){
        opt.flags <- paste(opt.flags, qq(queue_flag))
    }

    if (!is.null(memory)){
        opt.flags <- paste(opt.flags, qq(memory_flag))
    }
    if (!is.null(threads)){
        opt.flags <- paste(opt.flags, qq(threads_flag))
    }
    if (!is.null(io_saturation)){
        opt.flags <- paste(opt.flags, qq(io_saturation_flag))
    }

    if (length(commands) < 1)
        stop("Usage: gcluster.run2(..., command_list = NULL, opt.flags = \"\" max.jobs = 400, debug = FALSE)",
            call. = F)
    if (!length(system("which qsub", ignore.stderr = T, intern = T)))
        stop("gcluster.run2 must run on a host that supports Sun Grid Engine (qsub)",
            call. = F)
    .gcheckroot()
    tmp.dirname <- ""
    submitted.jobs <- c()
    tryCatch({
        tmp.dirname <- tempfile(pattern = "", tmpdir = paste(get("GROOT"),
            "/tmp", sep = ""))
        if (!dir.create(tmp.dirname, recursive = T, mode = "0777"))
            stop(sprintf("Failed to create a directory %s", tmp.dirname),
                call. = F)
        cat("Preparing for distribution...\n")
        save(.GLIBDIR, file = paste(tmp.dirname, "libdir", sep = "/"))
        vars <- ls(all.names = TRUE, envir = parent.frame())
        envir <- parent.frame()
        while (!identical(envir, .GlobalEnv)) {
            envir <- parent.env(envir)
            if (!isNamespace(envir)) {                
                vars <- union(vars, ls(all.names = TRUE, envir = envir))
            }
        }

        suppressWarnings(save(list = vars, file = paste(tmp.dirname, "envir",
            sep = "/"), envir = parent.frame()))
        .GSGECMD <- commands
        save(.GSGECMD, file = paste(tmp.dirname, "commands",
            sep = "/"))
        opts <- options()
        save(opts, file = paste(tmp.dirname, "opts", sep = "/"))
        if (!is.null(packages)){
            .GPACKAGES <- as.list(packages)
        } else {
            .GPACKAGES <- as.list(.packages())    
        }        
        save(.GPACKAGES, file = paste(tmp.dirname, "packages", sep = "/"))

        cat("Running the commands...\n")
        completed.jobs <- c()
        progress <- -1
        repeat {
            num.running.jobs <- length(submitted.jobs) - length(completed.jobs)
            if (length(submitted.jobs) < length(commands) &&
                num.running.jobs < max.jobs) {
                istart <- length(submitted.jobs) + 1
                iend <- min(length(commands), istart + (max.jobs -
                  num.running.jobs) - 1)
                for (i in istart:iend) {
                  out.file <- sprintf("%s/%d.out", tmp.dirname,
                    i)
                  err.file <- sprintf("%s/%d.err", tmp.dirname,
                    i)                  
                    if (!is.null(job_names)){
                        job.name <-job_names[i]
                    } else if (!is.null(jobs_title)){
                        job.name <- sprintf('%s_%s', jobs_title, i)
                    } else {
                        job.name <- sprintf('sgjob_%s', i)
                    }
                  command <- sprintf("qsub -terse -cwd -S /bin/bash -N %s -o %s -e %s -V %s %s %d '%s' '%s'",
                    job.name, out.file, err.file, opt.flags, script, i,
                    tmp.dirname, R)
                  jobid <- system(command, intern = TRUE)
                  if (length(jobid) != 1)
                    stop("Failed to run qsub", call. = FALSE)
                  if (debug)
                    cat(sprintf("\tSubmitted job %d (id: %s)\n",
                      i, jobid))
                  submitted.jobs <- c(submitted.jobs, jobid)
                }
            }
            Sys.sleep(3)
            running.jobs <- .gcluster.running.jobs(submitted.jobs)
            old.completed.jobs <- completed.jobs
            completed.jobs <- setdiff(submitted.jobs, running.jobs)
            if (debug) {
                delta.jobs <- setdiff(completed.jobs, old.completed.jobs)
                if (length(delta.jobs) > 0) {
                  for (jobid in delta.jobs) cat(sprintf("\tJob %d (id: %s) completed\n",
                    match(jobid, submitted.jobs), jobid))
                }
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- length(completed.jobs)
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("\t%d job(s) still in progress\n",
                    length(commands) - progress))
                }
            }
            else {
                if (!length(running.jobs) && length(submitted.jobs) ==
                  length(commands))
                  break
                new.progress <- as.integer(100 * length(completed.jobs)/length(commands))
                if (new.progress != progress) {
                  progress <- new.progress
                  cat(sprintf("%d%%...", progress))
                }
                else cat(".")
            }
        }
        if (!debug && progress != -1 && progress != 100)
            cat("100%\n")
    }, interrupt = function(interrupt) {
        cat("\n")
        stop("Command interrupted!", call. = FALSE)
    }, finally = {
        if (length(submitted.jobs) > 0) {
            running.jobs <- .gcluster.running.jobs(submitted.jobs)
            answer <- c()
            for (i in 1:length(commands)) {
                res <- list()
                res$exit.status <- NA
                res$retv <- NA
                res$stdout <- NA
                res$stderr <- NA
                if (submitted.jobs[i] %in% running.jobs)
                  res$exit.status <- "interrupted"
                else {
                  fname <- sprintf("%s/%d.retv", tmp.dirname,
                    i)
                  if (file.exists(fname)) {
                    load(fname)
                    res$exit.status <- "success"
                    res$retv <- retv
                  }
                  else res$exit.status <- "failure"
                }
                out.file <- sprintf("%s/%d.out", tmp.dirname,
                  i)
                if (file.exists(out.file)) {
                  f <- file(out.file, "rc")
                  res$stdout <- readChar(f, 1000000)
                  close(f)
                }
                err.file <- sprintf("%s/%d.err", tmp.dirname,
                  i)
                if (file.exists(err.file)) {
                  f <- file(err.file, "rc")
                  res$stderr <- readChar(f, 1000000)
                  close(f)
                }
                answer[[i]] <- res
            }
            for (job in running.jobs) system(sprintf("qdel %s",
                job), ignore.stderr = T, intern = T)
            unlink(tmp.dirname, recursive = TRUE)

            if (collapse_results){
                canswer <- tryCatch(
                    purrr::map_df(answer, function(x) x$retv),
                    error = function(e) {
                        message('returning original output due to an error. collapse your reults manually (are all the parts data frames?)'
                        )
                        return(NULL)
                    })

                if (!is.null(canswer)) {
                    return(canswer)
                }
            }

            return(answer)
        }
        unlink(tmp.dirname, recursive = TRUE)
    })
}
