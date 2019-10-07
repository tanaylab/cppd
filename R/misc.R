do_call_ellipsis <- function(f, additional_params=list(), ...){
    f_args <- names(as.list(args(f)))
    elipsis <- list(...)        
    if (!is.null(names(elipsis))){
        new_elipsis <- list()        
        for (x in names(elipsis)[names(elipsis) %in% f_args]) { 
            new_elipsis[[x]] <- elipsis[[x]] 
        }                      
        do.call(f, c(additional_params, new_elipsis))
    } else {
        do.call(f, additional_params)
    }
}


#
#' Set parallel threads
#'
#' @param thread_num number of threads. use '1' for non parallel behaviour
#'
#' @return None
#'
#' @examples
#' cppd.set_parallel(8)
#'
#' @export
cppd.set_parallel <- function(thread_num) {
    if (1 == thread_num) {
        options(cppd.parallel = FALSE)
    } else {
        doMC::registerDoMC(thread_num)
        options(cppd.parallel = TRUE)
        options(cppd.parallel.thread_num = thread_num)
    }
}

memfree <- function() {
    round(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)))
}

