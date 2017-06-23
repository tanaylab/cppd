.onLoad <- function(libname, pkgname) {        
    logging::basicConfig()
    cppd.set_parallel(parallel::detectCores() / 2)
}


