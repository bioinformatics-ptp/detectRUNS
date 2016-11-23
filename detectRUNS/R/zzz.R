
.onAttach <- function(libname, pkgname) {
  version = packageVersion("detectRUNS")
  packageStartupMessage(paste("Using detectRUNS", version))
}

.onUnload <- function (libpath) {
  library.dynam.unload("detectRUNS", libpath)
}
