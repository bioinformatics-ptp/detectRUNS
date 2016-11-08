
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Using development version for detectRUNS")
}

.onUnload <- function (libpath) {
  library.dynam.unload("detectRUNS", libpath)
}
