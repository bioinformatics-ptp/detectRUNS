
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Using detectRUNS v0.2.3")
}

.onUnload <- function (libpath) {
  library.dynam.unload("detectRUNS", libpath)
}
