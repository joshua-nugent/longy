.onLoad <- function(libname, pkgname) {
  # Ensure data.table knows about our package
  # Suppress the "Column X of item Y is missing" messages from rbindlist
  invisible()
}
