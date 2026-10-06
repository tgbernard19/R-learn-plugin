.onLoad <- function(libname, pkgname) {
  knitr::knit_engines$set(claude = eng_claude)
}
