# The knitr engine for {claude} chunks. Knitting never calls Claude: the
# answers are already in the document as Markdown. The engine just shows
# your question, styled to match the answer below it.

rlearn_css <- "<style>
.rlearn-question, .rlearn-response { border-left: 4px solid #d97757; padding: .5em 1em; margin: 0 0 1em; border-radius: 4px; }
.rlearn-question { background: rgba(217, 119, 87, .06); }
.rlearn-response { background: rgba(217, 119, 87, .12); }
.rlearn-response h4 { font-size: 1em; font-weight: 600; margin-top: 1em; }
</style>"

rlearn_hide_css <- "<style>.rlearn-question, .rlearn-response { display: none; }</style>"

#' knitr engine for {claude} chunks
#'
#' Registered automatically when the package loads, so put `library(rlearn)`
#' in your setup chunk. Set `options(rlearn.show = FALSE)` before knitting
#' to hide all questions and answers in HTML output.
#'
#' @param options Chunk options, supplied by knitr.
#' @export
eng_claude <- function(options) {
  html <- knitr::is_html_output()
  if (!isTRUE(getOption("rlearn.show", TRUE))) {
    return(if (html) rlearn_hide_css else "")
  }
  if (isFALSE(options$echo)) return("")

  mode <- options$mode %||% default_mode()
  question <- trim_blank(options$code)
  out <- c(
    "",
    if (html) rlearn_css,
    "",
    "::: {.rlearn-question}",
    sprintf("**You asked Claude** &middot; *%s*", mode),
    if (length(question)) c("", question),
    ":::",
    ""
  )
  paste(out, collapse = "\n")
}
