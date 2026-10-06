# Putting answers into the document. The same steps are used by the RStudio
# add-ins (which edit the open editor) and by ask_claude_file() (which edits
# a file on disk, for use outside RStudio).

# Edit that writes Claude's reply under the chunk ending at `chunk_end`,
# replacing an earlier reply if there is one.
response_edit <- function(lines, chunk_end, text, mode) {
  new <- format_response(text, mode, claude_model())
  existing <- response_after(lines, chunk_end)
  if (is.null(existing)) {
    list(start = chunk_end + 1L, end = chunk_end, text = new)
  } else {
    list(start = existing[1], end = existing[2], text = new)
  }
}

# Edit that adds an empty {claude} chunk below an R chunk (and below any
# answer already attached to it). `row` is a line inside the new chunk.
new_claude_chunk_edit <- function(lines, after) {
  resp <- response_after(lines, after)
  if (!is.null(resp)) after <- resp[2]
  list(
    start = after + 1L, end = after,
    text = c("", "```{claude}", sprintf("#| mode: %s", default_mode()), "```"),
    row = after + 3L
  )
}

# Claude can take a while to answer and RStudio lets you keep typing, so find
# the chunk again rather than trusting the old line numbers.
relocate_chunk_end <- function(lines, chunk_text, hint) {
  n <- length(chunk_text)
  starts <- which(lines == chunk_text[1])
  starts <- starts[vapply(starts, function(s) {
    s + n - 1L <= length(lines) && identical(lines[s:(s + n - 1L)], chunk_text)
  }, logical(1))]
  if (!length(starts)) {
    stop("The {claude} chunk changed while Claude was answering. Ask again.", call. = FALSE)
  }
  starts[which.min(abs(starts - hint))] + n - 1L
}

# Edit that copies the main code block from a Claude answer into a new R
# chunk below that answer. Uses the answer containing `row`, or the nearest
# one above it.
adopt_edit <- function(lines, row) {
  blocks <- find_responses(lines)
  inside <- which(blocks[, "start"] <= row & blocks[, "end"] >= row)
  above <- which(blocks[, "end"] < row)
  i <- if (length(inside)) inside[1] else if (length(above)) max(above) else NULL
  if (is.null(i)) stop("No Claude answer at or above the cursor.", call. = FALSE)

  b <- blocks[i, ]
  code <- Filter(function(x) x$lang %in% c("r", ""), code_blocks(lines[b["start"]:b["end"]]))
  if (!length(code)) stop("That answer has no R code to adopt.", call. = FALSE)
  # In "try" answers the first function definition is the revised version;
  # later blocks are usually the calls for checking it.
  defines <- vapply(code, function(x) any(grepl("function\\s*\\(", x$code)), logical(1))
  pick <- code[[if (any(defines)) which(defines)[1] else 1L]]

  list(start = b[["end"]] + 1L, end = b[["end"]],
       text = c("", "```{r}", pick$code, "```"))
}

#' Answer a {claude} chunk in a file on disk
#'
#' The same thing the "Ask Claude" add-in does, for use outside RStudio.
#' If `line` is inside an R chunk, a new {claude} chunk is added below it
#' first.
#'
#' @param path Path to an .Rmd or .qmd file.
#' @param line A line number inside the chunk.
#' @export
ask_claude_file <- function(path, line) {
  lines <- readLines(path, warn = FALSE)
  chunk <- chunk_at(find_chunks(lines), line)
  if (!is.null(chunk) && chunk$engine == "r") {
    edit <- new_claude_chunk_edit(lines, chunk$end)
    lines <- apply_edit(lines, edit)
    line <- edit$row
  }
  q <- prepare_question(lines, line)
  text <- claude_ask(q$system, q$user)
  writeLines(apply_edit(lines, response_edit(lines, q$chunk$end, text, q$mode)), path)
  invisible(path)
}
