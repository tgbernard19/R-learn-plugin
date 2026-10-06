# Line-based helpers for reading and editing an R Markdown document.
# Everything here works on a character vector of lines so it can be used
# both from the RStudio add-ins and on a file on disk.

response_start <- "<!-- rlearn:response -->"
response_end <- "<!-- /rlearn:response -->"

chunk_start_re <- "^\\s*(`{3,})\\s*\\{\\s*([A-Za-z0-9_.]+)([^}]*)\\}\\s*$"

# Find every executable chunk. Returns a data frame with one row per chunk:
# start/end line numbers (the fence lines), engine, and label.
find_chunks <- function(lines) {
  starts <- integer(); ends <- integer(); engines <- character(); labels <- character()
  i <- 1L
  n <- length(lines)
  while (i <= n) {
    m <- regmatches(lines[i], regexec(chunk_start_re, lines[i]))[[1]]
    if (length(m)) {
      fence <- m[2]
      close_re <- paste0("^\\s*", fence, "`*\\s*$")
      j <- i + 1L
      while (j <= n && !grepl(close_re, lines[j])) j <- j + 1L
      starts <- c(starts, i)
      ends <- c(ends, min(j, n))
      engines <- c(engines, tolower(m[3]))
      labels <- c(labels, header_label(m[4]))
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  data.frame(start = starts, end = ends, engine = engines, label = labels,
             stringsAsFactors = FALSE)
}

# The label in a header like `{r my-label, echo = FALSE}` is the first
# comma-separated token that isn't a key = value pair.
header_label <- function(rest) {
  rest <- trimws(sub("^,", "", trimws(rest)))
  if (!nzchar(rest)) return(NA_character_)
  first <- trimws(strsplit(rest, ",", fixed = TRUE)[[1]][1])
  if (grepl("=", first, fixed = TRUE) || !nzchar(first)) NA_character_ else first
}

chunk_at <- function(chunks, row) {
  hit <- which(chunks$start <= row & chunks$end >= row)
  if (length(hit)) chunks[hit[1], ] else NULL
}

# Lines strictly inside the fences.
chunk_inner <- function(lines, chunk) {
  if (chunk$end - chunk$start < 2) return(character())
  lines[(chunk$start + 1):(chunk$end - 1)]
}

# Quarto-style `#| key: value` options at the top of a chunk.
chunk_options <- function(lines, chunk) {
  inner <- chunk_inner(lines, chunk)
  opt <- inner[grepl("^\\s*#\\|", inner)]
  if (!length(opt)) return(list())
  parsed <- tryCatch(
    yaml::yaml.load(paste(sub("^\\s*#\\|\\s?", "", opt), collapse = "\n")),
    error = function(e) NULL
  )
  if (is.list(parsed)) parsed else list()
}

# Chunk contents without the `#|` option lines, trimmed of blank edges.
chunk_body <- function(lines, chunk) {
  inner <- chunk_inner(lines, chunk)
  trim_blank(inner[!grepl("^\\s*#\\|", inner)])
}

trim_blank <- function(x) {
  keep <- which(nzchar(trimws(x)))
  if (!length(keep)) return(character())
  x[min(keep):max(keep)]
}

chunk_name <- function(chunk, lines) {
  lab <- chunk_options(lines, chunk)$label
  if (is.null(lab)) lab <- chunk$label
  if (is.null(lab) || is.na(lab)) sprintf("unnamed chunk at line %d", chunk$start) else lab
}

# If a response block directly follows line `after` (blank lines allowed),
# return the line range it occupies, starting just after `after`.
response_after <- function(lines, after) {
  i <- after + 1L
  while (i <= length(lines) && !nzchar(trimws(lines[i]))) i <- i + 1L
  if (i > length(lines) || trimws(lines[i]) != response_start) return(NULL)
  end <- which(trimws(lines) == response_end)
  end <- end[end > i]
  if (!length(end)) return(NULL)
  c(after + 1L, end[1])
}

# All response blocks in the document, as a two-column matrix (start, end).
find_responses <- function(lines) {
  s <- which(trimws(lines) == response_start)
  e <- which(trimws(lines) == response_end)
  out <- matrix(integer(), ncol = 2, dimnames = list(NULL, c("start", "end")))
  for (si in s) {
    ei <- e[e > si]
    if (length(ei)) out <- rbind(out, c(si, ei[1]))
  }
  out
}

# An edit replaces lines start..end with `text` (a character vector of lines).
# end == start - 1 means a pure insertion before line `start`.
apply_edit <- function(lines, edit) {
  before <- if (edit$start > 1) lines[seq_len(edit$start - 1)] else character()
  after <- if (edit$end < length(lines)) lines[(edit$end + 1):length(lines)] else character()
  c(before, edit$text, after)
}

# Strip anything from Claude's reply that would break the document: chunk
# headers (which would become executable) and our own block markers.
sanitize_response <- function(text) {
  text <- gsub("```\\s*\\{\\s*([A-Za-z0-9_.]+)[^}\\n]*\\}", "```\\1", text, perl = TRUE)
  text <- gsub(response_start, "", text, fixed = TRUE)
  text <- gsub(response_end, "", text, fixed = TRUE)
  trim_blank(strsplit(text, "\n", fixed = TRUE)[[1]])
}

format_response <- function(text, mode, model) {
  c(
    "",
    response_start,
    "",
    "::: {.rlearn-response}",
    sprintf("**Claude** \u00b7 *%s* \u00b7 %s", mode, model),
    "",
    sanitize_response(text),
    ":::",
    "",
    response_end
  )
}

# Fenced code blocks inside a range of lines, as a list of
# list(lang, code) in document order.
code_blocks <- function(lines) {
  out <- list()
  i <- 1L
  while (i <= length(lines)) {
    m <- regmatches(lines[i], regexec("^\\s*(`{3,})\\s*([A-Za-z]*)\\s*$", lines[i]))[[1]]
    if (length(m)) {
      close_re <- paste0("^\\s*", m[2], "`*\\s*$")
      j <- i + 1L
      while (j <= length(lines) && !grepl(close_re, lines[j])) j <- j + 1L
      code <- if (j > i + 1L) lines[(i + 1L):(j - 1L)] else character()
      out[[length(out) + 1L]] <- list(lang = tolower(m[3]), code = code)
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  out
}
