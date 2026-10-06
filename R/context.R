# Building the message Claude sees: the notebook so far, the output of the
# R chunk being asked about, and the question itself.

default_mode <- function() getOption("rlearn.default_mode", "try")

# Work out everything needed to answer the {claude} chunk containing `row`.
prepare_question <- function(lines, row, run = TRUE) {
  chunks <- find_chunks(lines)
  chunk <- chunk_at(chunks, row)
  if (is.null(chunk) || chunk$engine != "claude") {
    stop("Put your cursor inside a {claude} chunk (or an R chunk to have it reviewed).",
         call. = FALSE)
  }
  opts <- chunk_options(lines, chunk)
  mode <- opts$mode %||% default_mode()
  question <- chunk_body(lines, chunk)

  target <- chunks[chunks$end < chunk$start & chunks$engine == "r", , drop = FALSE]
  target <- if (nrow(target)) target[nrow(target), ] else NULL
  output <- NULL
  if (!is.null(target) && !isFALSE(opts$run) && run) {
    output <- run_chunk(chunk_body(lines, target))
  }

  list(
    chunk = chunk,
    chunk_text = lines[chunk$start:chunk$end],
    mode = mode,
    system = system_prompt(mode),
    user = build_user_message(
      notebook = lines[seq_len(chunk$start - 1)],
      question = question,
      mode = mode,
      target_name = if (!is.null(target)) chunk_name(target, lines),
      output = output
    )
  )
}

build_user_message <- function(notebook, question, mode, target_name = NULL, output = NULL) {
  notebook <- paste(notebook, collapse = "\n")
  limit <- getOption("rlearn.max_context_chars", 200000)
  note <- ""
  if (nchar(notebook) > limit) {
    notebook <- substr(notebook, nchar(notebook) - limit + 1, nchar(notebook))
    note <- sprintf(" (only the last %d characters; the start was cut to fit rlearn.max_context_chars)", limit)
  }

  parts <- c(
    paste0("Here is my notebook up to my question", note,
           ". Earlier answers from you appear between ", response_start,
           " and ", response_end, " markers."),
    "<notebook>", notebook, "</notebook>"
  )
  if (!is.null(output)) {
    parts <- c(parts, "",
      sprintf("Output from re-running the R chunk directly above my question (%s) just now:", target_name),
      "<output>", if (length(output)) output else "(no printed output)", "</output>")
  }
  parts <- c(parts, "", sprintf("My question (mode: %s):", mode), "<question>",
             if (length(question)) question else
               "(No question written: respond according to the mode, about the R chunk directly above.)",
             "</question>")
  paste(parts, collapse = "\n")
}

# Re-run a chunk's code in a scratch environment whose parent is the global
# environment, so it sees your objects but can't overwrite them.
run_chunk <- function(code) {
  if (!length(code)) return(character())
  env <- new.env(parent = globalenv())
  res <- tryCatch(
    evaluate::evaluate(paste(code, collapse = "\n"), envir = env, stop_on_error = 1L),
    error = function(e) list(e)
  )
  out <- character()
  for (x in res) {
    if (inherits(x, "source")) next
    if (is.character(x)) {
      out <- c(out, x)
    } else if (inherits(x, "condition")) {
      kind <- if (inherits(x, "error")) "Error" else if (inherits(x, "warning")) "Warning" else NULL
      call <- conditionCall(x)
      # A top-level stop() reports evaluate's own eval() call; hide it.
      if (is.call(call) && identical(call[[1]], quote(eval))) call <- NULL
      msg <- conditionMessage(x)
      out <- c(out, if (is.null(kind)) msg else if (is.null(call)) paste0(kind, ": ", msg)
               else paste0(kind, " in ", paste(deparse(call), collapse = " "), ": ", msg))
    } else if (inherits(x, "recordedplot")) {
      out <- c(out, "[a plot was drawn]")
    }
  }
  out <- trim_blank(unlist(strsplit(out, "\n", fixed = TRUE)))
  max_lines <- getOption("rlearn.max_output_lines", 200)
  if (length(out) > max_lines) {
    out <- c(out[seq_len(max_lines)], sprintf("[... %d more lines not shown]", length(out) - max_lines))
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
