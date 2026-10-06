# RStudio add-ins. Bind them to keyboard shortcuts via
# Tools > Modify Keyboard Shortcuts (search for "Claude").

editor_context <- function(id = NULL) {
  if (!rstudioapi::isAvailable()) stop("This add-in needs RStudio.", call. = FALSE)
  rstudioapi::getSourceEditorContext(id)
}

cursor_row <- function(ctx) ctx$selection[[1]]$range$start[["row"]]

# Apply a line-based edit (see apply_edit()) to an open document.
editor_apply <- function(id, lines, edit) {
  text <- paste(edit$text, collapse = "\n")
  pos <- rstudioapi::document_position
  if (edit$end < edit$start) {
    if (edit$start > length(lines)) {
      rstudioapi::insertText(pos(length(lines), 100000L), paste0("\n", text), id)
    } else {
      rstudioapi::insertText(pos(edit$start, 1L), paste0(text, "\n"), id)
    }
  } else {
    rng <- rstudioapi::document_range(pos(edit$start, 1L), pos(edit$end, 100000L))
    rstudioapi::modifyRange(rng, text, id)
  }
}

#' Ask Claude about the chunk at the cursor
#'
#' With the cursor in a {claude} chunk, sends your question - along with
#' the notebook so far and the output of the R chunk above - and writes the
#' answer below the chunk. Asking again replaces the earlier answer; to
#' follow up instead, add a new {claude} chunk below it.
#'
#' With the cursor in an R chunk, adds an empty {claude} chunk below it and
#' asks about the chunk using the default mode.
#' @export
addin_ask <- function() {
  ctx <- editor_context()
  id <- ctx$id
  lines <- ctx$contents
  row <- cursor_row(ctx)

  chunk <- chunk_at(find_chunks(lines), row)
  if (!is.null(chunk) && chunk$engine == "r") {
    edit <- new_claude_chunk_edit(lines, chunk$end)
    editor_apply(id, lines, edit)
    lines <- apply_edit(lines, edit)
    row <- edit$row
  }

  q <- prepare_question(lines, row)
  message(sprintf("Asking Claude (%s mode)...", q$mode))
  text <- claude_ask(q$system, q$user)

  lines <- editor_context(id)$contents
  chunk_end <- relocate_chunk_end(lines, q$chunk_text, q$chunk$start)
  edit <- response_edit(lines, chunk_end, text, q$mode)
  editor_apply(id, lines, edit)
  rstudioapi::setCursorPosition(rstudioapi::document_position(edit$start + 1L, 1L), id)
  message("Done.")
  invisible()
}

#' Insert a {claude} chunk
#'
#' On a blank line, the chunk goes there; otherwise it goes below the
#' current line, or below the current chunk and its answer.
#' @export
addin_insert <- function() {
  ctx <- editor_context()
  lines <- ctx$contents
  row <- cursor_row(ctx)
  chunk_lines <- c("```{claude}", sprintf("#| mode: %s", default_mode()), "", "```")

  chunk <- chunk_at(find_chunks(lines), row)
  if (!is.null(chunk)) {
    edit <- new_claude_chunk_edit(lines, chunk$end)
    edit$text <- c("", chunk_lines)
  } else if (row <= length(lines) && !nzchar(trimws(lines[row]))) {
    edit <- list(start = row, end = row, text = chunk_lines)
  } else {
    edit <- list(start = row + 1L, end = row, text = chunk_lines)
  }
  editor_apply(ctx$id, lines, edit)
  question_row <- edit$start + which(edit$text == "")[sum(edit$text == "")] - 1L
  rstudioapi::setCursorPosition(rstudioapi::document_position(question_row, 1L), ctx$id)
  invisible()
}

#' Copy Claude's suggested code into a new R chunk
#'
#' Takes the main code block from the answer at (or nearest above) the
#' cursor and puts it in a new R chunk directly below the answer, so
#' adopting a suggestion is a deliberate step you can see in the notebook.
#' @export
addin_adopt <- function() {
  ctx <- editor_context()
  edit <- adopt_edit(ctx$contents, cursor_row(ctx))
  editor_apply(ctx$id, ctx$contents, edit)
  invisible()
}
