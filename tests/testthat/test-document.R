test_that("find_chunks finds engines, labels and ranges", {
  ch <- find_chunks(doc())
  expect_equal(ch$start, c(7, 12))
  expect_equal(ch$end, c(10, 15))
  expect_equal(ch$engine, c("r", "claude"))
  expect_equal(ch$label, c("add-one", NA))
})

test_that("plain code fences are not chunks", {
  lines <- c("```r", "x", "```", "```{r}", "y", "```")
  ch <- find_chunks(lines)
  expect_equal(nrow(ch), 1)
  expect_equal(ch$start, 4)
})

test_that("chunk options and body are separated", {
  lines <- doc()
  ch <- find_chunks(lines)[2, ]
  expect_equal(chunk_options(lines, ch), list(mode = "review"))
  expect_equal(chunk_body(lines, ch), "Is this ok for NA?")
})

test_that("responses are inserted, then replaced on re-ask", {
  lines <- doc()
  edit <- response_edit(lines, 15, "First answer.", "review")
  once <- apply_edit(lines, edit)
  expect_equal(once[1:15], lines[1:15])
  expect_true("First answer." %in% once)
  expect_equal(tail(once, 2), c("", "The end."))
  expect_equal(nrow(find_chunks(once)), 2)

  twice <- apply_edit(once, response_edit(once, 15, "Second answer.", "review"))
  expect_false("First answer." %in% twice)
  expect_true("Second answer." %in% twice)
  expect_equal(nrow(find_responses(twice)), 1)
  expect_equal(length(twice), length(once))
})

test_that("chunk headers in replies are defused", {
  out <- sanitize_response("Try:\n```{r, eval=TRUE}\nx\n```\n<!-- /rlearn:response -->")
  expect_equal(out, c("Try:", "```r", "x", "```"))
})

test_that("code blocks in a response don't confuse chunk detection", {
  lines <- apply_edit(doc(), response_edit(doc(), 15,
    "Here:\n```r\nf <- function(x) x\n```\n\n```\nplain\n```", "try"))
  expect_equal(find_chunks(lines)$engine, c("r", "claude"))
})

test_that("relocate_chunk_end follows a chunk that moved", {
  lines <- doc()
  text <- lines[12:15]
  moved <- c("new line", "another", lines)
  expect_equal(relocate_chunk_end(moved, text, 12), 17)
  expect_error(relocate_chunk_end(lines[-14], text, 12), "changed")
})

test_that("adopt_edit takes the function definition, not the check calls", {
  reply <- "#### Revised version\n```r\nf <- function(x) {\n  x + 1\n}\n```\n\n#### Check it\n```r\nf(1)\n```"
  lines <- apply_edit(doc(), response_edit(doc(), 15, reply, "try"))
  resp <- find_responses(lines)
  edit <- adopt_edit(lines, resp[1, "start"] + 3)
  out <- apply_edit(lines, edit)
  ch <- find_chunks(out)
  expect_equal(ch$engine, c("r", "claude", "r"))
  new <- ch[3, ]
  expect_equal(chunk_body(out, new), c("f <- function(x) {", "  x + 1", "}"))
  # Cursor below the answer also works.
  expect_equal(adopt_edit(lines, length(lines))$text, edit$text)
  expect_error(adopt_edit(doc(), 17), "No Claude answer")
})
