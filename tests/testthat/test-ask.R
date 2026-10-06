test_that("prepare_question includes notebook, output and question", {
  q <- prepare_question(doc(), 14)
  expect_equal(q$mode, "review")
  expect_match(q$system, "Mode: review")
  expect_match(q$user, "add_one <- function", fixed = TRUE)
  expect_match(q$user, "[1] 2 3 4", fixed = TRUE)
  expect_match(q$user, "(add-one)", fixed = TRUE)
  expect_match(q$user, "Is this ok for NA?", fixed = TRUE)
  # The {claude} chunk itself isn't repeated in the notebook part.
  expect_false(grepl("```{claude}", q$user, fixed = TRUE))
})

test_that("errors and warnings from the chunk are reported", {
  out <- run_chunk(c("x <- log(-1)", "stop('nope')", "print('not reached')"))
  expect_match(out[1], "^Warning in log\\(-1\\)")
  expect_match(out[2], "^Error: nope")
  expect_length(out, 2)
})

test_that("re-running a chunk doesn't change the global environment", {
  run_chunk("rlearn_test_object <- 1")
  expect_false(exists("rlearn_test_object", envir = globalenv()))
})

test_that("run: false skips re-running", {
  lines <- doc()
  lines[13] <- "#| run: false"
  q <- prepare_question(lines, 14)
  expect_false(grepl("<output>", q$user, fixed = TRUE))
  expect_equal(q$mode, "try")
})

test_that("cursor outside a {claude} chunk is an error", {
  expect_error(prepare_question(doc(), 5), "cursor")
  expect_error(prepare_question(doc(), 8), "cursor")
})

test_that("unknown modes are rejected", {
  lines <- doc()
  lines[13] <- "#| mode: wizard"
  expect_error(prepare_question(lines, 14), "Unknown mode")
})

test_that("ask_claude_file answers a {claude} chunk", {
  fake <- fake_chat()
  withr_opts <- options(rlearn.chat_fn = fake$fn)
  on.exit(options(withr_opts))
  path <- tempfile(fileext = ".Rmd")
  writeLines(doc(), path)

  ask_claude_file(path, 13)
  out <- readLines(path)
  expect_equal(fake$calls$n, 1)
  expect_true("Looks good." %in% out)
  expect_equal(nrow(find_responses(out)), 1)

  ask_claude_file(path, 13)
  expect_equal(nrow(find_responses(readLines(path))), 1)
})

test_that("ask_claude_file on an R chunk adds a {claude} chunk below it", {
  fake <- fake_chat()
  old <- options(rlearn.chat_fn = fake$fn)
  on.exit(options(old))
  path <- tempfile(fileext = ".Rmd")
  writeLines(doc(), path)

  ask_claude_file(path, 8)
  out <- readLines(path)
  ch <- find_chunks(out)
  expect_equal(ch$engine, c("r", "claude", "claude"))
  expect_equal(ch$start[2], 12)
  expect_match(fake$calls$system, "Mode: try")
  expect_match(fake$calls$user, "No question written", fixed = TRUE)
  # The answer sits directly under the new chunk.
  expect_equal(response_after(out, ch$end[2])[1], ch$end[2] + 1)
})

test_that("parse_response handles text, refusals and truncation", {
  ok <- list(stop_reason = "end_turn", content = list(
    list(type = "thinking", thinking = ""),
    list(type = "text", text = "a"),
    list(type = "text", text = "b")))
  expect_equal(parse_response(ok), "a\n\nb")
  expect_error(parse_response(list(stop_reason = "refusal", content = list())), "declined")
  cut <- list(stop_reason = "max_tokens", content = list(list(type = "text", text = "a")))
  expect_match(parse_response(cut), "cut off")
})
