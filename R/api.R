# Calling Claude. There's no official Anthropic SDK for R, so this talks to
# the Messages API directly.

claude_model <- function() getOption("rlearn.model", "claude-opus-5-5")

# Send one question and return Claude's reply as a single string.
# Set options(rlearn.chat_fn = function(system, user) ...) to swap in your own
# backend (or a fake one for testing).
claude_ask <- function(system, user) {
  fn <- getOption("rlearn.chat_fn")
  if (is.function(fn)) return(fn(system = system, user = user))

  key <- Sys.getenv("ANTHROPIC_API_KEY")
  if (!nzchar(key)) {
    stop("ANTHROPIC_API_KEY isn't set. Add it to ~/.Renviron ",
         "(usethis::edit_r_environ() opens the file) and restart R.", call. = FALSE)
  }

  body <- list(
    model = claude_model(),
    max_tokens = getOption("rlearn.max_tokens", 16000L),
    system = system,
    messages = list(list(role = "user", content = user)),
    output_config = list(effort = getOption("rlearn.effort", "medium")),
    # If a safety classifier declines, retry on Anthropic's recommended
    # fallback model instead of returning nothing.
    fallbacks = "default"
  )

  req <- httr2::request("https://api.anthropic.com/v1/messages")
  req <- httr2::req_headers(req,
    `x-api-key` = key,
    `anthropic-version` = "2023-06-01",
    `anthropic-beta` = "server-side-fallback-2026-07-01",
    .redact = "x-api-key"
  )
  req <- httr2::req_body_json(req, body, auto_unbox = TRUE)
  req <- httr2::req_timeout(req, getOption("rlearn.timeout", 600))
  req <- httr2::req_retry(req, max_tries = 3)
  req <- httr2::req_error(req, body = function(resp) {
    tryCatch(httr2::resp_body_json(resp)$error$message, error = function(e) NULL)
  })

  parse_response(httr2::resp_body_json(httr2::req_perform(req)))
}

parse_response <- function(resp) {
  if (identical(resp$stop_reason, "refusal")) {
    stop("Claude declined to answer this one. Try rephrasing the question.", call. = FALSE)
  }
  texts <- vapply(
    Filter(function(b) identical(b$type, "text"), resp$content),
    function(b) b$text, character(1)
  )
  text <- paste(texts, collapse = "\n\n")
  if (identical(resp$stop_reason, "max_tokens")) {
    text <- paste0(text, "\n\n*(Reply cut off at the length limit; raise options(rlearn.max_tokens).)*")
  }
  text
}
