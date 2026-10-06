doc <- function() {
  c(
    "---",                                   # 1
    "title: test",                           # 2
    "---",                                   # 3
    "",                                      # 4
    "Add one to things.",                    # 5
    "",                                      # 6
    "```{r add-one, echo = TRUE}",           # 7
    "add_one <- function(x) x + 1",          # 8
    "add_one(1:3)",                          # 9
    "```",                                   # 10
    "",                                      # 11
    "```{claude}",                           # 12
    "#| mode: review",                       # 13
    "Is this ok for NA?",                    # 14
    "```",                                   # 15
    "",                                      # 16
    "The end."                               # 17
  )
}

fake_chat <- function(reply = "#### Fine\n\nLooks good.\n\n```r\nf <- function(x) x\n```") {
  calls <- new.env()
  calls$n <- 0
  fn <- function(system, user) {
    calls$n <- calls$n + 1
    calls$system <- system
    calls$user <- user
    reply
  }
  list(fn = fn, calls = calls)
}
