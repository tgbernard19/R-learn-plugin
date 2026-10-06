# rlearn

Talk to Claude about your R code **inside your R Markdown notebook**, one function at a time.

You write a function, or an English description of one, in an R chunk. You ask a question in a `{claude}` chunk under it. Claude's answer is written into the notebook directly below your question. Claude sees your notebook so far and the output of your chunk, so there's no copying and pasting, and the conversation stays next to the code it's about.

The package is built around **try first, then get a thorough explanation**. Claude may write code, but the instructions require it to explain everything it writes: what changed from your version, why, which R concepts are involved, and how to check the result. The aim is to get good at specifying, reading and judging code, not at memorising syntax.

## What it looks like

````markdown
Takes a data frame, a grouping column and a numeric column. Returns one row
per group with the mean and the count of non-missing values. Ignores NAs.

```{r}
summarise_by_group <- function(df, group, value) {
  groups <- unique(df[[group]])
  out <- data.frame()
  for (g in groups) {
    x <- df[df[[group]] == g, value]
    out <- rbind(out, data.frame(group = g, mean = mean(x), n = length(x)))
  }
  out
}
summarise_by_group(airquality, "Month", "Ozone")
```

```{claude}
#| mode: try
Does this handle the NAs in Ozone the way I described?
```
````

Put the cursor in either chunk and press your "Ask Claude" shortcut. Claude gets the notebook, the chunk's re-run output (here, a column of `NA` means) and your question. Its answer is inserted below the chunk:

```markdown
<!-- rlearn:response -->

::: {.rlearn-response}
**Claude** · *try* · claude-opus-5-5

#### What you're aiming for
...
#### What needs to change
...
#### Revised version
...
#### Line by line
...
#### Check it
...
:::

<!-- /rlearn:response -->
```

In the editor the answer is plain Markdown. When you knit, your questions and Claude's answers show up as shaded callouts beneath the code. Knitting never calls Claude; the answers are already in the file.

## Setup

1. Install:

   ```r
   # install.packages("remotes")
   remotes::install_github("tgbernard19/R-learn-plugin")
   ```

2. Add your Anthropic API key to `~/.Renviron` (`usethis::edit_r_environ()` opens the file), then restart R:

   ```
   ANTHROPIC_API_KEY=sk-ant-...
   ```

3. Bind keyboard shortcuts: **Tools → Modify Keyboard Shortcuts**, search for "Claude", and pick any free key combinations.

   | Add-in | What it does |
   |---|---|
   | **Ask Claude** | Cursor in a `{claude}` chunk: answers it. Cursor in an R chunk: adds a `{claude}` chunk below it and asks about it. |
   | **Insert Claude chunk** | Inserts an empty `{claude}` chunk. |
   | **Adopt Claude's code** | Copies the main code block from the answer at the cursor into a new R chunk below the answer. |

4. Start a notebook from **File → New File → R Markdown → From Template → Function notebook (rlearn)**, or add `library(rlearn)` to the setup chunk of any `.Rmd`. Knitting needs that line so knitr knows about `{claude}` chunks.

## Modes

Set the mode with `#| mode:` in a `{claude}` chunk. The default is `try`; change it with `options(rlearn.default_mode = "review")`.

| Mode | When to use it | What you get |
|---|---|---|
| `try` | You've had a go and want to move on | What you were aiming for, what works, what's wrong (with an input that breaks it), a revised version, a line-by-line walkthrough of the changes, and calls to check it |
| `explain` | You don't understand some code, output or an error | A step-by-step trace of what R is doing, with intermediate values |
| `review` | It works and you want to know if it's *good* | Code-review findings ranked must-fix / should-fix / optional, plus testthat expectations to add |
| `hint` | You want to work it out yourself | A nudge only; no solution |
| `ask` | Anything else | A direct answer |

The wording behind each mode is in [`R/prompts.R`](R/prompts.R). That file shapes what you learn more than anything else in the package, so edit it to suit you.

## How it works

- **Follow-ups:** add another `{claude}` chunk below an answer. Claude sees the whole notebook above the question, including its earlier answers, so it reads like a conversation.
- **Re-asking:** running "Ask Claude" on the same chunk again replaces its answer.
- **Output:** before asking, the add-in re-runs the nearest R chunk above the question in a scratch environment (its parent is your global environment, so it sees your objects but can't overwrite them). Claude gets the printed output, warnings and errors; for plots it only learns that a plot was drawn. If that chunk has side effects you don't want repeated (writing files, querying a database, a slow model fit), add `#| run: false` to the `{claude}` chunk.
- **Knitting:** `{claude}` chunks render as "You asked Claude" callouts. Set `options(rlearn.show = FALSE)` before knitting to hide all questions and answers, e.g. for a clean report.
- **Outside RStudio:** `rlearn::ask_claude_file("analysis.Rmd", line = 42)` does what "Ask Claude" does, on a file on disk.

## Options

| Option | Default | |
|---|---|---|
| `rlearn.model` | `"claude-opus-5-5"` | Model to ask |
| `rlearn.effort` | `"medium"` | How hard Claude thinks: `"low"`, `"medium"`, `"high"`, `"xhigh"`, `"max"`. Lower is faster and cheaper. |
| `rlearn.default_mode` | `"try"` | Mode for chunks without `#| mode:` |
| `rlearn.max_tokens` | `16000` | Maximum reply length |
| `rlearn.max_context_chars` | `200000` | Notebook text sent with each question; longer notebooks are cut from the start, and Claude is told |
| `rlearn.max_output_lines` | `200` | Chunk output lines sent |
| `rlearn.show` | `TRUE` | Show questions and answers when knitting |
| `rlearn.chat_fn` | `NULL` | `function(system, user)` returning a string, to use a different backend |

## Things to know

- **Your notebook and chunk output are sent to Anthropic's API** with every question. Check your organisation's policy before using this on work code or data, and use `#| run: false` if the output itself is sensitive.
- Each question is a separate API call billed to your key. With the whole notebook as context, a question usually costs a few cents, more for long notebooks; `options(rlearn.effort = "low")` makes it cheaper and faster.
- R is blocked while Claude answers, usually for a few seconds to a minute. You can keep typing in the editor.
- RStudio's green "run chunk" button doesn't work on `{claude}` chunks; use the add-in.

## Development

```r
devtools::test()   # or testthat::test_local()
```

The add-ins are thin wrappers around line-based functions in `R/document.R` and `R/answer.R`, which are tested without RStudio or network access.
