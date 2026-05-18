# chunk_checker.R
# Static validator for R code chunks in this snRNAseq pipeline.
#
# USAGE:
#   source("code/chunk_checker.R")
#   # Select chunk text in the RStudio editor, then run in console:
#   check_chunk()
#   # Or pass code directly:
#   check_chunk(my_code_string)
#   # Hide OK confirmations:
#   check_chunk(show_ok = FALSE)
#
# CHECKS:
#   1. R syntax validity
#   2. Undefined global variables (via codetools)
#   3. Function argument names validated against formals()
#   4. File/directory path existence (evaluates path expressions)
#   5. Checkpoint .rds file existence (load_checkpoint / check_checkpoint calls)
#   6. pkg::fn() package availability

check_chunk <- function(code = NULL, env = .GlobalEnv, show_ok = TRUE) {

  # ── 0. Acquire code ──────────────────────────────────────────────────────
  if (is.null(code)) {
    if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable())
      stop("Supply `code` as a string, or call from inside RStudio.")
    sel <- rstudioapi::getActiveDocumentContext()$selection[[1]]$text
    if (!nzchar(trimws(sel))) stop("No text is selected in the editor.")
    code <- sel
  }

  issues <- list()
  .add <- function(lvl, msg)
    issues[[length(issues) + 1L]] <<- list(lvl = lvl, msg = msg)

  # ── 1. Syntax ─────────────────────────────────────────────────────────────
  parsed <- tryCatch(parse(text = code),
    error = function(e) { .add("ERROR", paste0("Syntax: ", conditionMessage(e))); NULL })
  if (is.null(parsed)) return(invisible(.cc_report(issues)))

  # ── 2. Undefined global variables (codetools) ─────────────────────────────
  fn_wrap <- tryCatch(
    eval(parse(text = paste0("function(){\n", code, "\n}"))),
    error = function(e) NULL
  )

  if (!is.null(fn_wrap) && requireNamespace("codetools", quietly = TRUE)) {
    g <- codetools::findGlobals(fn_wrap, merge = FALSE)

    skip_vars <- c("T", "F", "pi", ".Machine", "LETTERS", "letters",
                   "month.name", "month.abb", ".GlobalEnv", ".BaseNamespaceEnv")
    for (v in setdiff(g$variables, skip_vars)) {
      if (!exists(v, envir = env, inherits = TRUE)) {
        .add("ERROR", sprintf("Undefined variable: `%s`", v))
      } else if (show_ok) {
        val <- get(v, envir = env, inherits = TRUE)
        .add("OK", sprintf("`%s` found [%s]", v, paste(class(val), collapse = "/")))
      }
    }

    # Functions not found in search path
    built_ins <- c(
      "{", "(", "[", "[[", "<-", "->", "=", "<<-",
      "+", "-", "*", "/", "^", "%%", "%/%", "!", "&", "|", "&&", "||",
      "==", "!=", "<", ">", "<=", ">=", "$", "@", ":", "~", "%in%",
      "for", "if", "while", "repeat", "break", "next", "return", "function", "switch",
      "::", ":::", "c", "list", "vector", "character", "numeric", "integer", "logical",
      "data.frame", "matrix", "array", "factor",
      "paste", "paste0", "sprintf", "cat", "message", "warning", "stop", "print",
      "invisible", "tryCatch", "withCallingHandlers", "on.exit",
      "rm", "gc", "Reduce", "Filter", "Map",
      "lapply", "sapply", "vapply", "mapply", "tapply", "do.call",
      "match.arg", "missing", "Sys.time", "Sys.sleep", "proc.time",
      "format", "formatC", "nchar", "trimws", "toupper", "tolower",
      "grepl", "grep", "sub", "gsub", "strsplit", "regmatches", "regexpr",
      "file.exists", "dir.exists", "file.path", "dir.create", "file.rename", "file.copy",
      "readLines", "writeLines", "readRDS", "saveRDS",
      "seq", "seq_len", "seq_along", "length", "nrow", "ncol", "dim",
      "names", "colnames", "rownames", "head", "tail",
      "which", "sort", "order", "rev", "unique", "duplicated", "table", "cumsum",
      "sum", "prod", "mean", "median", "var", "sd", "min", "max", "range",
      "abs", "round", "floor", "ceiling",
      "is.null", "is.na", "is.character", "is.numeric", "is.integer",
      "is.logical", "is.list", "is.data.frame", "is.function",
      "as.character", "as.numeric", "as.integer", "as.logical",
      "as.data.frame", "as.matrix", "as.list",
      "Inf", "NaN", "NA", "NULL", "TRUE", "FALSE",
      "inherits", "class", "attr", "attributes",
      "environment", "parent.env", "globalenv", "baseenv",
      "match", "setdiff", "intersect", "union",
      "identical", "all.equal", "isTRUE", "isFALSE", "any", "all",
      "sqrt", "exp", "log", "log2", "log10",
      "tic", "toc"
    )
    for (f in setdiff(g$functions, built_ins)) {
      found <- tryCatch(
        is.function(get(f, envir = env, inherits = TRUE)),
        error = function(e) FALSE
      )
      if (!found)
        .add("WARNING", sprintf("Function `%s()` not found in search path", f))
    }
  } else if (is.null(fn_wrap)) {
    .add("WARNING", "Could not wrap code for symbol analysis (check syntax or complex quoting).")
  } else {
    .add("WARNING", "`codetools` not available — install it for symbol checks.")
  }

  # ── 3. pkg::fn() package availability ────────────────────────────────────
  .cc_walk(parsed, function(node) {
    if (!is.call(node)) return()
    fn <- node[[1]]
    if (!is.call(fn)) return()
    if (!(as.character(fn[[1]]) %in% c("::", ":::"))) return()
    pkg_name <- as.character(fn[[2]])
    if (!requireNamespace(pkg_name, quietly = TRUE))
      .add("ERROR", sprintf("Package `%s` is not installed.", pkg_name))
    else if (show_ok)
      .add("OK", sprintf("Package `%s` is available.", pkg_name))
  })

  # ── 4. Function argument names vs formals() ───────────────────────────────
  .cc_walk(parsed, function(node) {
    if (!is.call(node)) return()
    fn_sym <- node[[1]]

    if (is.call(fn_sym)) {
      if (!(as.character(fn_sym[[1]]) %in% c("::", ":::"))) return()
      fn_obj    <- tryCatch(eval(fn_sym, envir = env), error = function(e) NULL)
      fn_label  <- paste0(fn_sym[[2]], "::", fn_sym[[3]])
    } else if (is.symbol(fn_sym)) {
      fn_obj   <- tryCatch(get(as.character(fn_sym), envir = env, inherits = TRUE),
                           error = function(e) NULL)
      fn_label <- as.character(fn_sym)
    } else return()

    if (!is.function(fn_obj)) return()

    fmls <- names(formals(fn_obj))
    if ("..." %in% fmls) return()   # function accepts ..., can't validate named args

    node_names <- names(as.list(node)[-1])
    if (is.null(node_names)) return()
    named <- node_names[nzchar(node_names)]
    if (length(named) == 0) return()

    bad <- setdiff(named, fmls)
    if (length(bad) > 0)
      .add("ERROR", sprintf(
        "`%s()` — invalid arg(s): %s\n        Valid: %s",
        fn_label,
        paste(bad, collapse = ", "),
        paste(fmls, collapse = ", ")
      ))
    else if (show_ok)
      .add("OK", sprintf("`%s()` args OK: %s", fn_label, paste(named, collapse = ", ")))
  })

  # ── 5. File/directory path existence ──────────────────────────────────────
  path_build_fns <- c("file.path", "paste0", "paste")

  .cc_walk(parsed, function(node) {
    evaluated <- NULL

    if (is.character(node) && length(node) == 1L) {
      # String literal that looks like a path (contains "/" or starts with ./)
      if (grepl("^[./~]|/", node) && !grepl("^https?://|[\\^$*?]", node))
        evaluated <- node

    } else if (is.call(node) && is.symbol(node[[1]]) &&
               as.character(node[[1]]) %in% path_build_fns) {
      evaluated <- tryCatch({
        r <- eval(node, envir = env)
        if (is.character(r) && length(r) == 1L) r else NULL
      }, error = function(e) NULL, warning = function(w) NULL)
    }

    if (is.null(evaluated) || !nzchar(evaluated)) return()
    if (grepl("[*?]", evaluated)) return()     # skip glob/regex patterns
    if (!grepl("[./]", evaluated)) return()    # skip bare names like "clustering"

    parent <- dirname(evaluated)
    if (!dir.exists(parent) && !file.exists(evaluated)) {
      src <- if (is.character(node)) sprintf('"%s"', node) else deparse(node)[1L]
      .add("WARNING", sprintf(
        "Path parent dir not found: \"%s\"\n        Expression: %s", parent, src))
    } else if (show_ok && grepl("/", evaluated)) {
      .add("OK", sprintf("Path OK: \"%s\"", evaluated))
    }
  })

  # ── 6. Checkpoint file existence ──────────────────────────────────────────
  op  <- tryCatch(get("output_path",      envir = env, inherits = TRUE), error = function(e) NULL)
  pfx <- tryCatch(get("CHECKPOINT_PREFIX", envir = env, inherits = TRUE), error = function(e) "")

  if (is.null(op)) {
    .add("WARNING", "Cannot check checkpoints: `output_path` not found in environment.")
  } else {
    .cc_walk(parsed, function(node) {
      if (!is.call(node)) return()
      fn_sym <- node[[1]]
      if (!is.symbol(fn_sym)) return()
      if (!(as.character(fn_sym) %in% c("load_checkpoint", "check_checkpoint"))) return()

      args     <- as.list(node)[-1]
      name_arg <- if (!is.null(names(args)) && "name" %in% names(args))
        args[["name"]] else args[[1L]]

      if (!is.character(name_arg) || length(name_arg) != 1L) return()

      cp_path <- file.path(op, "RData", paste0(pfx, name_arg, ".rds"))
      if (file.exists(cp_path)) {
        if (show_ok)
          .add("OK", sprintf('Checkpoint "%s" exists.', name_arg))
      } else {
        .add("WARNING", sprintf(
          'Checkpoint "%s" NOT found.\n        Expected: %s', name_arg, cp_path))
      }
    })
  }

  invisible(.cc_report(issues))
}


# ── Internal helpers (prefixed .cc_ to avoid name collisions) ─────────────────

# Depth-first AST walker: calls fn(node) on every node.
.cc_walk <- function(x, fn) {
  if (is.null(x) || is.environment(x)) return()
  fn(x)
  if (is.recursive(x)) lapply(as.list(x), .cc_walk, fn = fn)
  invisible()
}

# Format and print the issue report; return the issues list invisibly.
.cc_report <- function(issues) {
  errors   <- Filter(function(i) i$lvl == "ERROR",   issues)
  warnings <- Filter(function(i) i$lvl == "WARNING", issues)
  oks      <- Filter(function(i) i$lvl == "OK",      issues)

  bar <- strrep("─", 62)
  cat("\n", bar, "\n CHUNK CHECKER\n", bar, "\n", sep = "")

  for (e in errors)   cat(" [ERROR]  ", e$msg, "\n", sep = "")
  for (w in warnings) cat(" [WARN]   ", w$msg, "\n", sep = "")
  for (o in oks)      cat(" [OK]     ", o$msg, "\n", sep = "")

  if (length(issues) == 0L) cat(" ✓ No issues found.\n")

  cat(bar, "\n", sep = "")
  cat(sprintf(" %d error(s)  |  %d warning(s)  |  %d OK\n",
              length(errors), length(warnings), length(oks)))
  cat(bar, "\n\n", sep = "")

  invisible(issues)
}
