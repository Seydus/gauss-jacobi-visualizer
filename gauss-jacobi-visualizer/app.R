library(shiny)
library(DT)
library(ggplot2)

# ---------- Helper Functions ----------

# Parse a space- or comma-separated string into a numeric vector
parse_vector <- function(text, field = "Input") {
  text <- gsub(",", " ", text)
  tokens <- strsplit(trimws(text), "\\s+")[[1]]
  vals <- suppressWarnings(as.numeric(tokens))
  if (any(is.na(vals)))
    stop(paste(field, "must contain only numeric values separated by spaces or commas."))
  vals
}

# Parse a multi-line string into an n x n numeric matrix
parse_matrix <- function(text, n) {
  lines <- strsplit(trimws(text), "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) != n)
    stop(paste("Matrix A must have exactly", n, "rows. You provided", length(lines), "."))
  mat <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
     row <- parse_vector(lines[i], paste("Row", i, "of A"))
    if (length(row) != n)
      stop(paste("Row", i, "of A must have exactly", n, "numeric values."))
    mat[i, ] <- row
  }
  mat
}

# Check strict diagonal dominance: |a_ii| > sum_{j!=i} |a_ij|
is_diagonally_dominant <- function(A) {
  n <- nrow(A)
  for (i in seq_len(n)) {
    if (abs(A[i, i]) <= sum(abs(A[i, -i]))) return(FALSE)
  }
  TRUE
}

# Gauss-Jacobi iterative solver for Ax = b
# Returns the solution, full iteration history, and human-readable steps.
gauss_jacobi <- function(A, b, x0, tol = 1e-6, max_iter = 100) {
  n <- length(b)
  x <- x0

  # Iteration history table
  history <- data.frame(matrix(NA, nrow = max_iter + 1, ncol = n + 2))
  colnames(history) <- c("Iteration", paste0("x", seq_len(n)), "Error")
  history[1, ] <- c(0, round(x0, 8), NA)

  steps_text <- character(0)
  err <- Inf
  final_k <- 0

  for (k in seq_len(max_iter)) {
    x_new <- numeric(n)
    step_lines <- paste0("Iteration ", k)

    for (i in seq_len(n)) {
      # Sum of off-diagonal terms using the OLD x (defining trait of Jacobi)
      off_terms <- A[i, -i] * x[-i]
      s <- sum(off_terms)
      x_new[i] <- (b[i] - s) / A[i, i]

      # Build a readable expression for the Steps tab
      expr <- paste(sprintf("(%.4f)*(%.6f)", A[i, -i], x[-i]), collapse = " + ")
      step_lines <- c(
        step_lines,
        sprintf("  x%d = ( %.4f - [ %s ] ) / %.4f = %.8f",
                i, b[i], expr, A[i, i], x_new[i])
      )
    }

    err <- max(abs(x_new - x))
    history[k + 1, ] <- c(k, round(x_new, 8), err)
    step_lines <- c(step_lines, sprintf("  Max |x_new - x_old| = %.3e", err))
    steps_text <- c(steps_text, paste(step_lines, collapse = "\n"))

    x <- x_new
    final_k <- k
    if (err < tol) break
  }

  # Trim unused rows
  history <- history[seq_len(final_k + 1), , drop = FALSE]

  list(
    solution   = x,
    history    = history,
    steps      = steps_text,
    converged  = (err < tol),
    iterations = final_k,
    final_err  = err
  )
}


# ---------- User Interface ----------

ui <- fluidPage(
  tags$head(
    HTML("<link href='https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap' rel='stylesheet'>"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('darkMode', function(value) {
        document.body.classList.toggle('dark-mode', value);
      });
      Shiny.addCustomMessageHandler('scrollToTop', function(value) {
        setTimeout(function() {
          window.scrollTo(0, 0);
        }, 100);
      });
    "))
  ),
  tags$style(HTML("
    html, body, input, button, select, textarea, table, .dataTables_wrapper {
      font-family: 'Roboto', Arial, sans-serif;
    }
    * { font-family: inherit; font-weight: 400; }
    .title { text-align: center; }
    .theme-toggle { text-align: center; margin: 0 auto 15px; }
    .theme-toggle .form-group { display: inline-block; margin-bottom: 0; }
    footer { text-align: center; padding: 20px 0; }
    .btn-default {
      color: white; background-color: #4a7fb8; border-color: transparent;
      margin: 0 auto; display: block;
    }
    #download_table { margin-top: 12px; }
    pre#answer, pre#error_out, pre#status {
      background: white; font-size: 14px;
    }
    .tab-content { padding-top: 20px; min-height: 650px; }
    h4 { font-weight: 400; line-height: 1.45; }
    .center { font-weight: bold; text-align: center; margin: 30px auto 0; }
    .description { padding: 10px 100px 20px; text-align: left; }
    .intro-section { max-width: 960px; margin: 0 auto 28px; }
    .intro-heading {
      color: #2c5f93; font-weight: 700; text-align: left;
      border-bottom: 2px solid #d8e6f3; padding-bottom: 8px; margin-top: 30px;
    }
    .intro-lead { font-size: 19px; line-height: 1.55; }
    .intro-note {
      background: #eef5fb; border-left: 5px solid #4a7fb8;
      padding: 14px 18px; margin: 16px auto; font-size: 18px; line-height: 1.5;
    }
    .worked-example {
      background: #fbfcfe; border: 1px solid #d8e6f3;
      padding: 18px 22px; margin: 18px auto; border-radius: 8px;
    }
    .worked-example h4 { margin-top: 8px; }
    .formula-line {
      background: #f7f9fc; border: 1px solid #e1e8f0;
      padding: 12px 14px; margin: 10px 0; border-radius: 6px;
      font-size: 17px; line-height: 1.45;
    }
    .example-dropdown {
      background: #ffffff; border: 2px solid #4a7fb8;
      border-radius: 8px; margin: 18px auto; padding: 0;
      overflow: hidden; box-shadow: 0 2px 8px rgba(44, 95, 147, 0.12);
    }
    .example-dropdown summary {
      cursor: pointer; padding: 14px 52px 14px 18px; font-size: 18px;
      font-weight: 700; color: #1f4f7d; background: #eef5fb;
      list-style: none; position: relative;
    }
    .example-dropdown summary::-webkit-details-marker { display: none; }
    .example-dropdown summary::after {
      content: ''; position: absolute; right: 20px; top: 50%;
      width: 0; height: 0; margin-top: -3px;
      border-left: 7px solid transparent; border-right: 7px solid transparent;
      border-top: 9px solid #1f4f7d;
    }
    .example-dropdown[open] summary::after {
      border-top: 0; border-bottom: 9px solid #1f4f7d;
    }
    .example-dropdown > div { padding: 16px 20px 20px; }
    .description ul { font-size: 18px; line-height: 1.55; margin: 12px auto 22px; max-width: 900px; }
    .description li { margin-bottom: 8px; }
    .intro-table {
      border-collapse: separate; border-spacing: 0; width: 100%;
      margin: 20px auto; overflow: hidden; border-radius: 8px;
      border: 1px solid #d8e6f3; font-size: 17px;
    }
    .intro-table th {
      background: #2c5f93; color: white; font-weight: 700;
      padding: 12px 16px; text-align: center;
    }
    .intro-table td {
      padding: 12px 16px; text-align: left; vertical-align: middle;
      border-top: 1px solid #d8e6f3;
    }
    .intro-table tr:nth-child(odd) td { background: #f7f9fc; }
    .intro-table .numeric-cell {
      font-family: 'Roboto', Arial, sans-serif; white-space: nowrap;
      text-align: center; font-weight: 500;
    }
    .matrix-grid {
      border-collapse: collapse; margin: 0 auto; display: inline-table;
      font-size: 16px; background: white;
    }
    .matrix-grid td {
      border: 1px solid #b9cbe0; padding: 6px 12px; min-width: 42px;
      text-align: center; font-weight: 500;
    }
    .vector-grid td { min-width: 54px; }
    .matrix-equation {
      display: flex; align-items: center; justify-content: center;
      gap: 14px; margin: 18px auto 8px; flex-wrap: wrap;
    }
    .matrix-operator {
      font-size: 24px; font-weight: 700; color: #1f4f7d;
    }
    .equation-pair { display: grid; gap: 18px; margin: 16px 0; }
    .equation-row {
      display: grid; grid-template-columns: 1fr 54px 1fr;
      align-items: stretch; gap: 12px;
    }
    .equation-arrow {
      display: flex; align-items: center; justify-content: center;
      color: #1f4f7d; font-size: 28px; font-weight: 700;
    }
    .equation-box {
      border: 1px solid #d8e6f3; border-radius: 8px; padding: 14px;
      background: #fbfcfe;
    }
    .equation-box strong {
      display: block; color: #1f4f7d; font-weight: 700; margin-bottom: 8px;
    }
    .calculation-stack {
      display: grid; gap: 14px; margin: 14px 0;
    }
    .calculation-step {
      border: 1px solid #d8e6f3; border-radius: 8px;
      background: #fbfcfe; padding: 14px 18px;
    }
    .calculation-step strong {
      display: block; color: #1f4f7d; font-weight: 700; margin-bottom: 8px;
    }
    .solution-grid {
      border-collapse: separate; border-spacing: 0; width: 100%;
      margin: 14px 0; border: 1px solid #d8e6f3; border-radius: 8px;
      overflow: hidden; font-size: 17px;
    }
    .solution-grid th {
      background: #eef5fb; color: #1f4f7d; font-weight: 700;
      padding: 10px 14px; border-bottom: 1px solid #d8e6f3;
    }
    .intro-table tr:first-child > *,
    .solution-grid tr:first-child > *,
    .matrix-grid tr:first-child > * {
      text-align: center;
    }
    .solution-grid td {
      padding: 10px 14px; text-align: center; border-top: 1px solid #d8e6f3;
    }
    .tolerance-highlight {
      display: inline-block; background: #fff4cc; color: #4f3900;
      border: 1px solid #e5c65f; border-radius: 6px; padding: 2px 8px;
      font-weight: 700; white-space: nowrap;
    }
    .calculator-jump {
      display: inline-block; margin-top: 10px; font-weight: 700;
      color: white; background: #2c5f93; border-color: #2c5f93;
    }
    .convergence-check {
      display: grid; grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 12px; margin: 14px 0;
    }
    .convergence-check .calculation-step {
      border-left: 5px solid #4a7fb8;
    }
    .dominance-grid {
      display: grid; grid-template-columns: 1fr auto 1fr auto 1fr;
      gap: 8px; align-items: center; margin-top: 10px;
    }
    .dominance-part {
      background: #eef5fb; border: 1px solid #d8e6f3; border-radius: 6px;
      padding: 8px 10px; text-align: center;
    }
    .dominance-label {
      display: block; color: #1f4f7d; font-size: 13px; font-weight: 700;
      margin-bottom: 4px;
    }
    .dominance-value { display: block; font-size: 18px; font-weight: 700; }
    .dominance-symbol { color: #1f4f7d; font-size: 22px; font-weight: 700; }
    .dominance-pass .dominance-value { color: #1d7a3a; }
    .dominance-table td,
    .dominance-table th {
      text-align: center;
      vertical-align: middle;
    }
    .dominance-table .comparison-cell {
      color: #1f4f7d;
      font-size: 22px;
      font-weight: 700;
    }
    .dominance-table .pass-cell {
      color: #1d7a3a;
      font-weight: 700;
    }
    .calculator-page {
      max-width: 1220px;
      margin: 0 auto;
    }
    .calculator-panel {
      background: #fbfcfe;
      border: 1px solid #d8e6f3;
      border-radius: 8px;
      padding: 18px;
      margin-bottom: 18px;
    }
    .calculator-panel h4 {
      color: #1f4f7d;
      font-weight: 700;
      margin-top: 0;
    }
    .calculator-section-label {
      color: #1f4f7d;
      font-weight: 700;
      margin: 18px 0 8px;
    }
    .calculator-grid-wrap {
      overflow-x: auto;
      padding-bottom: 4px;
      width: 100%;
    }
    .calculator-grid {
      display: grid;
      gap: 8px;
      min-width: max-content;
      margin-bottom: 10px;
      align-items: stretch;
    }
    .calculator-grid-cell,
    .calculator-grid-label {
      min-width: 86px;
    }
    .calculator-grid-label {
      background: #eef5fb;
      color: #1f4f7d;
      font-weight: 700;
      text-align: center;
      border: 1px solid #d8e6f3;
      border-radius: 6px;
      padding: 9px 10px;
      display: flex;
      align-items: center;
      justify-content: center;
      min-height: 44px;
    }
    .calculator-grid-cell {
      border: 1px solid #d8e6f3;
      border-radius: 6px;
      padding: 7px;
      text-align: center;
      background: white;
    }
    .calculator-grid .form-group {
      margin-bottom: 0;
    }
    .calculator-grid input.form-control {
      text-align: center;
      padding-left: 6px;
      padding-right: 6px;
      width: 100%;
      min-width: 72px;
      box-shadow: none;
    }
    .calculator-grid input.form-control:focus {
      border-color: #4a7fb8;
      box-shadow: 0 0 0 2px rgba(74, 127, 184, 0.16);
    }
    .grid-pair {
      display: grid;
      gap: 18px;
    }
    .calculator-actions {
      display: flex;
      gap: 10px;
      align-items: center;
      flex-wrap: wrap;
      margin-top: 14px;
    }
    .calculator-actions .btn {
      margin: 0;
      display: inline-block;
    }
    .result-summary {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 14px;
      margin-bottom: 18px;
    }
    .side-result-summary {
      grid-template-columns: 1fr;
      margin-top: 18px;
    }
    .result-card {
      background: #fbfcfe;
      border: 1px solid #d8e6f3;
      border-radius: 8px;
      padding: 14px;
      min-height: 130px;
    }
    .result-card h4 {
      color: #1f4f7d;
      font-weight: 700;
      margin-top: 0;
      font-size: 1.05em;
    }
    .calculator-result-card {
      background: #ffffff;
      color: #222;
      border: 1px solid #d8e6f3;
      border-radius: 8px;
      padding: 14px;
      margin: 0 0 16px;
      display: grid;
      grid-template-columns: minmax(340px, 1.2fr) minmax(220px, 0.8fr) minmax(280px, 1fr);
      gap: 12px;
      align-items: start;
    }
    .calculator-result-card h4 {
      color: #1f4f7d;
      font-weight: 700;
      margin: 0 0 8px;
      font-size: 1em;
    }
    .calculator-result-card pre {
      font-size: 13px;
      margin-bottom: 8px;
      white-space: pre-wrap;
    }
    .calculator-result-item {
      min-width: 0;
    }
    .result-value-grid {
      border-collapse: separate;
      border-spacing: 0;
      width: 100%;
      border: 1px solid #d8e6f3;
      border-radius: 8px;
      overflow: hidden;
      background: #fbfcfe;
    }
    .result-value-grid th {
      background: #eef5fb;
      color: #1f4f7d;
      font-weight: 700;
      padding: 8px 10px;
      text-align: center;
      border-bottom: 1px solid #d8e6f3;
    }
    .result-value-grid td {
      padding: 10px;
      text-align: center;
      font-weight: 500;
      border-left: 1px solid #d8e6f3;
    }
    .result-value-grid td:first-child {
      border-left: 0;
    }
    .result-number {
      display: block;
      background: #fbfcfe;
      border: 1px solid #d8e6f3;
      border-radius: 8px;
      padding: 12px 14px;
      font-size: 18px;
      font-weight: 500;
      text-align: center;
      overflow-wrap: anywhere;
    }
    @media (max-width: 1000px) {
      .result-summary { grid-template-columns: 1fr; }
      .grid-pair { grid-template-columns: 1fr; }
      .calculator-result-card { grid-template-columns: 1fr; }
    }
    @media (max-width: 800px) {
      .convergence-check { grid-template-columns: 1fr; }
      .dominance-grid { grid-template-columns: 1fr; }
      .dominance-symbol { text-align: center; }
    }
    .well {
      background-color: #2d3644; border: 1px solid #0000003d; color: white;
    }
    .formatted_calculations {
      background-color: transparent; border: none; width: 90%;
      counter-reset: section; margin: 0 auto; max-height: 650px; overflow: auto;
    }
    .calculations_box {
      text-align: center; background-color: white; border: 1px solid #0000003d;
      padding: 50px 20px 20px; border-radius: 10px;
      box-shadow: 2px 2px 5px rgba(0,0,0,0.2); position: relative; margin-bottom: 15px;
    }
    .calculations_box:before {
      position: absolute; left: 0px; top: 0px; width: 100%; height: 40px;
      background-color: #4a7fb8; color: white; border-radius: 9px 9px 0 0;
      font-size: 18px; padding: 8px;
      counter-increment: section; content: 'Step ' counter(section);
    }
    .calculations_box pre {
      font-family: 'Courier New', monospace; text-align: left;
      background: #f7f9fc; border: 1px solid #e0e6ed; padding: 12px;
      border-radius: 6px; font-size: 14px; white-space: pre-wrap;
    }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:focus,
    .nav-tabs > li.active > a:hover {
      color: white; cursor: default; background-color: #4a7fb8;
      border: 1px solid #ddd; border-bottom-color: transparent;
    }
    .status-ok { color: #1d7a3a; font-weight: bold; }
    .status-warn { color: #b8860b; font-weight: bold; }
    .status-bad { color: #b22222; font-weight: bold; }
    body.dark-mode {
      background-color: #161b22; color: #e6edf3;
    }
    body.dark-mode .well,
    body.dark-mode .tab-content,
    body.dark-mode .calculations_box,
    body.dark-mode pre,
    body.dark-mode table,
    body.dark-mode .dataTables_wrapper {
      background-color: #0d1117; color: #e6edf3; border-color: #30363d;
    }
    body.dark-mode .form-control,
    body.dark-mode .selectize-input {
      background-color: #161b22; color: #e6edf3; border-color: #30363d;
    }
    body.dark-mode .nav-tabs > li > a {
      color: #e6edf3;
    }
    body.dark-mode .intro-heading { color: #79c0ff; border-bottom-color: #30363d; }
    body.dark-mode .intro-note { background: #0d1117; border-left-color: #79c0ff; }
    body.dark-mode .worked-example,
    body.dark-mode .example-dropdown,
    body.dark-mode .formula-line,
    body.dark-mode .intro-table tr:nth-child(odd) td {
      background: #0d1117; border-color: #30363d;
    }
    body.dark-mode .example-dropdown { border-color: #79c0ff; }
    body.dark-mode .example-dropdown summary { background: #161b22; color: #79c0ff; }
    body.dark-mode .example-dropdown summary::after { border-top-color: #79c0ff; }
    body.dark-mode .example-dropdown[open] summary::after { border-bottom-color: #79c0ff; }
    body.dark-mode .intro-table { border-color: #30363d; }
    body.dark-mode .intro-table th { background: #1f6feb; }
    body.dark-mode .intro-table td { border-top-color: #30363d; }
    body.dark-mode .matrix-grid,
    body.dark-mode .solution-grid,
    body.dark-mode .equation-box,
    body.dark-mode .calculation-step { background: #0d1117; border-color: #30363d; }
    body.dark-mode .matrix-grid td,
    body.dark-mode .solution-grid td,
    body.dark-mode .solution-grid th { border-color: #30363d; }
    body.dark-mode .solution-grid th { background: #161b22; color: #79c0ff; }
    body.dark-mode .matrix-operator,
    body.dark-mode .equation-box strong,
    body.dark-mode .calculation-step strong,
    body.dark-mode .equation-arrow { color: #79c0ff; }
    body.dark-mode .tolerance-highlight {
      background: #3d2f00; color: #f2cc60; border-color: #8f741f;
    }
    body.dark-mode .calculator-jump { background: #1f6feb; border-color: #1f6feb; }
    body.dark-mode .dominance-part {
      background: #161b22; border-color: #30363d;
    }
    body.dark-mode .dominance-label,
    body.dark-mode .dominance-symbol { color: #79c0ff; }
    body.dark-mode .dominance-pass .dominance-value { color: #7ee787; }
    body.dark-mode .dominance-table .comparison-cell { color: #79c0ff; }
    body.dark-mode .dominance-table .pass-cell { color: #7ee787; }
    body.dark-mode .calculator-panel,
    body.dark-mode .result-card {
      background: #0d1117;
      border-color: #30363d;
    }
    body.dark-mode .calculator-result-card {
      background: #161b22;
      color: #e6edf3;
      border-color: #30363d;
    }
    body.dark-mode .calculator-panel h4,
    body.dark-mode .calculator-section-label,
    body.dark-mode .calculator-grid th,
    body.dark-mode .result-card h4,
    body.dark-mode .calculator-result-card h4 {
      color: #79c0ff;
    }
    body.dark-mode .result-value-grid,
    body.dark-mode .result-number {
      background: #0d1117;
      border-color: #30363d;
    }
    body.dark-mode .result-value-grid th {
      background: #161b22;
      color: #79c0ff;
      border-color: #30363d;
    }
    body.dark-mode .result-value-grid td {
      border-color: #30363d;
    }
    body.dark-mode .calculator-grid-label {
      background: #161b22;
      border-color: #30363d;
    }
    body.dark-mode .calculator-grid-cell {
      background: #0d1117;
      border-color: #30363d;
    }
    body.dark-mode .status-ok { color: #7ee787; }
    body.dark-mode .status-warn { color: #f2cc60; }
    body.dark-mode .status-bad { color: #ff7b72; }
    body.dark-mode pre#answer,
    body.dark-mode pre#error_out {
      background-color: #161b22;
      color: #e6edf3;
      border-color: #30363d;
    }
  ")),

  titlePanel(div("Gauss-Jacobi Iterative Method", class = "title")),
  div(class = "theme-toggle", checkboxInput("dark_mode", "Dark mode", value = FALSE)),

  tabsetPanel(id = "main_tabs",

    # ----- Introduction Tab -----
    tabPanel("Introduction",
      fluidRow(
        div(class = "description",

          div(class = "intro-section",
            h4(class = "intro-heading", "1. The Problem"),
            h4(class = "intro-lead", HTML("Many science, engineering, and mathematics problems lead to
                     a <strong>system of linear equations</strong>. The goal is to find the unknown
                     values, called the <strong>solution vector</strong>, that make all equations true
                     at the same time. The <strong>Gauss-Jacobi Method</strong>, also called the
                     <strong>Jacobi iterative method</strong>, solves this kind of problem by improving
                     an initial approximation step by step.")),
            div(class = "intro-note",
              HTML("The next sections explain why the method is useful, what formula it uses,
                    and how to tell whether the iterations are likely to settle toward a solution."))
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "2. Where It Is Used"),
            h4(HTML("Gauss-Jacobi is useful when a problem creates many linear equations and direct
                     hand-solving would be too slow or too crowded to manage.")),
            h4(HTML("<strong>Engineering systems.</strong>
                     Structural and circuit problems often produce large linear systems.
                     <br><br><strong>Computational physics.</strong>
                     Discretized Laplace and Poisson equations often become systems of linear equations.
                     <br><br><strong>Parallel computing.</strong>
                     Since Jacobi uses only values from the previous iteration, many variable updates
                     can be computed independently."))
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "3. Example Problem"),
            h4(HTML("Solve the following system of equations using the Gauss-Jacobi Method.
                     We will use it throughout the introduction so each new idea builds on the
                     previous one.")),
            withMathJax(h4(class = "center",
              "$$\\begin{aligned}
              10x_1 - x_2 + 2x_3 &= 6 \\\\
              -x_1 + 11x_2 - x_3 &= 25 \\\\
              2x_1 - x_2 + 10x_3 &= -11
              \\end{aligned}$$"
            )),
            h4(HTML("Here, <em>x<sub>1</sub></em>, <em>x<sub>2</sub></em>, and
                     <em>x<sub>3</sub></em> are the unknowns. Together, they form the
                     <strong>solution vector</strong>. The numbers multiplying the unknowns are
                     <strong>coefficients</strong>, and the numbers on the right side form the
                     <strong>right-hand side vector</strong>."))
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "4. Matrix Form"),
            h4(HTML("A linear system is usually written as <em>Ax = b</em>. This is a compact
                     way to separate the coefficient matrix, the solution vector, and the
                     right-hand side vector.")),
            HTML("<table class='intro-table'>
                    <tr><th>Part</th><th>Meaning</th><th>From the sample system</th></tr>
                    <tr><td>A</td><td>Coefficient matrix</td>
                        <td class='numeric-cell'>
                          <table class='matrix-grid'>
                            <tr><td>10</td><td>-1</td><td>2</td></tr>
                            <tr><td>-1</td><td>11</td><td>-1</td></tr>
                            <tr><td>2</td><td>-1</td><td>10</td></tr>
                          </table>
                        </td></tr>
                    <tr><td>x</td><td>Solution vector, containing the unknowns</td>
                        <td class='numeric-cell'>
                          <table class='matrix-grid vector-grid'>
                            <tr><td>x<sub>1</sub></td></tr>
                            <tr><td>x<sub>2</sub></td></tr>
                            <tr><td>x<sub>3</sub></td></tr>
                          </table>
                        </td></tr>
                    <tr><td>b</td><td>Right-hand side vector</td>
                        <td class='numeric-cell'>
                          <table class='matrix-grid vector-grid'>
                            <tr><td>6</td></tr>
                            <tr><td>25</td></tr>
                            <tr><td>-11</td></tr>
                          </table>
                        </td></tr>
                 </table>")
            ,
            h4("Combined matrix equation:"),
            HTML("<div class='matrix-equation'>
                    <table class='matrix-grid'>
                      <tr><td>10</td><td>-1</td><td>2</td></tr>
                      <tr><td>-1</td><td>11</td><td>-1</td></tr>
                      <tr><td>2</td><td>-1</td><td>10</td></tr>
                    </table>
                    <span class='matrix-operator'>x</span>
                    <table class='matrix-grid vector-grid'>
                      <tr><td>x<sub>1</sub></td></tr>
                      <tr><td>x<sub>2</sub></td></tr>
                      <tr><td>x<sub>3</sub></td></tr>
                    </table>
                    <span class='matrix-operator'>=</span>
                    <table class='matrix-grid vector-grid'>
                      <tr><td>6</td></tr>
                      <tr><td>25</td></tr>
                      <tr><td>-11</td></tr>
                    </table>
                  </div>")
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "5. Convergence"),
            h4(HTML("Now that the coefficient matrix is visible, we can check whether Jacobi is
                     expected to converge. Convergence means the repeated approximations settle
                     toward a stable solution. A common beginner-friendly test is
                     <strong>strict diagonal dominance</strong>: in every row, the diagonal coefficient
                     should be larger than the combined size of the other coefficients in that row.")),
            withMathJax(h4(class = "center",
              "$$|a_{ii}| > \\sum_{j \\neq i}|a_{ij}|$$"
            )),
            h4(HTML("The values below come directly from the coefficient matrix in section 4.")),
            HTML("<table class='intro-table dominance-table'>
                    <tr>
                      <th>Row</th>
                      <th>Diagonal coefficient</th>
                      <th>Other coefficients in the row</th>
                      <th>Comparison</th>
                      <th>Result</th>
                    </tr>
                    <tr>
                      <td>Row 1</td>
                      <td>|10| = 10</td>
                      <td>|-1| + |2| = 3</td>
                      <td class='comparison-cell'>10 &gt; 3</td>
                      <td class='pass-cell'>Pass</td>
                    </tr>
                    <tr>
                      <td>Row 2</td>
                      <td>|11| = 11</td>
                      <td>|-1| + |-1| = 2</td>
                      <td class='comparison-cell'>11 &gt; 2</td>
                      <td class='pass-cell'>Pass</td>
                    </tr>
                    <tr>
                      <td>Row 3</td>
                      <td>|10| = 10</td>
                      <td>|2| + |-1| = 3</td>
                      <td class='comparison-cell'>10 &gt; 3</td>
                      <td class='pass-cell'>Pass</td>
                    </tr>
                  </table>"),
            h4(HTML("Since every row passes, the sample system is strictly diagonally dominant.
                     This condition guarantees convergence for Jacobi. The worked iterations later
                     show this numerically because the error estimate decreases until it is below
                     the chosen tolerance."))
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "6. General Formula"),
            h4(HTML("The method repeatedly computes a new approximation from the previous one. In
                     common notation, the value of variable <em>x<sub>i</sub></em> at the next
                     iteration is written as:")),
            withMathJax(h4(class = "center",
              "$$x_i^{(k+1)} = \\frac{1}{a_{ii}} \\left( b_i - \\sum_{j \\neq i} a_{ij} \\, x_j^{(k)} \\right)$$"
            )),
            h4(HTML("The superscript <em>k</em> means the previous iteration, and
                     <em>k + 1</em> means the next iteration. The diagonal coefficient
                     <em>a<sub>ii</sub></em> is the coefficient of the variable being solved for.")),
            h4("This app estimates the error using the maximum absolute difference between successive iterates:"),
            withMathJax(h4(class = "center",
              "$$\\varepsilon^{(k+1)} = \\max_{1 \\le i \\le n} \\left| x_i^{(k+1)} - x_i^{(k)} \\right|$$"
            ))
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "7. Method"),
            h4(HTML("Now that the system is organized, Jacobi's idea is to rewrite each equation
                     so one unknown is alone on the left side. That gives us iteration formulas for
                     <em>x<sub>1</sub></em>, <em>x<sub>2</sub></em>, and <em>x<sub>3</sub></em>.")),
            withMathJax(HTML("<div class='equation-pair'>
                <div class='equation-row'>
                  <div class='equation-box'>
                    <strong>Original equation 1</strong>
                    \\[10x_1 - x_2 + 2x_3 = 6\\]
                  </div>
                  <div class='equation-arrow'>&rarr;</div>
                  <div class='equation-box'>
                    <strong>Rearranged iteration formula for x<sub>1</sub></strong>
                    \\[x_1 = \\dfrac{6 + x_2 - 2x_3}{10}\\]
                  </div>
                </div>
                <div class='equation-row'>
                  <div class='equation-box'>
                    <strong>Original equation 2</strong>
                    \\[-x_1 + 11x_2 - x_3 = 25\\]
                  </div>
                  <div class='equation-arrow'>&rarr;</div>
                  <div class='equation-box'>
                    <strong>Rearranged iteration formula for x<sub>2</sub></strong>
                    \\[x_2 = \\dfrac{25 + x_1 + x_3}{11}\\]
                  </div>
                </div>
                <div class='equation-row'>
                  <div class='equation-box'>
                    <strong>Original equation 3</strong>
                    \\[2x_1 - x_2 + 10x_3 = -11\\]
                  </div>
                  <div class='equation-arrow'>&rarr;</div>
                  <div class='equation-box'>
                    <strong>Rearranged iteration formula for x<sub>3</sub></strong>
                    \\[x_3 = \\dfrac{-11 - 2x_1 + x_2}{10}\\]
                  </div>
                </div>
              </div>")),
            tags$details(class = "example-dropdown",
              tags$summary("Show full rearrangement for each equation"),
              div(
                h4(HTML("Equation 1: solve for <em>x<sub>1</sub></em>.")),
                withMathJax(h4(class = "center",
                  "$$\\begin{array}{l}
                  10x_1 - x_2 + 2x_3 = 6 \\\\[6pt]
                  10x_1 = 6 + x_2 - 2x_3 \\\\[6pt]
                  \\dfrac{10x_1}{10} = \\dfrac{6 + x_2 - 2x_3}{10} \\\\[6pt]
                  x_1 = \\dfrac{6 + x_2 - 2x_3}{10}
                  \\end{array}$$"
                )),
                h4(HTML("Equation 2: solve for <em>x<sub>2</sub></em>.")),
                withMathJax(h4(class = "center",
                  "$$\\begin{array}{l}
                  -x_1 + 11x_2 - x_3 = 25 \\\\[6pt]
                  11x_2 = 25 + x_1 + x_3 \\\\[6pt]
                  \\dfrac{11x_2}{11} = \\dfrac{25 + x_1 + x_3}{11} \\\\[6pt]
                  x_2 = \\dfrac{25 + x_1 + x_3}{11}
                  \\end{array}$$"
                )),
                h4(HTML("Equation 3: solve for <em>x<sub>3</sub></em>.")),
                withMathJax(h4(class = "center",
                  "$$\\begin{array}{l}
                  2x_1 - x_2 + 10x_3 = -11 \\\\[6pt]
                  10x_3 = -11 - 2x_1 + x_2 \\\\[6pt]
                  \\dfrac{10x_3}{10} = \\dfrac{-11 - 2x_1 + x_2}{10} \\\\[6pt]
                  x_3 = \\dfrac{-11 - 2x_1 + x_2}{10}
                  \\end{array}$$"
                ))
              )
            ),
            h4(HTML("Unlike Gauss-Seidel, Jacobi updates every unknown using only the values from
                     the <strong>previous iterate</strong>. That previous iterate is used to compute
                     the <strong>next iterate</strong>.")),
            HTML("<ul>
                    <li>Choose an initial approximation, such as <strong>0 0 0</strong>.</li>
                    <li>Use the iteration formulas to compute the next iterate.</li>
                    <li>Compare the next iterate with the previous iterate.</li>
                    <li>Repeat until the error estimate is smaller than the chosen tolerance.</li>
                 </ul>"),
            div(class = "intro-note",
              HTML("After the original equations are rearranged, those rearranged formulas become
                    the formulas used in every Jacobi iteration. The next section shows how the
                    first new values are computed from them."))
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "8. Iteration"),
            HTML("<ul>
                    <li>This section demonstrates the second method step: <strong>use the iteration formulas to compute the next iterate</strong>.</li>
                    <li>Start with the initial approximation <em>x<sup>(0)</sup> = (0, 0, 0)</em>.</li>
                    <li>The entries of <em>x<sup>(0)</sup></em> are used during the first Jacobi iteration.</li>
                  </ul>"),
            tags$details(class = "example-dropdown",
              tags$summary("Show detailed iteration 1 and 2 solution"),
              div(
                h4("Reference: rearranged formulas used for the computations."),
                withMathJax(HTML("<div class='calculation-stack'>
                    <div class='calculation-step'>
                      \\[x_1 = \\dfrac{6 + x_2 - 2x_3}{10}\\]
                    </div>
                    <div class='calculation-step'>
                      \\[x_2 = \\dfrac{25 + x_1 + x_3}{11}\\]
                    </div>
                    <div class='calculation-step'>
                      \\[x_3 = \\dfrac{-11 - 2x_1 + x_2}{10}\\]
                    </div>
                  </div>")),
                h4("Step 1: Use the previous iterate in each iteration formula."),
                withMathJax(HTML("<div class='calculation-stack'>
                    <div class='calculation-step'>
                      <strong>Compute \\(x_1\\) for iteration 1</strong>
                      \\[x_1^{(1)} = \\dfrac{6 + 0 - 2(0)}{10}\\]
                      \\[x_1^{(1)} = \\dfrac{6}{10}\\]
                      \\[x_1^{(1)} = 0.6000\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute \\(x_2\\) for iteration 1</strong>
                      \\[x_2^{(1)} = \\dfrac{25 + 0 + 0}{11}\\]
                      \\[x_2^{(1)} = \\dfrac{25}{11}\\]
                      \\[x_2^{(1)} = 2.2727\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute \\(x_3\\) for iteration 1</strong>
                      \\[x_3^{(1)} = \\dfrac{-11 - 2(0) + 0}{10}\\]
                      \\[x_3^{(1)} = \\dfrac{-11}{10}\\]
                      \\[x_3^{(1)} = -1.1000\\]
                    </div>
                  </div>")),
                h4("Step 2: Place the computed values into the next iterate."),
                HTML("<table class='solution-grid'>
                        <tr><th>x<sub>1</sub><sup>(1)</sup></th><th>x<sub>2</sub><sup>(1)</sup></th><th>x<sub>3</sub><sup>(1)</sup></th></tr>
                        <tr><td>0.6000</td><td>2.2727</td><td>-1.1000</td></tr>
                      </table>"),
                h4("Step 3: Compute the error estimate."),
                withMathJax(HTML("<div class='calculation-step'>
                  \\[\\varepsilon^{(1)} =
                  \\max\\left(
                  |0.6000 - 0|,
                  |2.2727 - 0|,
                  |-1.1000 - 0|
                  \\right)\\]
                  \\[\\varepsilon^{(1)} =
                  \\max(0.6000, 2.2727, 1.1000)\\]
                  \\[\\varepsilon^{(1)} = 2.2727\\]
                </div>")),
                h4(HTML("This maximum absolute difference is the error estimate for iteration 1.
                         Since it is still larger than the tolerance, the process does not stop yet.
                         The newly computed vector becomes the previous iterate for iteration 2.")),
                h4("Step 4: Use iteration 1 as the previous iterate to compute iteration 2."),
                withMathJax(HTML("<div class='calculation-stack'>
                    <div class='calculation-step'>
                      <strong>Compute \\(x_1\\) for iteration 2</strong>
                      \\[x_1^{(2)} = \\dfrac{6 + x_2^{(1)} - 2x_3^{(1)}}{10}\\]
                      \\[x_1^{(2)} = \\dfrac{6 + 2.2727 - 2(-1.1000)}{10}\\]
                      \\[x_1^{(2)} = \\dfrac{10.4727}{10}\\]
                      \\[x_1^{(2)} = 1.0473\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute \\(x_2\\) for iteration 2</strong>
                      \\[x_2^{(2)} = \\dfrac{25 + x_1^{(1)} + x_3^{(1)}}{11}\\]
                      \\[x_2^{(2)} = \\dfrac{25 + 0.6000 + (-1.1000)}{11}\\]
                      \\[x_2^{(2)} = \\dfrac{24.5000}{11}\\]
                      \\[x_2^{(2)} = 2.2273\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute \\(x_3\\) for iteration 2</strong>
                      \\[x_3^{(2)} = \\dfrac{-11 - 2x_1^{(1)} + x_2^{(1)}}{10}\\]
                      \\[x_3^{(2)} = \\dfrac{-11 - 2(0.6000) + 2.2727}{10}\\]
                      \\[x_3^{(2)} = \\dfrac{-9.9273}{10}\\]
                      \\[x_3^{(2)} = -0.9927\\]
                    </div>
                  </div>")),
                h4("Step 5: Compute the iteration 2 error estimate."),
                withMathJax(HTML("<div class='calculation-step'>
                  \\[\\varepsilon^{(2)} =
                  \\max\\left(
                  |1.0473 - 0.6000|,
                  |2.2273 - 2.2727|,
                  |-0.9927 - (-1.1000)|
                  \\right)\\]
                  \\[\\varepsilon^{(2)} =
                  \\max(0.4473, 0.0454, 0.1073)\\]
                  \\[\\varepsilon^{(2)} = 0.4473\\]
                </div>")),
                h4(HTML("After iteration 2, the same cycle repeats: use the newest values in the
                         rearranged formulas, compute the next vector, then compute a new error
                         estimate. The full iteration table below shows these repeated updates until
                         the error estimate becomes smaller than the highlighted tolerance."))
              )
            ),
            h4(HTML("Full iteration table using tolerance <span class='tolerance-highlight'>1e-6</span>:")),
            HTML("<table class='solution-grid'>
                    <tr><th>Iteration</th><th>x<sub>1</sub></th><th>x<sub>2</sub></th><th>x<sub>3</sub></th><th>Error estimate</th></tr>
                    <tr><td>0</td><td>0.00000000</td><td>0.00000000</td><td>0.00000000</td><td>-</td></tr>
                    <tr><td>1</td><td>0.60000000</td><td>2.27272727</td><td>-1.10000000</td><td>2.27272727</td></tr>
                    <tr><td>2</td><td>1.04727273</td><td>2.22727273</td><td>-0.99272727</td><td>0.44727273</td></tr>
                    <tr><td>3</td><td>1.02127273</td><td>2.27768595</td><td>-1.08672727</td><td>0.09400000</td></tr>
                    <tr><td>4</td><td>1.04511405</td><td>2.26677686</td><td>-1.07648595</td><td>0.02384132</td></tr>
                    <tr><td>5</td><td>1.04197488</td><td>2.26987528</td><td>-1.08234512</td><td>0.00585917</td></tr>
                    <tr><td>6</td><td>1.04345655</td><td>2.26905725</td><td>-1.08140745</td><td>0.00148168</td></tr>
                    <tr><td>7</td><td>1.04318721</td><td>2.26927719</td><td>-1.08178559</td><td>0.00037814</td></tr>
                    <tr><td>8</td><td>1.04328484</td><td>2.26921833</td><td>-1.08170972</td><td>0.00009762</td></tr>
                    <tr><td>9</td><td>1.04326378</td><td>2.26923410</td><td>-1.08173513</td><td>0.00002541</td></tr>
                    <tr><td>10</td><td>1.04327044</td><td>2.26922988</td><td>-1.08172935</td><td>0.00000666</td></tr>
                    <tr><td>11</td><td>1.04326886</td><td>2.26923101</td><td>-1.08173110</td><td>0.00000175</td></tr>
                    <tr><td>12</td><td>1.04326932</td><td>2.26923071</td><td>-1.08173067</td><td>0.00000046</td></tr>
                  </table>"),
            h4("Stopping result:"),
            HTML("<table class='solution-grid'>
                    <tr><th>Stopping iteration</th><th>Tolerance</th><th>Final error estimate</th></tr>
                    <tr><td>12</td><td><span class='tolerance-highlight'>1e-6</span></td><td>0.00000046</td></tr>
                  </table>"),
            h4("Approximate solution:"),
            HTML("<table class='solution-grid'>
                    <tr><th>x<sub>1</sub></th><th>x<sub>2</sub></th><th>x<sub>3</sub></th></tr>
                    <tr><td>1.04326932</td><td>2.26923071</td><td>-1.08173067</td></tr>
                  </table>"),
            div(class = "intro-note",
              HTML("Conclusion: the error estimate decreases from <strong>2.27272727</strong>
                    at iteration 1 to <strong>0.00000046</strong> at iteration 12. Since
                    <strong>0.00000046</strong> is smaller than the tolerance
                    <span class='tolerance-highlight'>1e-6</span>, the method stops and the
                    iteration 12 values are used as the approximate solution."))
          ),

          div(class = "intro-section",
            h4(class = "intro-heading", "9. Try It In The Calculator"),
            h4(HTML("Now that the system, initial approximation, iteration formula, and stopping criterion are clear,
                     these are the exact inputs needed by the calculator.")),
            HTML("<table class='intro-table'>
                    <tr><th>Field</th><th>Format</th><th>Sample input</th></tr>
                    <tr><td>Number of unknowns</td><td>Positive integer</td><td class='numeric-cell'>3</td></tr>
                    <tr><td>Coefficient matrix A</td><td>One row per line, values separated by spaces or commas</td>
                        <td class='numeric-cell'>
                          <table class='matrix-grid'>
                            <tr><td>10</td><td>-1</td><td>2</td></tr>
                            <tr><td>-1</td><td>11</td><td>-1</td></tr>
                            <tr><td>2</td><td>-1</td><td>10</td></tr>
                          </table>
                        </td></tr>
                    <tr><td>Right-hand side vector b</td><td>Single line, values separated by spaces or commas</td>
                        <td class='numeric-cell'>
                          <table class='matrix-grid vector-grid'>
                            <tr><td>6</td></tr>
                            <tr><td>25</td></tr>
                            <tr><td>-11</td></tr>
                          </table>
                        </td></tr>
                    <tr><td>Initial guess x<sup>(0)</sup></td><td>Single line, one starting value for each unknown</td>
                        <td class='numeric-cell'>
                          <table class='matrix-grid vector-grid'>
                            <tr><td>0</td></tr>
                            <tr><td>0</td></tr>
                            <tr><td>0</td></tr>
                          </table>
                        </td></tr>
                    <tr><td>Tolerance</td><td>Positive real number</td><td class='numeric-cell'><span class='tolerance-highlight'>1e-6</span></td></tr>
                    <tr><td>Max iterations</td><td>Positive integer</td><td class='numeric-cell'>50</td></tr>
                 </table>"),
            actionButton("go_to_calculator", "Go to Calculator", class = "calculator-jump")
          )
        ),
        align = "center"
      )
    ),

    # ----- Calculate Tab -----
    tabPanel("Calculate",
      fluidRow(
        column(width = 3,
          sidebarPanel(width = 12,
            h4("Calculator"),
            numericInput("n", "Number of unknowns:", value = 3, min = 2, max = 8, step = 1),
            textAreaInput("matA", "Coefficient matrix A:",
                          value = "10 -1 2\n-1 11 -1\n2 -1 10",
                          rows = 4, resize = "vertical"),
            textInput("vecB", "Right-hand side vector b:", value = "6 25 -11"),
            textInput("x0",   HTML("Initial guess x<sup>(0)</sup>:"), value = "0 0 0"),
            numericInput("tol", "Tolerance:", value = 1e-6, min = 1e-15, step = 1e-6),
            numericInput("max_iter", "Max iterations:", value = 50, min = 1, max = 1000, step = 1),
            actionButton("calculate", "Calculate"),
            downloadButton("download_table", "Download iteration table")
          ),
        ),

        column(width = 8,
          div(class = "calculator-result-card",
            div(class = "calculator-result-item",
              h4("Approximate solution vector"),
              uiOutput("answer")
            ),
            div(class = "calculator-result-item",
              h4("Final approximation error"),
              uiOutput("error_out")
            ),
            div(class = "calculator-result-item",
              h4("Convergence status"),
              uiOutput("status")
            )
          ),
          tabsetPanel(
            tabPanel("Convergence Plot",
              plotOutput("plot", height = "500px"),
              h4(style = "font-size: 1.1em; text-align: justify;",
                 HTML(paste(
                   "<div style='max-width: 900px; margin: 30px 50px;'>",
                   "This plot shows the <strong>error at each iteration</strong>",
                   "(the maximum absolute change in any component of x) on a logarithmic scale.",
                   "A steady downward trend indicates the method is converging toward the true solution.",
                   "If the curve grows or oscillates, the system likely fails the diagonal-dominance condition.",
                   "</div>"
                 )))
            ),
            tabPanel("Solution Plot",
              plotOutput("solution_plot", height = "500px"),
              h4(style = "font-size: 1.1em; text-align: justify;",
                 HTML(paste(
                   "<div style='max-width: 900px; margin: 30px 50px;'>",
                   "This plot shows how each variable changes across iterations.",
                   "Lines that flatten out indicate that the variable values are stabilizing.",
                   "</div>"
                 )))
            ),
            tabPanel("Steps", uiOutput("formatted_calculations")),
            tabPanel("Iteration Table",
              dataTableOutput("table"),
              h4(style = "font-size: 1.1em; text-align: justify;",
                 HTML(paste(
                   "<div style='max-width: 900px; margin: 30px 20px;'>",
                   "Each row records one full Gauss-Jacobi iteration. The columns",
                   "<strong>x<sub>1</sub>, x<sub>2</sub>, ...</strong> are the components of the current estimate,",
                   "and <strong>Error</strong> is the maximum absolute change since the previous",
                   "iterate. Iteration stops when the error falls below the chosen tolerance",
                   "or the maximum iteration count is reached.",
                   "</div>"
                 )))
            )
          )
        )
      )
    )
  ),

  tags$footer(
    strong("Numerical Methods Final Activity"),
    br(),
    "Gauss-Jacobi Method - 2026"
  )
)


# ---------- Server Logic ----------

default_A_value <- function(i, j, n) {
  sample_A <- matrix(c(10, -1, 2,
                       -1, 11, -1,
                       2, -1, 10), nrow = 3, byrow = TRUE)
  if (n == 3) return(sample_A[i, j])
  if (i == j) return(1)
  0
}

default_b_value <- function(i, n) {
  sample_b <- c(6, 25, -11)
  if (n == 3) return(sample_b[i])
  0
}

default_x0_value <- function(i, n) {
  0
}

server <- function(input, output, session) {
  observe({
      session$sendCustomMessage("darkMode", isTRUE(input$dark_mode))
    })

  observeEvent(input$go_to_calculator, {
    updateTabsetPanel(session, "main_tabs", selected = "Calculate")
    session$sendCustomMessage("scrollToTop", TRUE)
  })

  output$matrix_grid_input <- renderUI({
    n <- as.integer(input$n)
    headers <- c(
      tags$div(class = "calculator-grid-label", ""),
      lapply(seq_len(n), function(j) {
        tags$div(class = "calculator-grid-label", HTML(paste0("x<sub>", j, "</sub>")))
      })
    )
    rows <- lapply(seq_len(n), function(i) {
      cells <- lapply(seq_len(n), function(j) {
        id <- paste0("a_", i, "_", j)
        value <- isolate(input[[id]])
        if (is.null(value)) value <- default_A_value(i, j, n)
        tags$div(class = "calculator-grid-cell",
          numericInput(id, label = NULL, value = value, step = 1, width = "100%")
        )
      })
      c(tags$div(class = "calculator-grid-label", paste("Row", i)), cells)
    })

    div(class = "calculator-grid-wrap",
      tags$div(
        class = "calculator-grid",
        style = paste0("grid-template-columns: repeat(", n + 1, ", minmax(86px, 1fr));"),
        headers,
        unlist(rows, recursive = FALSE)
      )
    )
  })

  output$b_grid_input <- renderUI({
    n <- as.integer(input$n)
    headers <- lapply(seq_len(n), function(i) {
      tags$div(class = "calculator-grid-label", HTML(paste0("b<sub>", i, "</sub>")))
    })
    cells <- lapply(seq_len(n), function(i) {
      id <- paste0("b_", i)
      value <- isolate(input[[id]])
      if (is.null(value)) value <- default_b_value(i, n)
      tags$div(class = "calculator-grid-cell",
        numericInput(id, label = NULL, value = value, step = 1, width = "100%")
      )
    })

    div(class = "calculator-grid-wrap",
      tags$div(
        class = "calculator-grid",
        style = paste0("grid-template-columns: repeat(", n, ", minmax(86px, 1fr));"),
        headers,
        cells
      )
    )
  })

  output$x0_grid_input <- renderUI({
    n <- as.integer(input$n)
    headers <- lapply(seq_len(n), function(i) {
      tags$div(class = "calculator-grid-label", HTML(paste0("x<sub>", i, "</sub><sup>(0)</sup>")))
    })
    cells <- lapply(seq_len(n), function(i) {
      id <- paste0("x0_", i)
      value <- isolate(input[[id]])
      if (is.null(value)) value <- default_x0_value(i, n)
      tags$div(class = "calculator-grid-cell",
        numericInput(id, label = NULL, value = value, step = 1, width = "100%")
      )
    })

    div(class = "calculator-grid-wrap",
      tags$div(
        class = "calculator-grid",
        style = paste0("grid-template-columns: repeat(", n, ", minmax(86px, 1fr));"),
        headers,
        cells
      )
    )
  })

  read_numeric_cell <- function(id, label) {
    value <- input[[id]]
    validate(
      need(!is.null(value), paste(label, "is not ready yet.")),
      need(!is.na(value), paste(label, "must be numeric."))
    )
    value
  }

  read_matrix_grid <- function(n) {
    matrix(
      unlist(lapply(seq_len(n), function(i) {
        lapply(seq_len(n), function(j) {
          read_numeric_cell(paste0("a_", i, "_", j), paste0("A[", i, ",", j, "]"))
        })
      })),
      nrow = n,
      byrow = TRUE
    )
  }

  read_vector_grid <- function(prefix, n, label) {
    unlist(lapply(seq_len(n), function(i) {
      read_numeric_cell(paste0(prefix, "_", i), paste0(label, "[", i, "]"))
    }))
  }

  # Validate inputs and run Gauss-Jacobi. Returns a list or stops with a clear message.
  result <- reactive({
    req(input$calculate > 0)

    validate(
      need(!is.null(input$n) && input$n >= 2, "n must be at least 2."),
      need(nzchar(input$matA), "Matrix A is empty."),
      need(nzchar(input$vecB), "Vector b is empty."),
      need(nzchar(input$x0),   "Initial guess x0 is empty."),
      need(input$tol > 0, "Tolerance must be positive."),
      need(input$max_iter >= 1, "Max iterations must be at least 1.")
    )

    n <- as.integer(input$n)

    A  <- tryCatch(parse_matrix(input$matA, n),
                   error = function(e) { validate(need(FALSE, e$message)); NULL })
    b  <- tryCatch(parse_vector(input$vecB, "Vector b"),
                   error = function(e) { validate(need(FALSE, e$message)); NULL })
    x0 <- tryCatch(parse_vector(input$x0, "Initial guess x0"),
                   error = function(e) { validate(need(FALSE, e$message)); NULL })
    validate(
      need(length(b)  == n, paste("Vector b must have", n, "values; got", length(b), ".")),
      need(length(x0) == n, paste("Initial guess must have", n, "values; got", length(x0), ".")),
      need(all(diag(A) != 0), "Diagonal entries of A must be non-zero (cannot divide by zero).")
    )

    res <- gauss_jacobi(A, b, x0, tol = input$tol, max_iter = input$max_iter)
    res$dom <- is_diagonally_dominant(A)
    res$A   <- A
    res$b   <- b
    res
  })

  # Final solution vector
  output$answer <- renderUI({
    res <- result()
    labels <- paste0("x<sub>", seq_along(res$solution), "</sub>")
    values <- format(round(res$solution, 8), nsmall = 6, trim = TRUE)

    HTML(paste0(
      "<table class='result-value-grid'><tr>",
      paste0("<th>", labels, "</th>", collapse = ""),
      "</tr><tr>",
      paste0("<td>", values, "</td>", collapse = ""),
      "</tr></table>"
    ))
  })

  # Final error
  output$error_out <- renderUI({
    res <- result()
    HTML(paste0("<span class='result-number'>", formatC(res$final_err, format = "e", digits = 6), "</span>"))
  })

  # Convergence status banner
  output$status <- renderUI({
    res <- result()
    parts <- list()

    if (res$converged) {
      parts <- append(parts, list(div(class = "status-ok",
        sprintf("Converged in %d iterations.", res$iterations))))
    } else {
      parts <- append(parts, list(div(class = "status-bad",
        sprintf("Did NOT converge within %d iterations (final error = %.3e).",
                res$iterations, res$final_err))))
    }

    if (res$dom) {
      parts <- append(parts, list(div(class = "status-ok",
        "Matrix A is strictly diagonally dominant - convergence is guaranteed.")))
    } else {
      parts <- append(parts, list(div(class = "status-warn",
        "Warning: Matrix A is NOT strictly diagonally dominant. Convergence is not guaranteed.")))
    }

    do.call(tagList, parts)
  })

  # Convergence plot (error vs. iteration on log scale)
  output$plot <- renderPlot({
    res <- result()
    df <- res$history[-1, , drop = FALSE]   # drop iteration 0 (NA error)

    ggplot(df, aes(x = Iteration, y = Error)) +
      geom_line(color = "#4a7fb8", size = 1.1) +
      geom_point(color = "#b22222", size = 2.5) +
      scale_y_log10() +
      geom_hline(yintercept = input$tol, linetype = "dashed", color = "#1d7a3a") +
      annotate("text",
               x = max(df$Iteration), y = input$tol,
               label = sprintf("tolerance = %.0e", input$tol),
               vjust = -0.5, hjust = 1, color = "#1d7a3a") +
      labs(title = "Convergence of Gauss-Jacobi Iteration",
           x = "Iteration",
           y = "Max |x_new - x_old|  (log scale)") +
      theme_minimal(base_size = 14)
  })

  # Solution plot (x values vs. iteration)
    output$solution_plot <- renderPlot({
      res <- result()
      x_cols <- grep("^x", colnames(res$history), value = TRUE)

      solution_df <- data.frame(
        Iteration = rep(res$history$Iteration, times = length(x_cols)),
        Variable = rep(x_cols, each = nrow(res$history)),
        Value = unlist(res$history[x_cols], use.names = FALSE)
      )

      ggplot(solution_df, aes(x = Iteration, y = Value, color = Variable)) +
        geom_line(size = 1.1) +
        geom_point(size = 2.2) +
        labs(title = "Solution Values per Iteration",
            x = "Iteration",
            y = "Approximate value") +
        theme_minimal(base_size = 14)
    })

  # Iteration table
  output$table <- renderDataTable({
    res <- result()
    datatable(res$history,
              options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE),
              rownames = FALSE) |>
      formatRound(columns = colnames(res$history)[-1], digits = 8)
  })

  # Download iteration table as CSV
  output$download_table <- downloadHandler(
    filename = function() {
      paste0("gauss_jacobi_iterations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res <- result()
      write.csv(res$history, file, row.names = FALSE)
    }
  )


  # Step-by-step calculations
  output$formatted_calculations <- renderUI({
    res <- result()
    boxes <- lapply(res$steps, function(txt) {
      div(class = "calculations_box", tags$pre(txt))
    })
    fluidRow(div(class = "formatted_calculations", do.call(tagList, boxes)))
  })
}


# ---------- Run the App ----------

shinyApp(ui, server)
