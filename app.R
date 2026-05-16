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
    tags$title("Gauss-Jacobi Iterative Method"),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('scrollToTop', function(value) {
        setTimeout(function() {
          window.scrollTo(0, 0);
        }, 100);
      });
      function setCalculationLoader(visible) {
        var loaders = [
          document.getElementById('calculator_loading_bar'),
          document.getElementById('calculation_results_loader')
        ].filter(Boolean);
        loaders.forEach(function(loader) {
          if (visible) {
            loader.dataset.shownAt = String(Date.now());
            loader.classList.add('is-visible');
            return;
          }
          var shownAt = Number(loader.dataset.shownAt || 0);
          var delay = Math.max(0, 550 - (Date.now() - shownAt));
          setTimeout(function() {
            loader.classList.remove('is-visible');
          }, delay);
        });
      }
      Shiny.addCustomMessageHandler('calculationFinished', function(value) {
        setCalculationLoader(false);
      });
      function setCalculatePageLoader(visible) {
        var loader = document.getElementById('calculate_page_loader');
        if (!loader) return;
        if (visible) {
          loader.dataset.shownAt = String(Date.now());
          loader.classList.add('is-visible');
          return;
        }
        var shownAt = Number(loader.dataset.shownAt || 0);
        var delay = Math.max(0, 450 - (Date.now() - shownAt));
        setTimeout(function() {
          loader.classList.remove('is-visible');
        }, delay);
      }
      Shiny.addCustomMessageHandler('calculatePageLoading', function(value) {
        setCalculatePageLoader(!!value);
      });
      document.addEventListener('click', function(event) {
        var tab = event.target && event.target.closest ? event.target.closest('a[data-value=\"Calculate\"]') : null;
        if (tab) setCalculatePageLoader(true);
      });
      document.addEventListener('pointerdown', function(event) {
        var button = event.target && event.target.closest ? event.target.closest('#calculate') : null;
        if (!button) return;
        setCalculationLoader(true);
      }, true);
      document.addEventListener('click', function(event) {
        var button = event.target && event.target.closest ? event.target.closest('#calculate') : null;
        if (!button) return;
        setCalculationLoader(true);
      }, true);
      document.addEventListener('input', function(event) {
        if (!window.Shiny || !event.target || !event.target.id) return;
        if (/^(a_\\d+_\\d+|b_\\d+|x0_\\d+|n|tol|max_iter)$/.test(event.target.id)) {
          Shiny.setInputValue('calculator_inputs_changed', {
            id: event.target.id,
            value: event.target.value,
            timestamp: Date.now()
          }, {priority: 'event'});
        }
      }, true);
      document.addEventListener('change', function(event) {
        if (!window.Shiny || !event.target || !event.target.id) return;
        if (/^(a_\\d+_\\d+|b_\\d+|x0_\\d+|n|tol|max_iter)$/.test(event.target.id)) {
          Shiny.setInputValue('calculator_inputs_changed', {
            id: event.target.id,
            value: event.target.value,
            timestamp: Date.now()
          }, {priority: 'event'});
        }
      }, true);
    "))
  ),
  tags$style(HTML("
    /* Self-contained UI framework styles. Keep this in app.R so the app can be submitted as one R file. */
    .tw-calculator-page{margin-left:auto;margin-right:auto;width:100%;max-width:1500px;background-color:rgb(248 250 252);padding:.5rem 1.5rem 2.5rem}
    .tw-calculator-hero{margin-bottom:1.25rem;display:flex;align-items:flex-end;justify-content:space-between;gap:1rem;border-radius:1rem;border:1px solid rgb(226 232 240);background-color:#fff;padding:1.25rem 1.5rem;box-shadow:0 1px 2px rgba(0,0,0,.05)}
    .tw-calculator-title{margin:0;font-size:1.875rem;line-height:2.25rem;font-weight:700;letter-spacing:0;color:rgb(2 6 23)}
    .tw-calculator-subtitle{margin:.5rem 0 0;max-width:56rem;font-size:1rem;line-height:1.75rem;color:rgb(71 85 105)}
    .tw-calculator-layout{display:grid;grid-template-columns:repeat(1,minmax(0,1fr));gap:1.25rem}
    @media (min-width:1280px){.tw-calculator-layout{grid-template-columns:minmax(560px,.9fr) minmax(520px,1.1fr);align-items:flex-start}}
    .tw-calculator-card{border-radius:1rem;border:1px solid rgb(226 232 240);background-color:#fff;padding:1.5rem;box-shadow:0 1px 2px rgba(0,0,0,.05)}
    .tw-calculator-input-card{min-width:0}
    .tw-calculator-card-title{margin:0;font-size:1.25rem;line-height:1.75rem;font-weight:700;color:rgb(2 6 23)}
    .tw-calculator-help{margin-top:.5rem;font-size:.875rem;line-height:1.5rem;color:rgb(100 116 139)}
    .tw-story-flow{margin-top:1.25rem;display:grid;gap:1.25rem}
    .tw-story-step{overflow:hidden;border-radius:1rem;border:1px solid rgb(226 232 240);background-color:rgb(248 250 252);padding:1rem}
    .tw-story-step-header{margin-bottom:.75rem;display:flex;align-items:flex-start;gap:.75rem}
    .tw-step-number{display:flex;height:2.25rem;width:2.25rem;flex-shrink:0;align-items:center;justify-content:center;border-radius:9999px;background-color:#284F78;font-size:.875rem;line-height:1.25rem;font-weight:700;color:#fff}
    .tw-story-step-title{margin:0;font-size:1.125rem;line-height:1.75rem;font-weight:700;color:rgb(2 6 23)}
    .tw-story-step-copy{margin:.25rem 0 0;font-size:.875rem;line-height:1.5rem;color:rgb(100 116 139)}
    .tw-settings-grid{display:grid;grid-template-columns:repeat(1,minmax(0,1fr));gap:.75rem}
    @media (min-width:1024px){.tw-settings-grid{grid-template-columns:repeat(3,minmax(0,1fr))}}
    .tw-equation-inputs{position:relative;display:flex;width:100%;min-width:0;max-width:100%;align-items:flex-end;gap:.75rem;overflow-x:auto;overflow-y:hidden;padding-bottom:.75rem}
    .tw-equation-sign{display:flex;flex-shrink:0;align-items:center;justify-content:center;align-self:stretch;padding-top:2.75rem;font-size:1.5rem;line-height:2rem;font-weight:700;color:rgb(100 116 139)}
    .tw-equation-sign-spacer,.tw-equation-spacer-label{display:none}
    .tw-equation-sign-symbol{display:flex;flex:1 1 0%;align-items:center;justify-content:center;padding-left:.25rem;padding-right:.25rem}
    .tw-equation-part,.tw-equation-vector-part{position:relative;display:flex;width:max-content;flex-shrink:0;flex-direction:column;align-items:center;padding-top:2.75rem}
    .tw-calculator-section-label{position:absolute;left:50%;top:1rem;margin:0;display:flex;height:1.25rem;transform:translateX(-50%);align-items:center;justify-content:center;gap:.5rem;white-space:nowrap;text-align:center;font-size:.75rem;line-height:1rem;font-weight:700;text-transform:uppercase;letter-spacing:.025em;color:#356291}
    .tw-calculator-grid{display:grid;width:100%;align-items:stretch;gap:.5rem}
    .tw-vector-table-wrap{max-width:none;overflow:visible;border-radius:.75rem;padding:.5rem;width:max-content;background-color:rgb(248 250 252)}
    .tw-vector-table{width:max-content;background-color:rgb(248 250 252);table-layout:fixed;border-collapse:collapse}
    .tw-vector-table td,.tw-vector-table th{box-sizing:border-box;border:2px solid rgb(203 213 225);padding:0;text-align:center;vertical-align:middle}
    .tw-vector-table th{background-color:rgb(241 245 249);padding:.75rem 1rem;font-size:1rem;line-height:1.5rem;font-weight:700;color:rgb(71 85 105)}
    .tw-vector-table .form-group{margin-bottom:0}
    .tw-vector-table .control-label{display:none}
    .tw-unknown-vector{width:max-content;border-collapse:collapse;background-color:rgb(248 250 252)}
    .tw-unknown-vector td{box-sizing:border-box;height:64px;min-width:60px;border:2px solid rgb(203 213 225);background-color:rgb(241 245 249);padding:0;text-align:center;vertical-align:middle;font-size:1.5rem;line-height:2rem;font-weight:700;color:rgb(100 116 139)}
    .tw-calculator-actions{display:flex;flex-wrap:wrap;align-items:center;gap:.75rem}
    .tw-solve-panel{display:grid;gap:1rem;border-radius:1rem;border:1px solid #D9E8F6;background-color:#EEF6FC;padding:1rem}
    @media (min-width:1024px){.tw-solve-panel{grid-template-columns:minmax(0,1fr) auto;align-items:center}}
    .tw-solve-title{margin:0;font-size:1.125rem;line-height:1.75rem;font-weight:700;color:rgb(2 6 23)}
    .tw-solve-copy{margin:.25rem 0 0;font-size:.875rem;line-height:1.5rem;color:rgb(71 85 105)}
    .tw-results-heading{margin:0 0 1rem;font-size:1.25rem;line-height:1.75rem;font-weight:700;color:rgb(2 6 23)}
    .tw-result-grid{margin-bottom:1.25rem;display:grid;grid-template-columns:repeat(1,minmax(0,1fr));gap:1rem}
    @media (min-width:1024px){.tw-result-grid{grid-template-columns:repeat(3,minmax(0,1fr))}}
    .tw-result-card{min-width:0;overflow:hidden;border-radius:1rem;border:1px solid rgb(226 232 240);background-color:#fff;padding:1.25rem;box-shadow:0 1px 2px rgba(0,0,0,.05)}
    .tw-result-card-title{margin:0 0 .5rem;font-size:.875rem;line-height:1.25rem;font-weight:700;text-transform:uppercase;letter-spacing:.025em;color:rgb(100 116 139)}
    .tw-calculator-tabs{border-radius:1rem;border:1px solid rgb(226 232 240);background-color:#fff;padding:1.25rem;box-shadow:0 1px 2px rgba(0,0,0,.05)}
    .tw-chart-note{margin-top:1rem;width:100%;max-width:none;border-radius:.75rem;background-color:rgb(248 250 252);padding:1rem;font-size:1rem;line-height:1.75rem;color:rgb(71 85 105)}
    html, body, input, button, select, textarea, table, .dataTables_wrapper {
      font-family: Arial, Helvetica, sans-serif;
    }
    * { font-family: inherit; font-weight: 400; }
    .title { text-align: center; font-weight: 700; }
    .title-credit {
      color: #475569;
      font-size: 16px;
      font-weight: 400;
      line-height: 1.6;
      margin-top: 8px;
    }
    .title-credit-label {
      color: #284F78;
      font-weight: 700;
      margin-right: 6px;
    }
    footer { text-align: center; padding: 20px 0; }
    .btn-default {
      color: white; background-color: #427AB5; border-color: transparent;
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
      color: #356291; font-weight: 700; text-align: left;
      border-bottom: 2px solid #D9E8F6; padding-bottom: 8px; margin-top: 30px;
    }
    .intro-lead { font-size: 19px; line-height: 1.55; }
    .intro-note {
      background: #EEF6FC; border-left: 5px solid #427AB5;
      padding: 14px 18px; margin: 16px auto; font-size: 18px; line-height: 1.5;
    }
    .worked-example {
      background: #F8FBFE; border: 1px solid #D9E8F6;
      padding: 18px 22px; margin: 18px auto; border-radius: 8px;
    }
    .worked-example h4 { margin-top: 8px; }
    .formula-line {
      background: #f7f9fc; border: 1px solid #e1e8f0;
      padding: 12px 14px; margin: 10px 0; border-radius: 6px;
      font-size: 17px; line-height: 1.45;
    }
    .example-dropdown {
      background: #ffffff; border: 2px solid #427AB5;
      border-radius: 8px; margin: 18px auto; padding: 0;
      overflow: hidden; box-shadow: 0 2px 8px rgba(66, 122, 181, 0.12);
    }
    .example-dropdown summary {
      cursor: pointer; padding: 14px 52px 14px 18px; font-size: 18px;
      font-weight: 700; color: #284F78; background: #EEF6FC;
      list-style: none; position: relative;
    }
    .example-dropdown summary::-webkit-details-marker { display: none; }
    .example-dropdown summary::after {
      content: ''; position: absolute; right: 20px; top: 50%;
      width: 0; height: 0; margin-top: -3px;
      border-left: 7px solid transparent; border-right: 7px solid transparent;
      border-top: 9px solid #284F78;
    }
    .example-dropdown[open] summary::after {
      border-top: 0; border-bottom: 9px solid #284F78;
    }
    .example-dropdown > div { padding: 16px 20px 20px; }
    .description ul { font-size: 18px; line-height: 1.55; margin: 12px auto 22px; max-width: 900px; }
    .description li { margin-bottom: 8px; }
    .intro-table {
      border-collapse: separate; border-spacing: 0; width: 100%;
      margin: 20px auto; overflow: hidden; border-radius: 8px;
      border: 1px solid #D9E8F6; font-size: 17px;
    }
    .intro-table th {
      background: #356291; color: white; font-weight: 700;
      padding: 12px 16px; text-align: center;
    }
    .intro-table td {
      padding: 12px 16px; text-align: left; vertical-align: middle;
      border-top: 1px solid #D9E8F6;
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
      font-size: 24px; font-weight: 700; color: #284F78;
    }
    .equation-pair { display: grid; gap: 18px; margin: 16px 0; }
    .equation-row {
      display: grid; grid-template-columns: 1fr 54px 1fr;
      align-items: stretch; gap: 12px;
    }
    .equation-arrow {
      display: flex; align-items: center; justify-content: center;
      color: #284F78; font-size: 28px; font-weight: 700;
    }
    .equation-box {
      border: 1px solid #D9E8F6; border-radius: 8px; padding: 14px;
      background: #F8FBFE;
    }
    .equation-box strong {
      display: block; color: #284F78; font-weight: 700; margin-bottom: 8px;
    }
    .calculation-stack {
      display: grid; gap: 14px; margin: 14px 0;
    }
    .calculation-step {
      border: 1px solid #D9E8F6; border-radius: 8px;
      background: #F8FBFE; padding: 14px 18px;
      overflow-x: auto;
      overflow-y: hidden;
      scrollbar-color: #94a3b8 #e2e8f0;
      scrollbar-width: thin;
    }
    .calculation-step::-webkit-scrollbar {
      height: 10px;
    }
    .calculation-step::-webkit-scrollbar-track {
      background: #e2e8f0;
      border-radius: 999px;
    }
    .calculation-step::-webkit-scrollbar-thumb {
      background: #94a3b8;
      border-radius: 999px;
    }
    .calculation-step strong {
      display: block; color: #284F78; font-weight: 700; margin-bottom: 8px;
      line-height: 1.4;
      white-space: normal;
    }
    .calculation-step .MathJax_Display,
    .calculation-step mjx-container[display='true'] {
      max-width: 100%;
      overflow-x: auto;
      overflow-y: hidden;
      padding-bottom: 4px;
    }
    .solution-grid {
      border-collapse: separate; border-spacing: 0; width: 100%;
      margin: 14px 0; border: 1px solid #D9E8F6; border-radius: 8px;
      overflow: hidden; font-size: 17px;
    }
    .solution-grid th {
      background: #EEF6FC; color: #284F78; font-weight: 700;
      padding: 10px 14px; border-bottom: 1px solid #D9E8F6;
    }
    .intro-table tr:first-child > *,
    .solution-grid tr:first-child > *,
    .matrix-grid tr:first-child > * {
      text-align: center;
    }
    .solution-grid td {
      padding: 10px 14px; text-align: center; border-top: 1px solid #D9E8F6;
    }
    .tolerance-highlight {
      display: inline-block; background: #fff4cc; color: #4f3900;
      border: 1px solid #e5c65f; border-radius: 6px; padding: 2px 8px;
      font-weight: 700; white-space: nowrap;
    }
    .calculator-jump {
      display: inline-block; margin-top: 10px; font-weight: 700;
      color: white; background: #356291; border-color: #356291;
    }
    .convergence-check {
      display: grid; grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 12px; margin: 14px 0;
    }
    .convergence-check .calculation-step {
      border-left: 5px solid #427AB5;
    }
    .dominance-grid {
      display: grid; grid-template-columns: 1fr auto 1fr auto 1fr;
      gap: 8px; align-items: center; margin-top: 10px;
    }
    .dominance-part {
      background: #EEF6FC; border: 1px solid #D9E8F6; border-radius: 6px;
      padding: 8px 10px; text-align: center;
    }
    .dominance-label {
      display: block; color: #284F78; font-size: 13px; font-weight: 700;
      margin-bottom: 4px;
    }
    .dominance-value { display: block; font-size: 18px; font-weight: 700; }
    .dominance-symbol { color: #284F78; font-size: 22px; font-weight: 700; }
    .dominance-pass .dominance-value { color: #356291; }
    .dominance-table td,
    .dominance-table th {
      text-align: center;
      vertical-align: middle;
    }
    .dominance-table .comparison-cell {
      color: #284F78;
      font-size: 22px;
      font-weight: 700;
    }
    .dominance-table .pass-cell {
      color: #356291;
      font-weight: 700;
    }
    .calculator-page {
      max-width: 1220px;
      margin: 0 auto;
    }
    .tw-calculator-page {
      position: relative;
    }
    .calculate-page-loader {
      align-items: center;
      background: rgba(248, 251, 254, 0.88);
      backdrop-filter: blur(3px);
      border-radius: 16px;
      bottom: 0;
      color: #284F78;
      display: none;
      flex-direction: column;
      gap: 12px;
      justify-content: center;
      left: 0;
      min-height: 220px;
      position: absolute;
      right: 0;
      top: 0;
      z-index: 20;
    }
    .calculate-page-loader.is-visible {
      display: flex;
    }
    .calculate-page-loader-bar {
      background: rgba(66, 122, 181, 0.14);
      border-radius: 999px;
      height: 8px;
      overflow: hidden;
      position: relative;
      width: min(360px, 72%);
    }
    .calculate-page-loader-bar::before {
      animation: calculator-loading-sweep 1s ease-in-out infinite;
      background: linear-gradient(90deg, transparent, #427AB5, transparent);
      content: '';
      height: 100%;
      left: -45%;
      position: absolute;
      top: 0;
      width: 45%;
    }
    .calculate-page-loader-label {
      font-size: 16px;
      font-weight: 700;
      letter-spacing: 0;
    }
    .calculator-hero {
      display: flex;
      align-items: flex-end;
      justify-content: space-between;
      gap: 18px;
      padding: 18px 0 8px;
      border-bottom: 1px solid #D9E8F6;
      margin-bottom: 18px;
    }
    .calculator-hero h3 {
      margin: 0 0 6px;
      color: #284F78;
      font-size: 28px;
      font-weight: 700;
    }
    .calculator-hero p {
      margin: 0;
      color: #566575;
      font-size: 18px;
      line-height: 1.5;
      max-width: 720px;
    }
    .calculator-badge {
      background: #EEF6FC;
      border: 1px solid #D9E8F6;
      color: #284F78;
      border-radius: 6px;
      padding: 8px 12px;
      font-size: 14px;
      font-weight: 700;
      white-space: nowrap;
    }
    .calculator-layout {
      gap: 18px;
      align-items: start;
    }
    .calculator-panel {
      border-radius: 8px;
      margin-bottom: 18px;
    }
    .calculator-panel h4 {
      color: #284F78;
      font-size: 20px;
      font-weight: 700;
      margin-top: 0;
    }
    .calculator-panel .help-block {
      color: #657487;
      font-size: 14px;
      margin: -4px 0 12px;
    }
    .calculator-control-row {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px;
    }
    .calculator-input-card {
      position: sticky;
      top: 14px;
    }
    .calculator-input-card .form-group label {
      color: #334155;
      font-size: 14px;
      font-weight: 700;
      letter-spacing: 0;
    }
    .tw-calculator-input-card .form-group label {
      color: #334155;
      font-size: 14px;
      font-weight: 700;
      letter-spacing: 0;
    }
    .tw-calculator-input-card .form-control {
      border: 1px solid #cbd5e1;
      border-radius: 10px;
      box-shadow: none;
      font-size: 18px;
      min-height: 46px;
    }
    .tw-calculator-input-card .form-control:focus {
      border-color: #356291;
      box-shadow: 0 0 0 3px rgba(66, 122, 181, 0.14);
    }
    .calculator-section-label {
      color: #284F78;
      font-size: 18px;
      font-weight: 700;
      margin: 18px 0 8px;
    }
    .calculator-grid-wrap {
      overflow-x: auto;
      padding-bottom: 4px;
      width: 100%;
      max-width: 100%;
    }
    .calculator-grid {
      display: grid;
      gap: 8px;
      width: max-content;
      max-width: 100%;
      margin-bottom: 10px;
      align-items: stretch;
    }
    .calculator-grid-cell,
    .calculator-grid-label {
      min-width: 0;
    }
    .calculator-grid-label {
      background: #EEF6FC;
      color: #284F78;
      font-size: 18px;
      font-weight: 700;
      text-align: center;
      border: 1px solid #D9E8F6;
      border-radius: 6px;
      padding: 9px 10px;
      display: flex;
      align-items: center;
      justify-content: center;
      min-height: 44px;
    }
    .calculator-grid-cell {
      border: 1px solid #D9E8F6;
      border-radius: 6px;
      padding: 7px;
      text-align: center;
      background: white;
    }
    .calculator-grid .form-group {
      margin-bottom: 0;
    }
    .calculator-grid .control-label {
      display: none;
    }
    .tw-calculator-grid .form-group {
      margin-bottom: 0;
    }
    .tw-calculator-grid .control-label {
      display: none;
    }
    .calculator-grid input.form-control {
      text-align: center;
      padding-left: 6px;
      padding-right: 6px;
      width: 100%;
      min-width: 0;
      box-shadow: none;
    }
    .calculator-grid input.form-control:focus {
      border-color: #427AB5;
      box-shadow: 0 0 0 2px rgba(74, 127, 184, 0.16);
    }
    .tw-calculator-grid input.form-control {
      text-align: center;
      padding-left: 6px;
      padding-right: 6px;
      width: 100%;
      min-width: 0;
      box-shadow: none;
      border-radius: 8px;
      font-size: 18px;
      font-weight: 500;
      min-height: 40px;
    }
    .tw-calculator-grid input.form-control:focus {
      border-color: #427AB5;
      box-shadow: 0 0 0 2px rgba(74, 127, 184, 0.16);
    }
    .tw-matrix-table-wrap {
      width: max-content;
      max-width: none;
      overflow: visible;
      border-radius: 14px;
      background: #f8fafc;
      padding: 8px;
    }
    .tw-matrix-table {
      border-collapse: collapse;
      width: max-content;
      table-layout: fixed;
      background: #f8fafc;
      border: 2px solid #b9cbe0;
    }
    .tw-matrix-table td {
      border: 2px solid #b9cbe0;
      padding: 0;
      text-align: center;
      vertical-align: middle;
    }
    .tw-matrix-table .form-group {
      margin-bottom: 0;
    }
    .tw-matrix-table .control-label {
      display: none;
    }
    .tw-matrix-table input.form-control {
      background: #f8fafc;
      border: 0;
      border-radius: 0;
      box-shadow: none;
      color: #1f2937;
      font-size: var(--cell-font-size, 24px);
      font-weight: 500;
      height: 100%;
      min-height: var(--cell-size, 60px);
      line-height: 1.2;
      padding: 0 4px;
      text-align: center;
      width: 100%;
    }
    .tw-matrix-table input.form-control[type='number'],
    .tw-vector-table input.form-control[type='number'] {
      -moz-appearance: textfield;
      appearance: textfield;
    }
    .tw-matrix-table input.form-control[type='number']::-webkit-inner-spin-button,
    .tw-matrix-table input.form-control[type='number']::-webkit-outer-spin-button,
    .tw-vector-table input.form-control[type='number']::-webkit-inner-spin-button,
    .tw-vector-table input.form-control[type='number']::-webkit-outer-spin-button {
      -webkit-appearance: none;
      margin: 0;
    }
    .tw-matrix-table input.form-control:focus {
      background: #ffffff;
      box-shadow: inset 0 0 0 3px rgba(66, 122, 181, 0.22);
      outline: none;
    }
    .tw-vector-row {
      display: grid;
      gap: 10px;
      grid-template-columns: repeat(auto-fit, minmax(88px, 1fr));
    }
    .tw-vector-cell {
      background: #ffffff;
      border: 1px solid #e2e8f0;
      border-radius: 12px;
      padding: 10px;
    }
    .tw-vector-cell-label {
      color: #475569;
      display: block;
      font-size: 14px;
      font-weight: 700;
      margin-bottom: 8px;
      text-align: center;
    }
    .tw-vector-cell .form-group {
      margin-bottom: 0;
    }
    .tw-vector-cell .control-label {
      display: none;
    }
    .tw-vector-cell input.form-control {
      border: 1px solid #cbd5e1;
      border-radius: 8px;
      box-shadow: none;
      font-size: 18px;
      text-align: center;
    }
    .tw-vector-cell input.form-control:focus {
      border-color: #356291;
      box-shadow: 0 0 0 3px rgba(66, 122, 181, 0.14);
    }
    .tw-vector-table-wrap {
      width: max-content;
      max-width: none;
      overflow: visible;
      border-radius: 12px;
      background: #f8fafc;
      padding: 8px;
    }
    .tw-vector-table {
      border-collapse: collapse;
      table-layout: fixed;
      width: max-content;
      background: #f8fafc;
    }
    .tw-vector-table th,
    .tw-vector-table td {
      border: 2px solid #cbd5e1;
      box-sizing: border-box;
      padding: 0;
      text-align: center;
      vertical-align: middle;
    }
    .tw-vector-table th {
      background: #f1f5f9;
      color: #475569;
      font-size: 14px;
      font-weight: 700;
      padding: 10px 14px;
    }
    .tw-vector-table .form-group {
      margin-bottom: 0;
    }
    .tw-vector-table .control-label {
      display: none;
    }
    .tw-vector-table input.form-control {
      background: #f8fafc;
      border: 0;
      border-radius: 0;
      box-sizing: border-box;
      box-shadow: none;
      color: #1f2937;
      font-size: var(--cell-font-size, 24px);
      font-weight: 500;
      height: 100%;
      min-height: var(--cell-size, 60px);
      line-height: 1.2;
      padding: 0 4px;
      text-align: center;
      width: 100%;
    }
    .tw-vector-table .tw-vector-cell-input {
      background: #f8fafc;
      border: 0;
      border-radius: 0;
      box-sizing: border-box;
      box-shadow: none;
      color: #1f2937;
      display: block;
      font-size: var(--cell-font-size, 24px);
      font-weight: 700;
      height: 100%;
      line-height: 1.2;
      margin: 0;
      min-height: var(--cell-size, 60px);
      padding: 0;
      text-align: center;
      width: 100%;
    }
    .tw-vector-table .tw-vector-cell-input:focus {
      background: #ffffff;
      box-shadow: inset 0 0 0 3px rgba(66, 122, 181, 0.22);
      outline: none;
    }
    .tw-vector-table .tw-vector-cell-input[type='number'] {
      -moz-appearance: textfield;
      appearance: textfield;
    }
    .tw-vector-table .tw-vector-cell-input[type='number']::-webkit-inner-spin-button,
    .tw-vector-table .tw-vector-cell-input[type='number']::-webkit-outer-spin-button {
      -webkit-appearance: none;
      margin: 0;
    }
    .tw-vector-table .tw-vector-cell-label-box {
      align-items: center;
      background: #eef2f7;
      box-sizing: border-box;
      color: #64748b;
      display: flex;
      font-size: var(--cell-font-size, 24px);
      font-weight: 700;
      height: 100%;
      justify-content: center;
      line-height: 1.2;
      min-height: var(--cell-size, 60px);
      padding: 0;
      text-align: center;
      width: 100%;
    }
    .tw-vector-table input.form-control:focus {
      background: #ffffff;
      box-shadow: inset 0 0 0 3px rgba(66, 122, 181, 0.22);
      outline: none;
    }
    .tw-vector-table-vertical th {
      min-width: 56px;
      width: 56px;
    }
    .tw-vector-table-vertical td {
      min-width: 76px;
      width: 76px;
      height: 64px;
    }
    .tw-unknown-vector {
      border-collapse: collapse;
      width: max-content;
      background: #f8fafc;
    }
    .tw-unknown-vector td {
      box-sizing: border-box;
      color: #64748b;
      font-size: var(--cell-font-size, 24px);
      font-weight: 700;
      line-height: 1.2;
      padding: 0;
    }
    .grid-pair {
      display: grid;
      gap: 18px;
    }
    .tw-story-step {
      overflow: hidden;
    }
    .tw-equation-inputs {
      display: flex;
      align-items: flex-end;
      gap: 12px;
      max-width: 100%;
      min-width: 0;
      overflow-x: auto;
      overflow-y: hidden;
      padding-bottom: 12px;
      position: relative;
      scrollbar-color: #94a3b8 #e2e8f0;
      scrollbar-width: thin;
      width: 100%;
    }
    .tw-equation-inputs::-webkit-scrollbar {
      height: 12px;
    }
    .tw-equation-inputs::-webkit-scrollbar-track {
      background: #e2e8f0;
      border-radius: 999px;
    }
    .tw-equation-inputs::-webkit-scrollbar-thumb {
      background: #94a3b8;
      border-radius: 999px;
    }
    .tw-vector-scroll-wrap {
      max-width: 100%;
      overflow-x: auto;
      overflow-y: hidden;
      scrollbar-color: #94a3b8 #e2e8f0;
      scrollbar-width: thin;
      width: 100%;
    }
    .tw-vector-scroll-wrap::-webkit-scrollbar {
      height: 10px;
    }
    .tw-vector-scroll-wrap::-webkit-scrollbar-track {
      background: #e2e8f0;
      border-radius: 999px;
    }
    .tw-vector-scroll-wrap::-webkit-scrollbar-thumb {
      background: #94a3b8;
      border-radius: 999px;
    }
    .tw-equation-part,
    .tw-equation-vector-part {
      align-items: center;
      display: flex;
      flex: 0 0 auto;
      flex-direction: column;
      padding-top: 44px;
      position: relative;
      width: max-content;
    }
    .tw-equation-sign {
      align-items: center;
      align-self: stretch;
      color: #475569;
      display: flex;
      flex: 0 0 auto;
      font-size: 24px;
      font-weight: 700;
      justify-content: center;
      padding-top: 44px;
    }
    .tw-equation-sign-spacer {
      display: none;
    }
    .tw-equation-sign-symbol {
      align-items: center;
      display: flex;
      flex: 1 1 auto;
      justify-content: center;
      padding: 0 4px;
    }
    .tw-equation-spacer-label {
      display: none;
    }
    .tw-calculator-section-label {
      align-items: center;
      color: #284F78;
      display: flex;
      height: 20px;
      justify-content: center;
      left: 50%;
      margin: 0;
      position: absolute;
      text-align: center;
      top: 16px;
      transform: translateX(-50%);
      white-space: nowrap;
      width: max-content;
    }
    .tw-equation-sign-spacer {
      display: none;
    }
    .calculator-actions {
      display: flex;
      gap: 10px;
      align-items: center;
      flex-wrap: wrap;
      margin-top: 0;
    }
    .calculator-actions .btn {
      margin: 0;
      display: inline-block;
    }
    .calculator-actions #calculate {
      background: #427AB5;
      border-color: #427AB5;
      color: #fff;
      font-weight: 700;
      min-width: 132px;
    }
    .tw-calculator-actions #calculate {
      background: #427AB5;
      border-color: #427AB5;
      color: #fff;
      font-size: 18px;
      font-weight: 700;
      min-height: 44px;
      min-width: 128px;
      margin: 0;
      border-radius: 10px;
      padding: 9px 18px;
    }
    .tw-calculator-actions #calculate:hover,
    .tw-calculator-actions #calculate:focus {
      background: #356291;
      border-color: #356291;
      color: #ffffff;
    }
    .calculator-actions #download_table {
      background: #ffffff;
      border-color: #b9cbe0;
      color: #284F78;
      font-weight: 700;
      margin-top: 0;
    }
    .tw-calculator-actions #download_table {
      background: #ffffff;
      border-color: #b9cbe0;
      color: #284F78;
      font-size: 18px;
      font-weight: 700;
      min-height: 44px;
      margin: 0;
      border-radius: 10px;
      padding: 9px 18px;
    }
    .calculator-outputs {
      min-width: 0;
      position: relative;
    }
    .calculation-results-loader {
      align-items: center;
      background: rgba(248, 251, 254, 0.86);
      backdrop-filter: blur(3px);
      border-radius: 16px;
      bottom: 0;
      color: #284F78;
      display: none;
      flex-direction: column;
      gap: 12px;
      justify-content: center;
      left: 0;
      min-height: 260px;
      position: absolute;
      right: 0;
      top: 34px;
      z-index: 15;
    }
    .calculation-results-loader.is-visible {
      display: flex;
    }
    .calculation-results-loader .calculator-loading-bar {
      display: block;
      max-width: 360px;
      width: min(360px, 72%);
    }
    .calculation-results-loader-label {
      font-size: 16px;
      font-weight: 700;
    }
    .calculator-result-grid {
      margin-bottom: 16px;
    }
    .calculator-result-grid .card {
      min-height: 145px;
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
      background: #F8FBFE;
      border: 1px solid #D9E8F6;
      border-radius: 8px;
      padding: 14px;
      min-height: 130px;
    }
    .result-card h4 {
      color: #284F78;
      font-weight: 700;
      margin-top: 0;
      font-size: 20px;
    }
    .calculator-result-card {
      background: #ffffff;
      color: #222;
      border: 1px solid #D9E8F6;
      border-radius: 8px;
      padding: 14px;
      margin: 0 0 16px;
      display: grid;
      grid-template-columns: minmax(340px, 1.2fr) minmax(220px, 0.8fr) minmax(280px, 1fr);
      gap: 12px;
      align-items: start;
    }
    .calculator-result-card h4 {
      color: #284F78;
      font-weight: 700;
      margin: 0 0 8px;
      font-size: 20px;
    }
    .calculator-result-card pre {
      font-size: 14px;
      margin-bottom: 8px;
      white-space: pre-wrap;
    }
    .calculator-result-item {
      min-width: 0;
    }
    .calculator-results-empty {
      color: #657487;
      font-size: 14px;
      line-height: 1.5;
    }
    .calculator-tabs {
      border-radius: 8px;
    }
    .calculator-tabs .tab-content {
      min-height: 540px;
      padding-top: 18px;
    }
    .chart-note {
      color: #475569;
      font-size: 16px;
      line-height: 1.55;
      margin: 16px 0 0;
      max-width: none;
      width: 100%;
      box-sizing: border-box;
    }
    .tw-chart-note {
      font-size: 16px;
      line-height: 1.65;
      max-width: none;
      width: 100%;
      box-sizing: border-box;
    }
    .calculator-loading-bar {
      background: rgba(66, 122, 181, 0.12);
      border-radius: 999px;
      display: none;
      grid-column: 1 / -1;
      height: 8px;
      overflow: hidden;
      position: relative;
      width: 100%;
    }
    .calculator-loading-bar.is-visible {
      display: block;
    }
    .calculator-loading-bar::before {
      animation: calculator-loading-sweep 1s ease-in-out infinite;
      background: linear-gradient(90deg, transparent, #427AB5, transparent);
      content: '';
      height: 100%;
      left: -45%;
      position: absolute;
      top: 0;
      width: 45%;
    }
    @keyframes calculator-loading-sweep {
      0% { left: -45%; }
      100% { left: 100%; }
    }
    .shiny-notification {
      border: 1px solid #D9E8F6;
      border-radius: 10px;
      box-shadow: 0 10px 30px rgba(66, 122, 181, 0.16);
      font-size: 16px;
      line-height: 1.45;
      padding: 14px 16px;
      right: 22px;
    }
    .shiny-notification-warning {
      background: #EEF6FC;
      color: #284F78;
    }
    .result-value-grid {
      border-collapse: separate;
      border-spacing: 0;
      min-width: 100%;
      width: max-content;
      border: 1px solid #e2e8f0;
      border-radius: 12px;
      overflow: hidden;
      background: #ffffff;
    }
    .result-value-grid th {
      background: #f8fafc;
      color: #475569;
      font-size: 14px;
      font-weight: 700;
      padding: 8px 10px;
      text-align: center;
      border-bottom: 1px solid #e2e8f0;
      white-space: nowrap;
    }
    .result-value-grid td {
      padding: 10px;
      text-align: center;
      font-size: 18px;
      font-weight: 700;
      color: #0f172a;
      border-left: 1px solid #e2e8f0;
      min-width: 96px;
      white-space: nowrap;
    }
    .result-value-grid td:first-child {
      border-left: 0;
    }
    .result-number {
      display: block;
      background: #f8fafc;
      border: 1px solid #e2e8f0;
      border-radius: 12px;
      padding: 12px 14px;
      font-size: 20px;
      font-weight: 700;
      color: #0f172a;
      text-align: center;
      overflow-wrap: anywhere;
    }
    .result-scroll {
      max-width: 100%;
      overflow-x: auto;
      overflow-y: hidden;
      padding-bottom: 4px;
      scrollbar-color: #94a3b8 #e2e8f0;
      scrollbar-width: thin;
      width: 100%;
    }
    .result-scroll::-webkit-scrollbar {
      height: 10px;
    }
    .result-scroll::-webkit-scrollbar-track {
      background: #e2e8f0;
      border-radius: 999px;
    }
    .result-scroll::-webkit-scrollbar-thumb {
      background: #94a3b8;
      border-radius: 999px;
    }
    .stale-result-notice {
      background: #EEF6FC;
      border: 1px solid #D9E8F6;
      border-radius: 12px;
      color: #284F78;
      font-size: 16px;
      font-weight: 700;
      line-height: 1.45;
      padding: 14px 16px;
    }
    .stale-result-note {
      color: #566575;
      display: block;
      font-size: 14px;
      font-weight: 500;
      margin-top: 4px;
    }
    .tw-result-card .calculator-results-empty {
      color: #64748b;
      font-size: 14px;
      line-height: 1.55;
    }
    .tw-calculator-tabs .nav-tabs {
      border-bottom: 1px solid #e2e8f0;
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-bottom: 14px;
    }
    .tw-calculator-tabs .nav-tabs > li {
      margin-bottom: -1px;
    }
    .tw-calculator-tabs .nav-tabs > li > a {
      border: 0;
      border-radius: 10px 10px 0 0;
      color: #475569;
      font-size: 18px;
      font-weight: 700;
      line-height: 1.2;
      padding: 10px 14px;
    }
    .tw-calculator-tabs .nav-tabs > li.active > a,
    .tw-calculator-tabs .nav-tabs > li.active > a:focus,
    .tw-calculator-tabs .nav-tabs > li.active > a:hover {
      background: #356291;
      border: 0;
      color: #ffffff;
    }
    .tw-calculator-tabs .tab-content {
      min-height: 560px;
      padding-top: 8px;
    }
    @media (max-width: 1000px) {
      .result-summary { grid-template-columns: 1fr; }
      .grid-pair { grid-template-columns: 1fr; }
      .calculator-input-card { position: static; }
      .calculator-result-card { grid-template-columns: 1fr; }
    }
    @media (max-width: 800px) {
      .calculator-hero { display: block; }
      .calculator-badge { display: inline-block; margin-top: 12px; }
      .calculator-control-row { grid-template-columns: 1fr; }
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
      background-color: #427AB5; color: white; border-radius: 9px 9px 0 0;
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
      color: white; cursor: default; background-color: #427AB5;
      border: 1px solid #ddd; border-bottom-color: transparent;
    }
    .nav-tabs > li > a {
      color: #356291;
    }
    .nav-tabs > li > a:hover,
    .nav-tabs > li > a:focus {
      color: #284F78;
      background-color: #EEF6FC;
      border-color: #D9E8F6;
    }
    .tw-calculator-tabs .nav-tabs > li.active > a,
    .tw-calculator-tabs .nav-tabs > li.active > a:focus,
    .tw-calculator-tabs .nav-tabs > li.active > a:hover {
      background: #427AB5;
      border-color: #427AB5;
      color: #ffffff;
    }
    .status-ok { color: #356291; font-weight: bold; }
    .status-warn { color: #b8860b; font-weight: bold; }
    .status-bad { color: #b22222; font-weight: bold; }
    .tw-calculator-card-title { font-size: 19px; line-height: 1.3; }
    .tw-story-step-title { font-size: 19px; line-height: 1.35; }
    .tw-result-card-title { font-size: 14px; line-height: 1.4; letter-spacing: .04em; }
    .tw-results-heading { font-size: 19px; line-height: 1.3; }
    .tw-solve-title { font-size: 19px; line-height: 1.35; }
    .tw-calculator-subtitle { font-size: 14px; line-height: 1.55; }
    .tw-calculator-help { font-size: 14px; line-height: 1.5; }
    .tw-story-step-copy { font-size: 14px; line-height: 1.55; }
    .tw-solve-copy { font-size: 14px; line-height: 1.55; }
    .tw-calculator-section-label { font-size: 12px; line-height: 1.4; }
  ")),

  titlePanel(
    div(class = "title",
      div("Gauss-Jacobi Iterative Method"),
      div(class = "title-credit",
        span(class = "title-credit-label", "Prepared by"),
        span("Catanpatan · Apos · Clarit · Vicen · Capoy")
      )
    )
  ),

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
                      <strong>Compute x<sub>1</sub> for iteration 1</strong>
                      \\[x_1^{(1)} = \\dfrac{6 + 0 - 2(0)}{10}\\]
                      \\[x_1^{(1)} = \\dfrac{6}{10}\\]
                      \\[x_1^{(1)} = 0.6000\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute x<sub>2</sub> for iteration 1</strong>
                      \\[x_2^{(1)} = \\dfrac{25 + 0 + 0}{11}\\]
                      \\[x_2^{(1)} = \\dfrac{25}{11}\\]
                      \\[x_2^{(1)} = 2.2727\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute x<sub>3</sub> for iteration 1</strong>
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
                      <strong>Compute x<sub>1</sub> for iteration 2</strong>
                      \\[x_1^{(2)} = \\dfrac{6 + x_2^{(1)} - 2x_3^{(1)}}{10}\\]
                      \\[x_1^{(2)} = \\dfrac{6 + 2.2727 - 2(-1.1000)}{10}\\]
                      \\[x_1^{(2)} = \\dfrac{10.4727}{10}\\]
                      \\[x_1^{(2)} = 1.0473\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute x<sub>2</sub> for iteration 2</strong>
                      \\[x_2^{(2)} = \\dfrac{25 + x_1^{(1)} + x_3^{(1)}}{11}\\]
                      \\[x_2^{(2)} = \\dfrac{25 + 0.6000 + (-1.1000)}{11}\\]
                      \\[x_2^{(2)} = \\dfrac{24.5000}{11}\\]
                      \\[x_2^{(2)} = 2.2273\\]
                    </div>
                    <div class='calculation-step'>
                      <strong>Compute x<sub>3</sub> for iteration 2</strong>
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
      div(class = "tw-calculator-page",
        div(id = "calculate_page_loader", class = "calculate-page-loader",
          div(class = "calculate-page-loader-bar"),
          div(class = "calculate-page-loader-label", "Loading calculator components...")
        ),
        div(class = "tw-calculator-hero",
          div(
            h3(class = "tw-calculator-title", "Build and Solve Ax = b"),
            p(class = "tw-calculator-subtitle", "Follow the setup from the size of the system, to the coefficient matrix, to the starting estimate, then inspect how the Jacobi method converges.")
          )
        ),
        div(class = "tw-calculator-layout",
          div(class = "tw-calculator-card tw-calculator-input-card",
            h4(class = "tw-calculator-card-title", "System Setup"),
            p(class = "tw-calculator-help", "The sample system is loaded by default. Change the size first, then fill in the matching matrix and vectors."),
            div(class = "tw-story-flow",
              div(class = "tw-story-step",
                div(class = "tw-story-step-header",
                  span(class = "tw-step-number", "1"),
                  div(
                    h5(class = "tw-story-step-title", "Choose the system size and stopping rules"),
                    p(class = "tw-story-step-copy", "The matrix, right-hand side vector, and initial guess all resize from this value.")
                  )
                ),
                div(class = "tw-settings-grid",
                  numericInput("n", "Matrix size (n x n)", value = 3, min = 2, step = 1),
                  numericInput("max_iter", "Max iterations", value = 50, min = 1, max = 1000, step = 1),
                  numericInput("tol", "Tolerance", value = 1e-6, min = 1e-15, step = 1e-6)
                )
              ),

              div(class = "tw-story-step",
                div(class = "tw-story-step-header",
                  span(class = "tw-step-number", "2"),
                  div(
                    h5(class = "tw-story-step-title", "Enter the linear system"),
                    p(class = "tw-story-step-copy", HTML("Fill the coefficient matrix <strong>A</strong> and the right-hand side vector <strong>b</strong>."))
                  )
                ),
                div(class = "tw-equation-inputs",
                  div(class = "tw-equation-part",
                    div(class = "tw-calculator-section-label", "Coefficient matrix A"),
                    uiOutput("matrix_grid_input")
                  ),
                  div(class = "tw-equation-vector-part",
                    div(class = "tw-equation-spacer-label"),
                    uiOutput("unknown_vector")
                  ),
                  div(class = "tw-equation-sign",
                    div(class = "tw-equation-sign-spacer"),
                    div(class = "tw-equation-sign-symbol", "=")
                  ),
                  div(class = "tw-equation-part",
                    div(class = "tw-calculator-section-label", "Right-hand side vector b"),
                    uiOutput("b_grid_input")
                  )
                )
              ),

              div(class = "tw-story-step",
                div(class = "tw-story-step-header",
                  span(class = "tw-step-number", "3"),
                  div(
                    h5(class = "tw-story-step-title", HTML("Set the starting estimate x<sup>(0)</sup>")),
                    p(class = "tw-story-step-copy", "These are the first values used by the Jacobi iteration.")
                  )
                ),
                uiOutput("x0_grid_input")
              ),

              div(class = "tw-solve-panel",
                div(
                  h5(class = "tw-solve-title", "Run the iteration"),
                  p(class = "tw-solve-copy", "After calculation, the results below show the solution, error, status, plots, steps, and table.")
                ),
                div(class = "tw-calculator-actions",
                  actionButton("calculate", "Calculate"),
                  downloadButton("download_table", "Download table")
                ),
                div(id = "calculator_loading_bar", class = "calculator-loading-bar")
              )
            )
          ),
          div(class = "calculator-outputs",
            div(id = "calculation_results_loader", class = "calculation-results-loader",
              div(class = "calculator-loading-bar is-visible"),
              div(class = "calculation-results-loader-label", "Calculating results...")
            ),
            h4(class = "tw-results-heading", "Results"),
            div(class = "tw-result-grid",
              div(class = "tw-result-card",
                h4(class = "tw-result-card-title", "Solution Vector"),
                uiOutput("answer")
              ),
              div(class = "tw-result-card",
                h4(class = "tw-result-card-title", "Final Error"),
                uiOutput("error_out")
              ),
              div(class = "tw-result-card",
                h4(class = "tw-result-card-title", "Status"),
                uiOutput("status")
              )
            ),
            div(class = "tw-calculator-tabs",
              tabsetPanel(
                tabPanel("Convergence Plot",
                plotOutput("plot", height = "500px"),
                div(class = "tw-chart-note",
                  HTML("This plot shows the <strong>error at each iteration</strong> on a logarithmic scale. A steady downward trend indicates the method is converging toward the true solution.")
                )
              ),
                tabPanel("Solution Plot",
                plotOutput("solution_plot", height = "500px"),
                div(class = "tw-chart-note",
                  "This plot shows how each variable changes across iterations. Lines that flatten out indicate that the values are stabilizing."
                )
              ),
                tabPanel("Steps", uiOutput("formatted_calculations")),
                tabPanel("Iteration Table",
                dataTableOutput("table"),
                div(class = "tw-chart-note",
                  HTML("Each row records one full Gauss-Jacobi iteration. The <strong>Error</strong> column is the maximum absolute change since the previous iterate.")
                )
                )
              )
            )
          )
        )
      )
    )
  ),

  tags$footer(
    strong("Numerical Methods Final Activity"),
    br(),
    strong("Gauss-Jacobi Iterative Method"),
    " - 2026"
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

matrix_ui_metrics <- function(n) {
  if (n <= 3) return(list(cell = 64, vector = 76, unknown = 76, font = 24))
  if (n <= 5) return(list(cell = 56, vector = 68, unknown = 68, font = 22))
  if (n <= 7) return(list(cell = 50, vector = 62, unknown = 62, font = 20))
  list(cell = 48, vector = 60, unknown = 60, font = 18)
}

server <- function(input, output, session) {
  observeEvent(input$go_to_calculator, {
    updateTabsetPanel(session, "main_tabs", selected = "Calculate")
    session$sendCustomMessage("scrollToTop", TRUE)
  })

  observeEvent(input$main_tabs, {
    if (!identical(input$main_tabs, "Calculate")) return()
    session$sendCustomMessage("calculatePageLoading", TRUE)
    session$onFlushed(function() {
      session$sendCustomMessage("calculatePageLoading", FALSE)
    }, once = TRUE)
  }, ignoreInit = TRUE)

  output$matrix_grid_input <- renderUI({
    n <- as.integer(input$n)
    metrics <- matrix_ui_metrics(n)
    cell_size <- metrics$cell
    cell_style <- paste0("width:", cell_size, "px; min-width:", cell_size, "px; height:", cell_size, "px;")

    rows <- lapply(seq_len(n), function(i) {
      cells <- lapply(seq_len(n), function(j) {
        id <- paste0("a_", i, "_", j)
        value <- isolate(input[[id]])
        if (is.null(value)) value <- default_A_value(i, j, n)
        tags$td(style = cell_style,
          numericInput(id, label = NULL, value = value, step = 1, width = "100%")
        )
      })
      tags$tr(cells)
    })

    div(class = "tw-matrix-table-wrap", style = paste0("--cell-font-size:", metrics$font, "px; --cell-size:", metrics$cell, "px;"),
      tags$table(class = "tw-matrix-table",
        do.call(tags$tbody, rows)
      )
    )
  })

  output$b_grid_input <- renderUI({
    n <- as.integer(input$n)
    metrics <- matrix_ui_metrics(n)
    cell_size <- metrics$vector
    cell_style <- paste0("width:", cell_size, "px; min-width:", cell_size, "px; height:", metrics$cell, "px;")
    rows <- lapply(seq_len(n), function(i) {
      id <- paste0("b_", i)
      value <- isolate(input[[id]])
      if (is.null(value)) value <- default_b_value(i, n)
      tags$tr(
        tags$td(style = cell_style,
          tags$input(
            id = id,
            class = "tw-vector-cell-input",
            type = "number",
            value = value,
            step = 1
          )
        )
      )
    })

    table_width <- metrics$vector

    div(class = "tw-vector-table-wrap", style = paste0("--cell-font-size:", metrics$font, "px; --cell-size:", metrics$cell, "px;"),
      tags$table(class = "tw-vector-table tw-vector-table-vertical", style = paste0("width:", table_width, "px;"),
        do.call(tags$tbody, rows)
      )
    )
  })

  output$unknown_vector <- renderUI({
    n <- as.integer(input$n)
    metrics <- matrix_ui_metrics(n)
    cell_style <- paste0(
      "height:", metrics$cell, "px;",
      "width:", metrics$vector, "px;",
      "min-width:", metrics$vector, "px;"
    )
    rows <- lapply(seq_len(n), function(i) {
      tags$tr(
        tags$td(style = cell_style,
          tags$div(class = "tw-vector-cell-label-box", HTML(paste0("x<sub>", i, "</sub>")))
        )
      )
    })

    table_width <- metrics$vector

    div(class = "tw-vector-table-wrap", style = paste0("--cell-font-size:", metrics$font, "px; --cell-size:", metrics$cell, "px;"),
      tags$table(class = "tw-vector-table tw-vector-table-vertical tw-unknown-vector", style = paste0("width:", table_width, "px;"),
        do.call(tags$tbody, rows)
      )
    )
  })

  output$x0_grid_input <- renderUI({
    n <- as.integer(input$n)
    metrics <- matrix_ui_metrics(n)
    cell_size <- metrics$vector
    cell_style <- paste0("width:", cell_size, "px; min-width:", cell_size, "px; height:", metrics$cell, "px;")
    headers <- lapply(seq_len(n), function(i) {
      tags$th(style = paste0("width:", cell_size, "px; min-width:", cell_size, "px;"),
        HTML(paste0("x<sub>", i, "</sub><sup>(0)</sup>"))
      )
    })
    cells <- lapply(seq_len(n), function(i) {
      id <- paste0("x0_", i)
      value <- isolate(input[[id]])
      if (is.null(value)) value <- default_x0_value(i, n)
      tags$td(style = cell_style,
        numericInput(id, label = NULL, value = value, step = 1, width = "100%")
      )
    })

    div(class = "tw-vector-table-wrap tw-vector-scroll-wrap", style = paste0("--cell-font-size:", metrics$font, "px; --cell-size:", metrics$cell, "px;"),
      tags$table(class = "tw-vector-table",
        tags$thead(tags$tr(headers)),
        tags$tbody(tags$tr(cells))
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

  current_input_signature <- reactive({
    n <- as.integer(input$n)
    if (is.null(n) || is.na(n) || n < 2) return(NULL)

    values <- c(input$n, input$tol, input$max_iter)
    values <- c(values, unlist(lapply(seq_len(n), function(i) {
      unlist(lapply(seq_len(n), function(j) input[[paste0("a_", i, "_", j)]]))
    }), use.names = FALSE))
    values <- c(values, unlist(lapply(seq_len(n), function(i) input[[paste0("b_", i)]]), use.names = FALSE))
    values <- c(values, unlist(lapply(seq_len(n), function(i) input[[paste0("x0_", i)]]), use.names = FALSE))

    paste(ifelse(is.null(values), "", values), collapse = "|")
  })

  result_store <- reactiveVal(NULL)
  result_signature <- reactiveVal(NULL)
  result_is_stale <- reactiveVal(FALSE)

  mark_result_stale_if_needed <- function() {
    req(!is.null(result_store()))
    sig <- isolate(current_input_signature())
    saved_sig <- isolate(result_signature())
    if (!is.null(saved_sig) && !identical(sig, saved_sig)) {
      result_is_stale(TRUE)
    }
  }

  observeEvent(input$calculator_inputs_changed, {
    req(!is.null(result_store()))
    result_is_stale(TRUE)
  }, ignoreInit = TRUE)

  observeEvent(list(input$n, input$tol, input$max_iter), {
    mark_result_stale_if_needed()
  }, ignoreInit = TRUE)

  stale_result_notice <- function() {
    div(class = "stale-result-notice",
      "Variables updated.",
      span(class = "stale-result-note", "Click Calculate again to refresh the results.")
    )
  }

  require_fresh_result <- function() {
    validate(
      need(!is.null(result_store()), "Run the calculator to show this output."),
      need(!isTRUE(result_is_stale()), "Variables updated. Click Calculate again to refresh the results.")
    )
  }

  # Validate inputs and run Gauss-Jacobi only when the Calculate button is clicked.
  observeEvent(input$calculate, {
    on.exit(session$onFlushed(function() {
      session$sendCustomMessage("calculationFinished", TRUE)
    }, once = TRUE), add = TRUE)

    validate(
      need(!is.null(input$n) && input$n >= 2, "n must be at least 2."),
      need(input$tol > 0, "Tolerance must be positive."),
      need(input$max_iter >= 1, "Max iterations must be at least 1.")
    )

    n <- as.integer(input$n)

    A <- read_matrix_grid(n)
    b <- read_vector_grid("b", n, "b")
    x0 <- read_vector_grid("x0", n, "Initial guess")
    validate(
      need(all(diag(A) != 0), "Diagonal entries of A must be non-zero (cannot divide by zero).")
    )

    res <- gauss_jacobi(A, b, x0, tol = input$tol, max_iter = input$max_iter)
    res$dom <- is_diagonally_dominant(A)
    res$A   <- A
    res$b   <- b
    res$tol <- input$tol

    result_store(res)
    result_signature(current_input_signature())
    result_is_stale(FALSE)
  })

  result <- reactive({
    req(result_store())
    result_store()
  })

  # Final solution vector
  output$answer <- renderUI({
    if (is.null(result_store())) {
      return(div(class = "calculator-results-empty", "Run the calculator to show the solution vector."))
    }
    if (isTRUE(result_is_stale())) {
      return(stale_result_notice())
    }
    res <- result()
    labels <- paste0("x<sub>", seq_along(res$solution), "</sub>")
    values <- format(round(res$solution, 8), nsmall = 6, trim = TRUE)

    HTML(paste0(
      "<div class='result-scroll'><table class='result-value-grid'><tr>",
      paste0("<th>", labels, "</th>", collapse = ""),
      "</tr><tr>",
      paste0("<td>", values, "</td>", collapse = ""),
      "</tr></table></div>"
    ))
  })

  # Final error
  output$error_out <- renderUI({
    if (is.null(result_store())) {
      return(div(class = "calculator-results-empty", "Waiting for a calculation."))
    }
    if (isTRUE(result_is_stale())) {
      return(stale_result_notice())
    }
    res <- result()
    HTML(paste0("<span class='result-number'>", formatC(res$final_err, format = "e", digits = 6), "</span>"))
  })

  # Convergence status banner
  output$status <- renderUI({
    if (is.null(result_store())) {
      return(div(class = "calculator-results-empty", "Convergence details will appear here."))
    }
    if (isTRUE(result_is_stale())) {
      return(stale_result_notice())
    }
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
    require_fresh_result()
    res <- result()
    df <- res$history[-1, , drop = FALSE]   # drop iteration 0 (NA error)

    ggplot(df, aes(x = Iteration, y = Error)) +
      geom_line(color = "#427AB5", size = 1.1) +
      geom_point(color = "#b22222", size = 2.5) +
      scale_y_log10() +
      geom_hline(yintercept = res$tol, linetype = "dashed", color = "#356291") +
      annotate("text",
               x = max(df$Iteration), y = res$tol,
               label = sprintf("tolerance = %.0e", res$tol),
               vjust = -0.5, hjust = 1, color = "#356291") +
      labs(title = "Convergence of Gauss-Jacobi Iteration",
           x = "Iteration",
           y = "Max |x_new - x_old|  (log scale)") +
      theme_minimal(base_size = 14)
  })

  # Solution plot (x values vs. iteration)
    output$solution_plot <- renderPlot({
      require_fresh_result()
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
    require_fresh_result()
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
    if (is.null(result_store())) {
      return(div(class = "calculator-results-empty", "Run the calculator to show the calculation steps."))
    }
    if (isTRUE(result_is_stale())) {
      return(stale_result_notice())
    }
    res <- result()
    boxes <- lapply(res$steps, function(txt) {
      div(class = "calculations_box", tags$pre(txt))
    })
    fluidRow(div(class = "formatted_calculations", do.call(tagList, boxes)))
  })
}


# ---------- Run the App ----------

shinyApp(ui, server)
