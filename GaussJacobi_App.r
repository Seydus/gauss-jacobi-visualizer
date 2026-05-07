library(shiny)
library(DT)
library(ggplot2)

# ---------- Helper Functions ----------

# Parse a space- or comma-separated string into a numeric vector
parse_vector <- function(text) {
  text <- gsub(",", " ", text)
  vals <- suppressWarnings(as.numeric(strsplit(trimws(text), "\\s+")[[1]]))
  vals[!is.na(vals)]
}

# Parse a multi-line string into an n x n numeric matrix
parse_matrix <- function(text, n) {
  lines <- strsplit(trimws(text), "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) != n)
    stop(paste("Matrix A must have exactly", n, "rows. You provided", length(lines), "."))
  mat <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    row <- parse_vector(lines[i])
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
    HTML("<link href='https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap' rel='stylesheet'>")
  ),
  tags$style(HTML("
    * { font-family: Roboto, Arial, sans-serif; font-weight: 400; }
    .title { text-align: center; }
    footer { text-align: center; padding: 20px 0; }
    .btn-default {
      color: white; background-color: #4a7fb8; border-color: transparent;
      margin: 0 auto; display: block;
    }
    pre#answer, pre#error_out, pre#status {
      background: white; font-size: 14px;
    }
    .tab-content { padding-top: 20px; min-height: 650px; }
    h4 { font-weight: 400; }
    .center { font-weight: bold; text-align: center; margin: 30px auto 0; }
    .description { padding: 10px 100px 20px; text-align: left; }
    .description table td, .description table th { padding: 10px 20px; text-align: center; }
    .description table th { font-weight: bold; }
    .description table { margin: 20px auto; }
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
  ")),

  titlePanel(div("Gauss-Jacobi Iterative Method", class = "title")),

  tabsetPanel(

    # ----- Introduction Tab -----
    tabPanel("Introduction",
      fluidRow(
        div(class = "description",

          h4(class = "center", "Definition"),
          h4(HTML("The <strong>Gauss-Jacobi Method</strong> is an iterative numerical technique
                   used to solve a system of linear equations of the form
                   <em>Ax = b</em>. It begins with an initial guess for the unknowns and
                   refines that guess iteratively using the equations of the system,
                   until the solution converges within a chosen tolerance.")),

          h4("For an n x n linear system, each unknown is updated using the formula:"),
          withMathJax(h4(class = "center",
            "$$x_i^{(k+1)} = \\frac{1}{a_{ii}} \\left( b_i - \\sum_{j \\neq i} a_{ij} \\, x_j^{(k)} \\right)$$"
          )),

          h4("The error between successive iterations is measured as:"),
          withMathJax(h4(class = "center",
            "$$\\varepsilon^{(k)} = \\max_{1 \\le i \\le n} \\left| x_i^{(k+1)} - x_i^{(k)} \\right|$$"
          )),

          h4(HTML("A sufficient condition for convergence is that <strong>A</strong> is
                   <strong>strictly diagonally dominant</strong>, meaning for every row:
                   <em>|a<sub>ii</sub>| &gt; &Sigma;<sub>j&ne;i</sub> |a<sub>ij</sub>|</em>.")),

          h4(class = "center", "Applications"),
          h4(HTML("<strong>Engineering Systems.</strong>
                   Solving large, sparse linear systems arising in structural and circuit analysis.
                   <br><br><strong>Computational Physics.</strong>
                   Discretized partial differential equations (e.g. Laplace, Poisson) yield
                   linear systems that are well-suited to Jacobi iteration.
                   <br><br><strong>Parallel Computing.</strong>
                   Each component x<sub>i</sub><sup>(k+1)</sup> depends only on the previous
                   iterate, so updates are independent and can be computed in parallel.")),

          h4(class = "center", "Input Format Rules"),
          HTML("<table border='1'>
                  <tr><th>Field</th><th>Format</th><th>Example</th></tr>
                  <tr><td>Matrix A</td><td>One row per line, values separated by spaces or commas</td>
                      <td>10 -1 2<br>-1 11 -1<br>2 -1 10</td></tr>
                  <tr><td>Vector b</td><td>Single line, space- or comma-separated</td>
                      <td>6 25 -11</td></tr>
                  <tr><td>Initial Guess</td><td>Single line, space- or comma-separated</td>
                      <td>0 0 0</td></tr>
                  <tr><td>Tolerance</td><td>Positive real number</td><td>1e-6</td></tr>
                  <tr><td>Max Iterations</td><td>Positive integer</td><td>50</td></tr>
               </table>")
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
            numericInput("n", "Number of unknowns (n):", value = 3, min = 2, max = 8, step = 1),
            textAreaInput("matA", "Matrix A (one row per line):",
                          value = "10 -1 2\n-1 11 -1\n2 -1 10",
                          rows = 4, resize = "vertical"),
            textInput("vecB", "Vector b:", value = "6 25 -11"),
            textInput("x0",   "Initial guess x0:", value = "0 0 0"),
            numericInput("tol", "Tolerance:", value = 1e-6, min = 1e-15, step = 1e-6),
            numericInput("max_iter", "Max iterations:", value = 50, min = 1, max = 1000, step = 1)
          ),
          column(width = 12,
            mainPanel(width = 12,
              h4(style = "font-size: 1.2em;", HTML("Approximated solution vector:")),
              verbatimTextOutput("answer"),
              h4(style = "font-size: 1.2em;", HTML("Final approximation error:")),
              verbatimTextOutput("error_out"),
              h4(style = "font-size: 1.2em;", HTML("Convergence status:")),
              uiOutput("status")
            )
          )
        ),

        column(width = 8,
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
            tabPanel("Steps", uiOutput("formatted_calculations")),
            tabPanel("Iteration Table",
              dataTableOutput("table"),
              h4(style = "font-size: 1.1em; text-align: justify;",
                 HTML(paste(
                   "<div style='max-width: 900px; margin: 30px 20px;'>",
                   "Each row records one full Gauss-Jacobi iteration. The columns",
                   "<strong>x1, x2, ...</strong> are the components of the current estimate,",
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

server <- function(input, output, session) {

  # Validate inputs and run Gauss-Jacobi. Returns a list or stops with a clear message.
  result <- reactive({
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
    b  <- parse_vector(input$vecB)
    x0 <- parse_vector(input$x0)

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
  output$answer <- renderPrint({
    res <- result()
    setNames(round(res$solution, 8), paste0("x", seq_along(res$solution)))
  })

  # Final error
  output$error_out <- renderPrint({
    res <- result()
    res$final_err
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

  # Iteration table
  output$table <- renderDataTable({
    res <- result()
    datatable(res$history,
              options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE),
              rownames = FALSE) |>
      formatRound(columns = colnames(res$history)[-1], digits = 8)
  })

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
