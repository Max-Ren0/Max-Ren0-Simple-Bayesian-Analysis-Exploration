library(shiny)
library(ggplot2)
library(dplyr)
library(Bolstad)

ui <- fluidPage(
    titlePanel("Simple Bayesian Analysis Exploration"),
    
    sidebarLayout(
        sidebarPanel(
            selectInput("model_type", "Select Model:", choices = c("Binomial", "Poisson", "Normal")),
            uiOutput("model_ui"),
            
            div(style = "display:flex; justify-content: space-between;",
                actionButton("clear", "Clear Drawing"),
                actionButton("update_posterior", "Compute Posterior Distribution")
            ),
            br(),
            wellPanel(h4("Summary of Posterior"), verbatimTextOutput("posteriorSummary"))
        ),
        
        mainPanel(
            plotOutput("drawPlot", click = "plot_click"),
            plotOutput("posteriorPlot")
        )
    )
)

server <- function(input, output, session) {
    extract_quantile <- function(result, q) {
        if (!is.null(result$quantiles)) {
            return(result$quantiles[as.character(q)])
        } else {
            return(result$quantileFun(q))
        }
    }
   
    values <- reactiveValues(
        data = data.frame(x = numeric(), y = numeric()),
        compute_posterior = FALSE,
        prior_vals = NULL,
        posterior_manual = NULL,
        posterior_summary = NULL
    )
    
    observeEvent(input$plot_click, {
        new_x <- input$plot_click$x
        new_y <- input$plot_click$y
        
        valid_click <- switch(input$model_type,
                              "Binomial" = new_x >= 0 && new_x <= 1,
                              "Poisson" = new_x >= input$lambda_min && new_x <= input$lambda_max,
                              "Normal"  = new_x >= input$mu_min && new_x <= input$mu_max
        )
        
        if (!is.null(new_x) && !is.null(new_y) && valid_click) {
            values$data <- rbind(values$data, data.frame(x = new_x, y = new_y))
        }
    })
    
    reset_values <- function() {
        values$data <- data.frame(x = numeric(), y = numeric())
        values$compute_posterior <- FALSE
        values$posterior_summary <- NULL
        values$posterior_manual <- NULL
        values$prior_vals <- NULL
    }
    
    observeEvent(input$clear, { reset_values() })
    observeEvent(input$model_type, { reset_values() })
    
    
    output$model_ui <- renderUI({
        switch(input$model_type,
               "Binomial" = tagList(
                   numericInput("x", "Number of Successes (x):", value = 25, min = 0),
                   numericInput("N", "Total Trials (N):", value = 50, min = 1),
                   checkboxInput("use_beta", "Use Beta Prior", value = FALSE),
                   conditionalPanel(
                       condition = "input.use_beta == true",
                       numericInput("a", "Beta Alpha (a):", value = 1, min = 0.1),
                       numericInput("b", "Beta Beta (b):", value = 1, min = 0.1)
                   )
               ),
               "Poisson" = tagList(
                   numericInput("lambda_min", "Lambda Range Min:", value = 0.01),
                   numericInput("lambda_max", "Lambda Range Max:", value = 50),
                   numericInput("sum_x", "Sum of Observed Counts (∑x):", value = 20, min = 0),
                   numericInput("n", "Number of Observations (n):", value = 10, min = 1),
                   checkboxInput("use_gamma", "Use Parametric Gamma Prior", value = FALSE),
                   conditionalPanel(
                       condition = "input.use_gamma == true",
                       numericInput("a", "Gamma Alpha (a):", value = 2, min = 0.1),
                       numericInput("b", "Gamma Beta (b):", value = 1, min = 0.1)
                   )
               ),
               "Normal" = tagList(
                   numericInput("mu_min", "Mu Range Min:", value = -20),
                   numericInput("mu_max", "Mu Range Max:", value = 20),
                   numericInput("sum_x", "Sum of Observations (∑x):", value = 100),
                   numericInput("n", "Number of Observations (n):", value = 10),
                   numericInput("sigma", "Known SD (σ):", value = 2, min = 0.1),
                   checkboxInput("use_normal", "Use Parametric Normal Prior", value = FALSE),
                   conditionalPanel(
                       condition = "input.use_normal == true",
                       numericInput("mu0", "Prior Mean (μ₀):", value = 0),
                       numericInput("tau", "Prior SD (τ):", value = 5, min = 0.1)
                   )
               )
        )
    })
    
    
    output$posteriorSummary <- renderText({
        if (input$model_type == "Normal" && input$mu_min >= input$mu_max) {
            return("Error: Mu Range Min must be less than Mu Range Max.")
        }
        if (!values$compute_posterior) return("Click 'Compute Posterior Distribution' to generate summary.")
        values$posterior_summary
    })
    
    observeEvent(input$update_posterior, {
        if (input$model_type == "Normal" && input$mu_min >= input$mu_max) {
            showModal(modalDialog("Error: Mu Range Min must be less than Mu Range Max.", easyClose = TRUE))
            return()
        }
        values$compute_posterior <- TRUE
        
        if (input$model_type == "Binomial") {
            if (input$x < 0 || input$N <= 0 || input$x > input$N) {
                showModal(modalDialog("Error: x must be ≥ 0, N must be > 0, and x ≤ N.", easyClose = TRUE))
                return()
            }
            theta_vals <- seq(0, 1, length.out = 500)
            if (input$use_beta) {
                if (input$a <= 0 || input$b <= 0) {
                    showModal(modalDialog("Error: Beta Alpha (a) and Beta Beta (b) must be greater than 0.", easyClose = TRUE))
                    return()
                }
                values$prior_vals <- dbeta(theta_vals, input$a, input$b)
                result <- binobp(input$x, input$N, a = input$a, b = input$b, pi = theta_vals, plot = FALSE)
            } else {
                df <- values$data
                if (nrow(df) <= 1) {
                    showModal(modalDialog("Please draw a valid prior first.", easyClose = TRUE))
                    return()
                }
                df <- df %>% arrange(x)
                prior_vals <- approx(df$x, df$y / sum(df$y), xout = theta_vals, rule = 2)$y
                prior_vals <- pmax(prior_vals, 0)
                prior_vals <- prior_vals / sum(prior_vals)
                values$prior_vals <- prior_vals
                result <- binodp(input$x, input$N, pi = theta_vals, pi.prior = prior_vals, plot = FALSE)
            }
            values$posterior_manual <- result$posterior
            values$posterior_summary <- paste(c(
                paste("Mean:", round(result$mean, 3)),
                paste("SD:", round(sqrt(result$var), 6)),
                paste("Median:", round(extract_quantile(result, 0.5), 3)),
                paste("Mode:", round(theta_vals[which.max(result$posterior)], 3)),
                paste0("95% CI: [", 
                       round(extract_quantile(result, 0.025), 3), ", ", 
                       round(extract_quantile(result, 0.975), 3), "]")
            ), collapse = "\n")
            
        }
        
        if (input$model_type == "Poisson") {
            if (input$lambda_min < 0 || input$lambda_max <= 0 || input$lambda_min >= input$lambda_max || input$n <= 0 || input$sum_x < 0) {
                showModal(modalDialog("Please enter valid values: λ min ≥ 0, λ max > 0, λ min < λ max, ∑x ≥ 0, and n > 0.", easyClose = TRUE))
                return()
            }
            lambda_vals <- seq(input$lambda_min, input$lambda_max, length.out = 500)
            gap <- diff(lambda_vals)[1]
            if (input$use_gamma) {
                prior_vals <- dgamma(lambda_vals, input$a, input$b)
                prior_vals <- prior_vals / (gap * sum(prior_vals)) 
            } else {
                df <- values$data
                if (nrow(df) <= 1) {
                    showModal(modalDialog("Please draw a valid prior first.", easyClose = TRUE))
                    return()
                }
                
                df <- df %>% arrange(x)
                prior_vals <- approx(df$x, df$y / sum(df$y), xout = lambda_vals, rule = 2)$y
                prior_vals <- pmax(prior_vals, 0)
                prior_vals <- prior_vals / (gap * sum(prior_vals))
            }
            
            likelihood <- lambda_vals^input$sum_x * exp(-input$n * lambda_vals)
            
            posterior <- prior_vals * likelihood
            posterior <- posterior / (gap * sum(posterior))
            
            values$prior_vals <- prior_vals
            values$posterior_manual <- posterior
            
            mean_post <- sum(lambda_vals * posterior) * gap
            var_post <- sum((lambda_vals - mean_post)^2 * posterior) * gap
            cdf_post <- cumsum(posterior) * gap
            median_post <- lambda_vals[which.min(abs(cdf_post - 0.5))]
            mode_post <- lambda_vals[which.max(posterior)]
            ci_lower <- lambda_vals[which.min(abs(cdf_post - 0.025))]
            ci_upper <- lambda_vals[which.min(abs(cdf_post - 0.975))]
            
            values$posterior_summary <- paste(c(
                paste("Mean:", round(mean_post, 3)),
                paste("SD:", round(sqrt(var_post), 6)),
                paste("Median:", round(median_post, 3)),
                paste("Mode:", round(mode_post, 3)),
                paste0("95% CI: [", round(ci_lower, 3), ", ", round(ci_upper, 3), "]")
            ), collapse = "\n")
        }
        
        if (input$model_type == "Normal") {
            
            
            mu_vals <- seq(input$mu_min, input$mu_max, length.out = 500)
            gap <- diff(mu_vals)[1]
            x_bar <- input$sum_x / input$n
            se <- input$sigma / sqrt(input$n)
            if (input$use_normal) {
                values$prior_vals <- dnorm(mu_vals, input$mu0, input$tau)
                likelihood <- dnorm(mu_vals, x_bar, se)
            } else {
                df <- values$data
                if (nrow(df) <= 1) {
                    showModal(modalDialog("Please draw a valid prior first.", easyClose = TRUE))
                    return()
                }
                df <- df %>% arrange(x)
                prior_vals <- approx(df$x, df$y / sum(df$y), xout = mu_vals, rule = 2)$y
                prior_vals <- pmax(prior_vals, 0)
                prior_vals <- prior_vals / (gap * sum(prior_vals))
                values$prior_vals <- prior_vals
                likelihood <- dnorm(mu_vals, x_bar, se)
            }
            posterior <- values$prior_vals * likelihood
            posterior <- posterior / (gap * sum(posterior))
            values$posterior_manual <- posterior
            mean_post <- sum(mu_vals * posterior) * gap
            var_post <- sum((mu_vals - mean_post)^2 * posterior) * gap
            cdf_post <- cumsum(posterior) * gap
            values$posterior_summary <- paste(c(
                paste("Mean:", round(mean_post, 3)),
                paste("SD:", round(sqrt(var_post), 6)),
                paste("Median:", round(mu_vals[which.min(abs(cdf_post - 0.5))], 3)),
                paste("Mode:", round(mu_vals[which.max(posterior)], 3)),
                paste0("95% CI: [", 
                       round(mu_vals[which.min(abs(cdf_post - 0.025))], 3), ", ",
                       round(mu_vals[which.min(abs(cdf_post - 0.975))], 3), "]")
            ), collapse = "\n")
            
        }
    })
    
    output$drawPlot <- renderPlot({
        df <- values$data
        x_range <- switch(input$model_type,
                          "Binomial" = c(0, 1),
                          "Poisson" = c(input$lambda_min, input$lambda_max),
                          "Normal" = c(input$mu_min, input$mu_max))
        
        x_label <- switch(input$model_type,
                          "Binomial" = expression(theta),
                          "Poisson" = expression(lambda),
                          "Normal" = expression(mu))
        
        y_range <- c(0, 1)
        
        ggplot() +
            geom_blank(data = data.frame(x = x_range, y = y_range), aes(x, y)) +
            scale_x_continuous(
                limits = x_range,
                breaks = unique(c(x_range[1], pretty(x_range, n = 5), x_range[2])),
                expand = c(0, 0)
            )+
            coord_cartesian(ylim = y_range)+
            labs(title = "User-Drawn Prior Distribution", x = x_label, y = "Density") +
            theme_minimal(base_size = 16) +
            theme(
                plot.margin = margin(10, 10, 10, 10),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, face = "bold")
            ) +
            {
                if (nrow(df) > 1) {
                    df <- df %>%
                        arrange(x) %>%
                        group_by(x) %>%
                        summarise(y = mean(y), .groups = "drop")
                    list(
                        geom_line(data = df, aes(x, y), color = "blue", linewidth = 1.2),
                        geom_point(data = df, aes(x, y), color = "black", size = 2.5)
                    )
                }
            }
    })
    
    
    output$posteriorPlot <- renderPlot({
        if (is.null(values$posterior_manual)) return(NULL)
        
        x_vals <- switch(input$model_type,
                         "Binomial" = seq(0, 1, length.out = 500),
                         "Poisson"  = seq(input$lambda_min, input$lambda_max, length.out = 500),
                         "Normal"   = seq(input$mu_min, input$mu_max, length.out = 500)
        )
        
        df <- data.frame(
            x = rep(x_vals, 2),
            density = c(values$prior_vals, values$posterior_manual),
            dist = rep(c("Prior", "Posterior"), each = length(x_vals))
        )
        
        x_label <- switch(input$model_type,
                          "Binomial" = expression(theta),
                          "Poisson" = expression(lambda),
                          "Normal" = expression(mu))
        
        ggplot(df, aes(x, density, color = dist)) +
            scale_x_continuous(
                limits = range(x_vals),
                breaks = unique(c(x_vals[1], pretty(x_vals, n = 5), x_vals[length(x_vals)])),
                expand = c(0, 0)
            ) +
            geom_line(size = 1.5) +
            labs(title = "Prior & Posterior Distributions", x = x_label, y = "Density") +
            theme_minimal(base_size = 16) +
            theme(
                legend.position = c(0.95, 0.95),
                legend.text = element_text(size = 16),
                legend.title = element_blank(),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, face = "bold")
            )
        
        
    })
}

shinyApp(ui = ui, server = server)