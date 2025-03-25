library(shiny)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
    titlePanel("Bayesian Posterior: Poisson-Gamma Model"),
    
    sidebarLayout(
        sidebarPanel(
            h3("Click to draw your Gamma prior (λ ∈ [0, 50])"),
            br(),
            
            numericInput("sum_x", "Sum of Observed Counts (∑x):", min = 0, value = 20),
            numericInput("n", "Number of Observations (n):", min = 1, value = 10),
            
            checkboxInput("use_gamma", "Use Parametric Gamma Prior", value = FALSE),
            conditionalPanel(
                condition = "input.use_gamma == true",
                numericInput("a", "Gamma Alpha (a):", min = 0.1, value = 2),
                numericInput("b", "Gamma Beta (b):", min = 0.1, value = 1)
            ),
            
          
            div(
                style = "display: flex; justify-content: space-between;",
                actionButton("clear", "Clear Drawing"),
                actionButton("update_posterior", "Compute Posterior Distribution")
            ),
            wellPanel(h4("Summary of Posterior"), verbatimTextOutput("posteriorSummary"))
        ),
        
        mainPanel(
            plotOutput("drawPlot", click = "plot_click"),
            plotOutput("posteriorPlot")
        )
    )
)

server <- function(input, output, session) {
    values <- reactiveValues(
        data = data.frame(x = numeric(), y = numeric()),
        compute_posterior = FALSE,
        posterior_summary = "",
        prior_vals = NULL,
        posterior_manual = NULL
    )
    
    observeEvent(input$plot_click, {
        new_x <- input$plot_click$x
        new_y <- input$plot_click$y
        
        if (!is.null(new_x) && !is.null(new_y) && new_x >= 0 && new_x <= 50) {
            values$data <- rbind(values$data, data.frame(x = new_x, y = new_y))
        }
    })
    
    observeEvent(input$clear, {
        values$data <- data.frame(x = numeric(), y = numeric())
    })
    
    output$drawPlot <- renderPlot({
        df <- values$data
        
        p <- ggplot() +
            ggtitle("User-Drawn Gamma Prior") +
            xlab("λ") + ylab("Density") +
            theme_minimal() +
            xlim(0, 50) + ylim(0, 1) +
            geom_blank()
        
        if (nrow(df) > 1) {
            df <- df %>% arrange(x)
            p <- p +
                geom_line(data = df, aes(x, y), color = "blue", size = 1) +
                geom_point(data = df, aes(x, y), color = "black", size = 2)
        }
        
        p
    })
    
    observeEvent(input$update_posterior, {
        values$compute_posterior <- TRUE
        
        if (input$sum_x < 0 || input$n <= 0) {
            showModal(modalDialog("∑x ≥ 0 and n > 0 required.", easyClose = TRUE))
            return()
        }
        
        lambda_vals <- seq(0, 50, length.out = 500)
        gap <- lambda_vals[2] - lambda_vals[1]
        
        if (input$use_gamma) {
            if (input$a <= 0 || input$b <= 0) {
                showModal(modalDialog("a and b must be positive.", easyClose = TRUE))
                return()
            }
            
            values$prior_vals <- dgamma(lambda_vals, shape = input$a, rate = input$b)
            a_post <- input$a + input$sum_x
            b_post <- input$b + input$n
        } else {
            df <- values$data %>% arrange(x)
            
            if (nrow(df) <= 1) {
                showModal(modalDialog("Draw a valid prior with multiple points.", easyClose = TRUE))
                return()
            }
            
            df$y <- df$y / sum(df$y)
            prior_vals <- approx(df$x, df$y, xout = lambda_vals, rule = 2)$y
            prior_vals <- prior_vals / (gap * sum(prior_vals))
            values$prior_vals <- prior_vals
            
            mean_prior <- sum(lambda_vals * prior_vals)
            var_prior <- sum((lambda_vals - mean_prior)^2 * prior_vals)
            
            a_prior <- (mean_prior^2) / var_prior
            b_prior <- mean_prior / var_prior
            
            a_post <- a_prior + input$sum_x
            b_post <- b_prior + input$n
        }
        
        likelihood_vals <- lambda_vals^input$sum_x * exp(-input$n * lambda_vals)
        posterior <- values$prior_vals * likelihood_vals
        posterior <- posterior / (gap * sum(posterior))
        values$posterior_manual <- posterior
        
        mean_post <- sum(lambda_vals * values$posterior_manual) / sum(values$posterior_manual)
        var_post <- sum((lambda_vals - mean_post)^2 * values$posterior_manual) / sum(values$posterior_manual)
        
        sd_post <- sqrt(var_post)
        cdf_post <- cumsum(posterior) / sum(posterior)
        
        median_post <- lambda_vals[which.min(abs(cdf_post - 0.5))]
        mode_post <- lambda_vals[which.max(posterior)]
        lower <- lambda_vals[which.min(abs(cdf_post - 0.025))]
        upper <- lambda_vals[which.min(abs(cdf_post - 0.975))]
        
        values$posterior_summary <- paste0(
            "Posterior Alpha (a): ", round(a_post, 3), "\n",
            "Posterior Beta (b): ", round(b_post, 3), "\n",
            "Mean: ", round(mean_post, 3), "\n",
            "SD: ", round(sd_post, 6), "\n",
            "Median: ", round(median_post, 3), "\n",
            "Mode: ", round(mode_post, 3), "\n",
            "95% Credible Interval: [", round(lower, 3), ", ", round(upper, 3), "]"
        )
    })
    
    output$posteriorPlot <- renderPlot({
        if (!values$compute_posterior || is.null(values$prior_vals) || is.null(values$posterior_manual)) return(NULL)
        
        lambda_vals <- seq(0, 50, length.out = 500)
        df <- data.frame(
            lambda = rep(lambda_vals, 2),
            density = c(values$prior_vals, values$posterior_manual),
            distribution = rep(c("Prior", "Posterior"), each = length(lambda_vals))
        )
        
        ggplot(df, aes(x = lambda, y = density, color = distribution)) +
            geom_line(size = 1.2) +
            scale_color_manual(values = c("blue", "red")) +
            labs(title = "Prior & Posterior Distributions", x = "λ", y = "Density") +
            theme_minimal() +
            coord_cartesian(xlim = c(0, 50)) +
            ylim(0, max(df$density) * 1.1) +
            theme(legend.position = c(0.8, 0.9))
    })
    
    output$posteriorSummary <- renderText({
        if (!values$compute_posterior) return("Click 'Compute Posterior' to see summary.")
        values$posterior_summary
    })
}

shinyApp(ui = ui, server = server)
