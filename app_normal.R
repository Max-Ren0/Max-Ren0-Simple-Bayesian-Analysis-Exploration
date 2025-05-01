library(shiny)
library(ggplot2)
library(dplyr)
library(Bolstad)

ui <- fluidPage(
    titlePanel("Bayesian Posterior: Hand-Drawn or Parametric Prior with Custom μ-axis"),
    
    sidebarLayout(
        sidebarPanel(
            h3("Click to draw your Normal prior (μ ∈ [μ min, μ max])"),
            numericInput("mu_min", "μ Axis Min:", value = -20),
            numericInput("mu_max", "μ Axis Max:", value = 20),
            numericInput("sum_x", "Sum of Observations (∑x):", value = 100),
            numericInput("n", "Number of Observations (n):", value = 10),
            numericInput("sigma", "Known Standard Deviation (σ):", value = 2, min = 0.1),
            
            checkboxInput("use_normal", "Use Parametric Normal Prior", value = FALSE),
            conditionalPanel(
                condition = "input.use_normal == true",
                numericInput("mu0", "Prior Mean (μ₀):", value = 0),
                numericInput("tau", "Prior SD (τ):", value = 5, min = 0.1)
            ),
            
            div(style = "display:flex; justify-content: space-between;",
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
        compute_posterior = FALSE
    )
    
    # Use reactiveVal to store posterior results and summary text
    posterior_data <- reactiveVal(NULL)
    summary_data <- reactiveVal("Click 'Compute Posterior' to see summary.")
    
    # Capture user-drawn prior points
    observeEvent(input$plot_click, {
        if (!is.null(input$plot_click$x) &&
            input$plot_click$x >= input$mu_min &&
            input$plot_click$x <= input$mu_max) {
            values$data <- rbind(values$data, data.frame(x = input$plot_click$x, y = input$plot_click$y))
        }
    })
    
    # Clear drawing
    observeEvent(input$clear, {
        values$data <- data.frame(x = numeric(), y = numeric())
    })
    
    # Render user-drawn prior curve
    output$drawPlot <- renderPlot({
        df <- values$data
        p <- ggplot() +
            ggtitle("User-Drawn Prior") +
            xlab("μ") + ylab("Density") +
            xlim(input$mu_min, input$mu_max) + ylim(0, 1) +
            theme_minimal()
        
        if (nrow(df) > 1) {
            df <- df %>% arrange(x)
            p <- p +
                geom_line(data = df, aes(x, y), color = "blue", linewidth = 1) +
                geom_point(data = df, aes(x, y), color = "black", size = 2)
        }
        
        p
    })
    
    # Compute posterior from prior and data
    observeEvent(input$update_posterior, {
        values$compute_posterior <- TRUE
        
        n <- as.numeric(input$n)
        sum_x <- as.numeric(input$sum_x)
        sigma <- as.numeric(input$sigma)
        
        if (is.na(n) || is.na(sum_x) || is.na(sigma) || n <= 0 || sigma <= 0 || input$mu_min >= input$mu_max) {
            showModal(modalDialog("Please enter valid numeric values for n, ∑x, σ, and μ-axis range.", easyClose = TRUE))
            return()
        }
        
        x_bar <- sum_x / n
        
        if (input$use_normal) {
            mu0 <- as.numeric(input$mu0)
            tau <- as.numeric(input$tau)
            
            mu_vals <- seq(input$mu_min, input$mu_max, length.out = 500)
            gap <- mu_vals[2] - mu_vals[1]
            
            # Parametric normal prior
            prior_vals <- dnorm(mu_vals, mean = mu0, sd = tau)
            
            # Likelihood based on sample mean and known SE
            se <- sigma / sqrt(n)
            likelihood_vals <- dnorm(mu_vals, mean = x_bar, sd = se)
            
            # Posterior ∝ prior × likelihood
            posterior_vals <- prior_vals * likelihood_vals
            posterior_vals <- posterior_vals / (gap * sum(posterior_vals))  # Normalize
            
            # Posterior statistics
            mean_post <- sum(mu_vals * posterior_vals) * gap
            sd_post <- sqrt(sum((mu_vals - mean_post)^2 * posterior_vals) * gap)
            cdf_post <- cumsum(posterior_vals) * gap
            median_post <- mu_vals[which.min(abs(cdf_post - 0.5))]
            lower <- mu_vals[which.min(abs(cdf_post - 0.025))]
            upper <- mu_vals[which.min(abs(cdf_post - 0.975))]
            mode_post <- mu_vals[which.max(posterior_vals)]
        } else {
            mu_vals <- seq(input$mu_min, input$mu_max, length.out = 500)
            gap <- mu_vals[2] - mu_vals[1]
            
            df <- values$data
            if (nrow(df) <= 1) {
                showModal(modalDialog("Draw a valid prior with multiple points.", easyClose = TRUE))
                return()
            }
            
            df <- df[!duplicated(df$x), ]
            df <- df[order(df$x), ]
            df$y[is.na(df$y)] <- 0
            
            if (sum(df$y) == 0) {
                showModal(modalDialog("Prior values must contain non-zero weights.", easyClose = TRUE))
                return()
            }
            
            # Interpolate drawn prior to match mu_vals
            prior_vals <- approx(df$x, df$y, xout = mu_vals, rule = 2)$y
            prior_vals <- prior_vals / (gap * sum(prior_vals))
            
            # Compute likelihood and posterior
            likelihood_vals <- dnorm(mu_vals, mean = x_bar, sd = sigma / sqrt(n))
            posterior_vals <- prior_vals * likelihood_vals
            posterior_vals <- posterior_vals / (gap * sum(posterior_vals))
            
            # Posterior statistics
            mean_post <- sum(mu_vals * posterior_vals) * gap
            sd_post <- sqrt(sum((mu_vals - mean_post)^2 * posterior_vals) * gap)
            cdf_post <- cumsum(posterior_vals) * gap
            median_post <- mu_vals[which.min(abs(cdf_post - 0.5))]
            lower <- mu_vals[which.min(abs(cdf_post - 0.025))]
            upper <- mu_vals[which.min(abs(cdf_post - 0.975))]
            mode_post <- mu_vals[which.max(posterior_vals)]
        }
        
        # Store posterior curves and summary
        posterior_data(list(
            mu = mu_vals,
            prior = prior_vals,
            posterior = posterior_vals
        ))
        
        summary_lines <- c(
            paste0("Mean: ", round(as.numeric(mean_post), 3)),
            paste0("SD: ", round(as.numeric(sd_post), 6)),
            paste0("Median: ", round(as.numeric(median_post), 3)),
            paste0("Mode: ", round(as.numeric(mode_post), 3)),
            paste0("95% Credible Interval: [", round(as.numeric(lower), 3), ", ", round(as.numeric(upper), 3), "]")
        )
        summary_data(paste(summary_lines, collapse = "\n"))
    })
    
    # Plot posterior distribution
    output$posteriorPlot <- renderPlot({
        post <- posterior_data()
        if (is.null(post)) return(NULL)
        
        df <- data.frame(
            mu = rep(post$mu, 2),
            density = c(post$prior, post$posterior),
            distribution = rep(c("Prior", "Posterior"), each = length(post$mu))
        )
        
        ggplot(df, aes(x = mu, y = density, color = distribution)) +
            geom_line(linewidth = 1.2) +
            scale_color_manual(values = c("red", "blue")) +
            labs(title = "Prior & Posterior Distributions", x = "μ", y = "Density") +
            theme_minimal() +
            coord_cartesian(xlim = range(post$mu)) +
            ylim(0, max(df$density) * 1.1) +
            theme(legend.position.inside = c(0.8, 0.9))
    })
    
    # Show summary of posterior
    output$posteriorSummary <- renderText({
        summary_data()
    })
}

shinyApp(ui = ui, server = server)
