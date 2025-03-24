library(shiny)
library(ggplot2)
library(dplyr)

# UI Layout
ui <- fluidPage(
    titlePanel("Draw Prior & Compute Bayesian Posterior"),
    
    sidebarLayout(
        sidebarPanel(
            h3("Click to draw your prior distribution"),
            actionButton("clear", "Clear Drawing"),
            br(),
            
            # Input fields for observed data
            numericInput("x", "Number of Successes (x):", min = 1, max = 50, value = 25),
            numericInput("N", "Total Trials (N):", min = 1, max = 100, value = 50),
            
            # Checkbox to switch between drawn prior and Beta prior
            checkboxInput("use_beta", "Use Beta Prior (instead of drawn prior)", value = FALSE),
            
            # Beta prior parameters
            numericInput("a", "Prior Alpha (a):", min = 0.1, max = 10, value = 1),
            numericInput("b", "Prior Beta (b):", min = 0.1, max = 10, value = 1),
            
            # Button to update posterior distribution
            actionButton("update_posterior", "Compute Posterior Distribution"),
            
            # Summary moved to sidebarPanel
            wellPanel(h4("Summary of Posterior"), verbatimTextOutput("posteriorSummary"))
        ),
        
        mainPanel(
            plotOutput("drawPlot", click = "plot_click"),
            plotOutput("posteriorPlot")
        )
    )
)

# Server Logic
server <- function(input, output, session) {
    values <- reactiveValues(data = data.frame(x = numeric(), y = numeric()), compute_posterior = FALSE, posterior_summary = "", prior_vals = NULL, posterior_manual = NULL)
    
    observeEvent(input$plot_click, {
        new_x <- input$plot_click$x
        new_y <- input$plot_click$y
        
        if (!is.null(new_x) && !is.null(new_y) && new_x >= 0 && new_x <= 1) {  
            values$data <- rbind(values$data, data.frame(x = new_x, y = new_y))
        }
    })
    
    observeEvent(input$clear, {
        values$data <- data.frame(x = numeric(), y = numeric())
    })
    
    output$drawPlot <- renderPlot({
        df <- values$data
        
        p <- ggplot() +
            ggtitle("User-Drawn Prior Distribution") +
            xlab("θ") + ylab("Density") +
            theme_minimal() +
            xlim(0, 1) + ylim(0, 1) +  
            geom_blank()
        
        if (nrow(df) > 0) {
            df <- df %>% arrange(x)
            p <- p + 
                geom_line(data = df, aes(x, y), color = "blue", size = 1) + 
                geom_point(data = df, aes(x, y), color = "black", size = 2)
        }
        
        print(p)
    })
    
    observeEvent(input$update_posterior, {
        values$compute_posterior <- TRUE
        
        if (input$x <= 0 || input$N <= 0) {
            showModal(modalDialog("Error: x and N must be positive numbers.", easyClose = TRUE))
            return()
        }
        if (input$x > input$N) {
            showModal(modalDialog("Error: Number of successes (x) cannot be greater than total trials (N).", easyClose = TRUE))
            return()
        }
        
        
        theta_vals <- seq(0, 1, length.out = 500)
        gap <- theta_vals[2] - theta_vals[1]
        
        if (input$use_beta) {
            if (input$a <= 0 || input$b <= 0) {
                showModal(modalDialog("Error: a and b must be positive numbers.", easyClose = TRUE))
                return()
            }
            values$prior_vals <- dbeta(theta_vals, shape1 = input$a, shape2 = input$b)
            a_post <- input$a + input$x
            b_post <- input$b + input$N - input$x
            
        } else {
            df <- values$data %>% arrange(x)
            if (nrow(df) > 1) {
                df$y <- df$y / sum(df$y)
                values$prior_vals <- approx(df$x, df$y, xout = theta_vals, rule = 2)$y
                values$prior_vals <- values$prior_vals / (gap * sum(values$prior_vals))
            } else {
                showModal(modalDialog("Error: Please draw a valid prior distribution first.", easyClose = TRUE))
                return()
            }
            
            mean_prior <- sum(theta_vals * values$prior_vals) / sum(values$prior_vals)
            var_prior <- sum((theta_vals - mean_prior)^2 * values$prior_vals) / sum(values$prior_vals)
            
            a_post <- (mean_prior * (1 - mean_prior) / var_prior - 1) * mean_prior
            b_post <- (mean_prior * (1 - mean_prior) / var_prior - 1) * (1 - mean_prior)
        }
        
        likelihood_vals <- dbinom(input$x, size = input$N, prob = theta_vals)
        values$posterior_manual <- values$prior_vals * likelihood_vals
        values$posterior_manual <- values$posterior_manual / (gap * sum(values$posterior_manual))
        
        mean_post <- sum(theta_vals * values$posterior_manual) / sum(values$posterior_manual)
        var_post <- sum((theta_vals - mean_post)^2 * values$posterior_manual) / sum(values$posterior_manual)
        
        
        post_cdf <- cumsum(values$posterior_manual) / sum(values$posterior_manual)
        posterior_median <- theta_vals[which.min(abs(post_cdf - 0.5))]
        
        
        posterior_mode <- theta_vals[which.max(values$posterior_manual)]
        
        
        lower_bound <- theta_vals[which.min(abs(post_cdf - 0.025))]
        upper_bound <- theta_vals[which.min(abs(post_cdf - 0.975))]
        
        
        values$posterior_summary <- paste0(
            "Posterior Alpha (a): ", round(a_post, 3), "\n",
            "Posterior Beta (b): ", round(b_post, 3), "\n",
            "Mean: ", round(mean_post, 3), "\n",
            "Variance: ", round(var_post, 6), "\n",
            "Median: ", round(posterior_median, 3), "\n",
            "Mode: ", round(posterior_mode, 3), "\n",
            "95% Credible Interval: [", round(lower_bound, 3), ", ", round(upper_bound, 3), "]"
        )
        
    })
    
    output$posteriorPlot <- renderPlot({
        if (!values$compute_posterior || is.null(values$prior_vals) || is.null(values$posterior_manual)) return(NULL)
        
        theta_vals <- seq(0, 1, length.out = 500)
        df <- data.frame(
            theta_vals = rep(theta_vals, 2),
            density = c(values$prior_vals, values$posterior_manual),
            distribution = rep(c("Prior", "Posterior"), each = length(theta_vals))
        )
        
        ggplot(df, aes(x = theta_vals, y = density, color = distribution)) +
            geom_line(size = 1.2) +
            scale_color_manual(values = c("blue", "red")) +
            labs(title = "Prior & Posterior Distributions", x = "θ", y = "Density") +
            theme_minimal() +
            coord_cartesian(xlim = c(0, 1)) +  
            scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +  
            ylim(0, max(df$density) * 1.1) +  
            theme(legend.position = c(0.8,0.9))
    })
    
    output$posteriorSummary <- renderText({
        if (!values$compute_posterior) return("Click 'Compute Posterior Distribution' to generate summary.")
        values$posterior_summary
    })
}

shinyApp(ui = ui, server = server)
