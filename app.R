library(shiny)
library(dplyr)
library(lme4)
library(ggplot2)
library(boot)

########## UI ##########
ui <- fluidPage(
  titlePanel("CIDER - Clinically Important Difference Estimation and Rating using anchor-based approach"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("upload", "Upload CSV File", accept = ".csv"),
      
      radioButtons(
        inputId = "PROM_type",
        label = "If your PROM does not have a linear relationship with the anchor, consider treating the PROM as a categorical variable in the analysis.",
        choices = c("Continuous" = "continuous", "Categorical" = "categorical")
      ),
      
      radioButtons(
        inputId = "Anchor_type",
        label = "Is your anchor a static status (e.g. PGI-S) or global impression of change scale (e.g. PGI-C)?",
        choices = c("Static" = "static", "Change" = "change")
      ),
      
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      plotOutput("scatter_plot"),
      uiOutput("results_ui"),
      plotOutput("forest_plot")
    )
  )
)

########## Server ##########
server <- function(input, output) {
  
  raw_data <- reactive({
    req(input$upload)
    read.csv(input$upload$datapath)
  })
  
  processed_data <- eventReactive(input$run, {
    df <- raw_data()
    df %>%
      mutate(
        SSID = as.factor(SSID),
        time = as.integer(time),
        PGI_C = as.integer(PGI_C),
        PROM = as.integer(PROM)
      ) %>%
      group_by(SSID) %>%
      mutate(
        baseline = PROM[time == 0],
        PROM_chg_from_bl = PROM - baseline
      ) %>%
      ungroup()
  })
  
  model_fit <- reactive({
    req(processed_data(), input$PROM_type)
    df <- processed_data()
    
    if (input$PROM_type == "continuous") {
      lmer(PROM_chg_from_bl ~ PGI_C + (1 | SSID), data = df)
    } else {
      lmer(PROM_chg_from_bl ~ as.factor(PGI_C) + (1 | SSID), data = df)
    }
  })
  
  output$scatter_plot <- renderPlot({
    req(processed_data())
    df <- processed_data()
    
    ggplot(df, aes(x = PGI_C, y = PROM_chg_from_bl)) +
      geom_point(color = "grey30") +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(x = "Patient Global Impression of Change (PGI-C)",
           y = "Change in PROM from Baseline")
  })
  
  output$results_ui <- renderUI({
    req(input$run, input$PROM_type == "continuous")
    tagList(
      h4("Estimated Coefficient and 95% CI for PGI-C"),
      tableOutput("ci_table")
    )
  })
  
  output$ci_table <- renderTable({
    req(input$PROM_type == "continuous")
    model <- model_fit()
    coef_md1 <- summary(model)
    coefci_md1 <- confint(model, method = "boot")
    
    data.frame(
      coefficient = coef_md1$coefficients["PGI_C", "Estimate"],
      lower_ci = coefci_md1["PGI_C", 1],
      upper_ci = coefci_md1["PGI_C", 2]
    )
  }, rownames = FALSE)
  
  output$forest_plot <- renderPlot({
    req(input$run)
    req(model_fit())
    
    model <- model_fit()
    
    if (input$PROM_type == "continuous") {
      new_data <- data.frame(PGI_C = 1:7, SSID = NA)
    } else {
      new_data <- data.frame(PGI_C = as.factor(1:7), SSID = NA)
    }
    
    boot_fun <- function(fit) {
      predict(fit, newdata = new_data, re.form = NA)
    }
    
    set.seed(123)
    boot_res <- bootMer(model, boot_fun, nsim = 1000)
    cis <- apply(boot_res$t, 2, quantile, probs = c(0.025, 0.975))
    
    new_data$lower_ci <- cis[1, ]
    new_data$upper_ci <- cis[2, ]
    new_data$predicted <- predict(model, newdata = new_data, re.form = NA)
    
    ggplot(new_data, aes(x = predicted, y = as.factor(PGI_C))) +
      geom_point(color = "grey10", size = 2) +
      geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "grey30") +
      labs(
        x = "Predicted PROM Change",
        y = "PGI-C",
        title = paste("Forest Plot of Expected PROM Change by PGI-C (", input$PROM_type, ")")
      ) +
      theme_bw() +
      theme(
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.title = element_text(face = "bold")
      )
  })
}

shinyApp(ui, server)