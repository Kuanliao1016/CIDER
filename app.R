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
      
      radioButtons(
        inputId = "Comparison",
        label = "Is your anchor comparing each PROM response to the baseline or to the previous time?",
        choices = c("Baseline" = "baseline", "Previous time" = "last")
      ),
      
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview", tableOutput("headdf")),
        tabPanel("Scatter Plot", plotOutput("scatter_plot")),
        tabPanel("Forest Plot + CI Table",
                 fluidRow(
                   column(6, plotOutput("forest_plot")),
                   column(6, tableOutput("ci_table"))
                 )
        )
      )
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
    df <- df %>%
      mutate(
        SSID = as.factor(SSID),
        time = as.integer(time),
        PROM = as.numeric(PROM)
      )
    
    if (input$Anchor_type == "static") {
      # Only process PGI_S columns
      if (!"PGI_S" %in% names(df)) stop("Column 'PGI_S' not found in your data.")
      df <- df %>%
        mutate(PGI_S = as.numeric(PGI_S)) %>%
        arrange(SSID, time) %>%
        group_by(SSID) %>%
        mutate(
          baseline_PROM = first(PROM),
          baseline_PGI_S = first(PGI_S),
          PROM_chg = PROM - lag(PROM),
          PROM_chg_from_bl = PROM - baseline_PROM,
          PGI_S_chg = PGI_S - lag(PGI_S),
          PGI_S_chg_from_bl = PGI_S - baseline_PGI_S
        ) %>%
        ungroup()
    } else if (input$Anchor_type == "change") {
      # Only process PGI_C columns
      if (!"PGI_C" %in% names(df)) stop("Column 'PGI_C' not found in your data.")
      df <- df %>%
        mutate(PGI_C = as.numeric(PGI_C)) %>%
        arrange(SSID, time) %>%
        group_by(SSID) %>%
        mutate(
          baseline_PROM = first(PROM),
          PROM_chg = PROM - lag(PROM),
          PROM_chg_from_bl = PROM - baseline_PROM,
          PGI_C_chg = PGI_C - lag(PGI_C)
        ) %>%
        ungroup()
    } else {
      stop("Unknown anchor type.")
    }
    df
  })
  
  
  
  model_fit <- reactive({
    req(processed_data(), input$PROM_type, input$Anchor_type, input$Comparison)
    df <- processed_data()
    

    if (input$PROM_type == "continuous" & input$Anchor_type == "change"& input$Comparison == "baseline") {
      lmer(PROM_chg_from_bl ~ PGI_C + (1 | SSID), data = df)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "change" & input$Comparison == "last") {
      lmer(PROM_chg ~ PGI_C + (1 | SSID), data = df)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change"& input$Comparison == "baseline") {
      lmer(PROM_chg_from_bl ~ as.factor(PGI_C) + (1 | SSID), data = df)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change"& input$Comparison == "last") {
      lmer(PROM_chg ~ as.factor(PGI_C) + (1 | SSID), data = df)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      lmer(PROM_chg_from_bl ~ PGI_S_chg_from_bl + (1 | SSID), data = df)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "last") {
      lmer(PROM_chg ~ PGI_S_chg + (1 | SSID), data = df)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      lmer(PROM_chg_from_bl ~ as.factor(PGI_S_chg_from_bl) + (1 | SSID), data = df)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "last") {
      lmer(PROM_chg ~ as.factor(PGI_S_chg) + (1 | SSID), data = df)
    } else {
      NULL
    }
  })
  output$headdf <- renderTable({
    req(processed_data())
    head(processed_data())
  })
  output$scatter_plot <- renderPlot({
    req(processed_data(), input$PROM_type, input$Anchor_type, input$Comparison)
    
    if (input$Anchor_type == "change"& input$Comparison == "baseline"){
      df <- processed_data()
      
      ggplot(df, aes(x = PGI_C, y = PROM_chg_from_bl)) +
        geom_point(color = "grey30") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Patient Global Impression of Change (PGI-C)",
             y = "Change in PROM from Baseline")
    } else if (input$Anchor_type == "change"& input$Comparison == "last"){
      df <- processed_data()
      
      ggplot(df, aes(x = PGI_C, y = PROM_chg)) +
        geom_point(color = "grey30") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Patient Global Impression of Change (PGI-C)",
             y = "Change in PROM from last response")      
    }
      else if (input$Anchor_type == "static" & input$Comparison == "baseline"){
      df <- processed_data()
      
      ggplot(df, aes(x = PGI_S_chg_from_bl, y = PROM_chg_from_bl)) +
        geom_point(color = "grey30") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Patient Global Impression of Severity (PGI-S)",
             y = "Change in PROM from Baseline")
      }else{
        df <- processed_data()
        
        ggplot(df, aes(x = PGI_S_chg, y = PROM_chg)) +
          geom_point(color = "grey30") +
          geom_smooth(method = "lm", se = TRUE, color = "red") +
          labs(x = "Patient Global Impression of Severity (PGI-S)",
               y = "Change in PROM from last response")  
    }

  })
  
  output$results_ui <- renderUI({
    req(input$run)
    tagList(
      
      tableOutput("ci_table")
    )
  })
  

  
  output$forest_plot <- renderPlot({
    req(input$run)
    req(model_fit(), input$PROM_type, input$Anchor_type)
    
    model <- model_fit()
    
    # Build new_data based on input
    new_data <- if (input$PROM_type == "continuous" & input$Anchor_type == "change") {
      data.frame(PGI_C = 1:7, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change") {
      data.frame(PGI_C = as.factor(1:7), SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = 1:7, SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = 1:7, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = as.factor(1:7), SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = as.factor(1:7), SSID = NA)
    }
    
    # Define bootstrapping function
    boot_fun <- function(fit) {
      predict(fit, newdata = new_data, re.form = NA)
    }
    
    set.seed(123)
    boot_res <- bootMer(model, boot_fun, nsim = 1000)
    cis <- apply(boot_res$t, 2, quantile, probs = c(0.025, 0.975))
    
    new_data$lower_ci <- cis[1, ]
    new_data$upper_ci <- cis[2, ]
    new_data$predicted <- predict(model, newdata = new_data, re.form = NA)
    
    plot_data <- if (input$Anchor_type == "change") {
      new_data %>% mutate(PGI = PGI_C)
    } else if (input$Anchor_type == "static" & input$Comparison == "baseline") {
      new_data %>% mutate(PGI = PGI_S_chg_from_bl)
    } else if (input$Anchor_type == "static" & input$Comparison == "last") {
      new_data %>% mutate(PGI = PGI_S_chg)
    }
    
    ggplot(plot_data, aes(x = predicted, y = as.factor(PGI))) +
      geom_point(color = "grey10", size = 2) +
      geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "grey30") +
      labs(
        x = "Predicted PROM Change",
        y = ifelse(input$Anchor_type == "change", "PGI-C", "PGI-S"),
        title = paste("Forest Plot of Expected PROM Change by",
                      ifelse(input$Anchor_type == "change", "PGI-C", "PGI-S"),
                      "(", input$PROM_type, ")")
      ) +
      theme_bw() +
      theme(
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        plot.title = element_text(face = "bold")
      )
  })
  
  output$ci_table <- renderTable({
    req(model_fit())
    model <- model_fit()
    
    # Define new_data as before
    new_data <- if (input$PROM_type == "continuous" & input$Anchor_type == "change") {
      data.frame(PGI_C = 1:7, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change") {
      data.frame(PGI_C = as.factor(1:7), SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = 1:7, SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = 1:7, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = as.factor(1:7), SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = as.factor(1:7), SSID = NA)
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
    
    if (input$Anchor_type == "change") {
      new_data %>%
        mutate(Expected_PROM = predicted) %>%
        arrange(desc(as.numeric(PGI_C))) %>%
        select(PGI_C, Expected_PROM, lower_ci, upper_ci)
    } else if (input$Anchor_type == "static" & input$Comparison == "baseline") {
      new_data %>%
        mutate(PGI_S = PGI_S_chg_from_bl,
               Expected_PROM = predicted) %>%
        arrange(desc(as.numeric(PGI_S))) %>%
        select(PGI_S, Expected_PROM, lower_ci, upper_ci)
    } else if (input$Anchor_type == "static" & input$Comparison == "last") {
      new_data %>%
        mutate(PGI_S = PGI_S_chg,
               Expected_PROM = predicted) %>%
        arrange(desc(as.numeric(PGI_S))) %>%
        select(PGI_S, Expected_PROM, lower_ci, upper_ci)
    }
  }, striped = TRUE)
}

shinyApp(ui, server)