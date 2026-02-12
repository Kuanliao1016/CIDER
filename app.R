
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
      
      numericInput("reference", "Reference level for 'no change'", value = 4, min = -999, max = 999, step = 1),
      numericInput("improve_dir", "Direction of improvement (1 or -1)", value = 1, min = -1, max = 1, step = 2),
      
      actionButton("run", "Run Analysis"),
      
      # Added section for more info
      tags$hr(),
      tags$p("For more information, please visit: ",
             tags$a(href = "https://www.cider-bar.com/", "https://www.cider-bar.com/", target = "_blank"))
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview", tableOutput("headdf")),
        tabPanel("Scatter Plot", plotOutput("scatter_plot")),
        tabPanel("CID Within-group",
                 fluidRow(
                   column(6, plotOutput("forest_plot_within")),
                   column(6, tableOutput("ci_table_within"))
                 )),
        tabPanel("CID Between-Group",
                 fluidRow(
                   column(6, plotOutput("forest_plot_between")),
                   column(6, tableOutput("ci_table_between"))
                 )),
        tabPanel("Correlation", tableOutput("corr_table"))
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
    
    # Basic required columns
    if (!all(c("SSID", "time", "PROM") %in% names(df))) {
      stop("Required columns: SSID, time, PROM.")
    }
    
    df <- df %>%
      mutate(
        SSID = as.factor(SSID),
        time = as.integer(time),
        PROM = as.numeric(PROM)
      )
    
    if (input$Anchor_type == "static") {
      if (!"PGI_S" %in% names(df)) stop("Column 'PGI_S' not found in your data.")
      df <- df %>%
        mutate(PGI_S = as.numeric(PGI_S)) %>%
        arrange(SSID, time) %>%
        group_by(SSID) %>%
        mutate(
          baseline_PROM = coalesce(first(PROM), 0),
          baseline_PGI_S = coalesce(first(PGI_S), 0),
          PROM_chg = PROM - lag(PROM),
          PROM_chg_from_bl = PROM - baseline_PROM,
          PGI_S_chg = PGI_S - lag(PGI_S),
          PGI_S_chg_from_bl = PGI_S - baseline_PGI_S,
          CID_level = (PGI_S_chg_from_bl - input$reference) * input$improve_dir
        ) %>%
        ungroup()
      
    } else if (input$Anchor_type == "change") {
      if (!"PGI_C" %in% names(df)) stop("Column 'PGI_C' not found in your data.")
      df <- df %>%
        mutate(PGI_C = as.numeric(PGI_C)) %>%
        arrange(SSID, time) %>%
        group_by(SSID) %>%
        mutate(
          baseline_PROM = coalesce(first(PROM), 0),
          baseline_PGI_C = coalesce(first(PGI_C), 0),
          PROM_chg = PROM - lag(PROM),
          PROM_chg_from_bl = PROM - baseline_PROM,
          PGI_C_chg = PGI_C - lag(PGI_C),
          PGI_C_chg_from_bl = PGI_C - baseline_PGI_C,
          CID_level = (PGI_C_chg_from_bl - input$reference) * input$improve_dir
        ) %>%
        ungroup()
      
    } else {
      stop("Unknown anchor type.")
    }
    
    # Exclude baseline rows (where lag is NA)
    df <- df %>% filter(!is.na(PROM_chg))
    
    if (nrow(df) == 0) {
      stop("No follow-up observations found after baseline exclusion.")
    }
    
    df
  })
  
  model_fit <- reactive({
    req(processed_data(), input$PROM_type, input$Anchor_type, input$Comparison)
    df <- processed_data()
    
    if (input$PROM_type == "continuous" & input$Anchor_type == "change" & input$Comparison == "baseline") {
      lmer(PROM_chg_from_bl ~ PGI_C + (1 | SSID), data = df)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "change" & input$Comparison == "last") {
      lmer(PROM_chg ~ PGI_C + (1 | SSID), data = df)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change" & input$Comparison == "baseline") {
      lmer(PROM_chg_from_bl ~ as.factor(PGI_C) + (1 | SSID), data = df)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change" & input$Comparison == "last") {
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
    df <- processed_data()
    
    if (input$Anchor_type == "change" & input$Comparison == "baseline") {
      ggplot(df, aes(x = PGI_C, y = PROM_chg_from_bl)) +
        geom_point(color = "grey30") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Patient Global Impression of Change (PGI-C)",
             y = "Change in PROM from Baseline")
    } else if (input$Anchor_type == "change" & input$Comparison == "last") {
      ggplot(df, aes(x = PGI_C, y = PROM_chg)) +
        geom_point(color = "grey30") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Patient Global Impression of Change (PGI-C)",
             y = "Change in PROM from last response")
    } else if (input$Anchor_type == "static" & input$Comparison == "baseline") {
      ggplot(df, aes(x = PGI_S_chg_from_bl, y = PROM_chg_from_bl)) +
        geom_point(color = "grey30") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Patient Global Impression of Severity (PGI-S)",
             y = "Change in PROM from Baseline")
    } else {
      ggplot(df, aes(x = PGI_S_chg, y = PROM_chg)) +
        geom_point(color = "grey30") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = "Patient Global Impression of Severity (PGI-S)",
             y = "Change in PROM from last response")
    }
  })
  
  output$corr_table <- renderTable({
    req(processed_data())
    df <- processed_data()
    cor_data <- data.frame(
      CID_level = df$CID_level,
      PROM_change = if (input$Comparison == "baseline") df$PROM_chg_from_bl else df$PROM_chg
    )
    cor_val <- cor(cor_data$CID_level, cor_data$PROM_change, use = "pairwise.complete.obs")
    data.frame(
      Variable = "CID_level vs PROM_change",
      Correlation = round(cor_val, 4)
    )
  }, striped = TRUE)
  
  output$forest_plot_within <- renderPlot({
    req(input$run, model_fit(), input$PROM_type, input$Anchor_type)
    model <- model_fit()
    
    new_data <- if (input$PROM_type == "continuous" & input$Anchor_type == "change") {
      data.frame(PGI_C = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change") {
      data.frame(PGI_C = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = 1:4, SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = as.factor(1:4), SSID = NA)
    }
    
    boot_fun <- function(fit) predict(fit, newdata = new_data, re.form = NA)
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
    } else {
      new_data %>% mutate(PGI = PGI_S_chg)
    }
    
    ggplot(plot_data, aes(x = predicted, y = as.factor(PGI))) +
      geom_point(color = "grey10", size = 2) +
      geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "grey30") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        x = "Predicted PROM Change (Within-Group)",
        y = ifelse(input$Anchor_type == "change", "PGI-C", "PGI-S"),
        title = paste("Within-Group Forest Plot (", input$PROM_type, ")", sep = "")
      ) +
      theme_bw()
  })
  
  output$forest_plot_between <- renderPlot({
    req(input$run, model_fit(), input$PROM_type, input$Anchor_type)
    model <- model_fit()
    
    new_data <- if (input$PROM_type == "continuous" & input$Anchor_type == "change") {
      data.frame(PGI_C = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change") {
      data.frame(PGI_C = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = 1:4, SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = as.factor(1:4), SSID = NA)
    }
    
    boot_fun <- function(fit) predict(fit, newdata = new_data, re.form = NA)
    set.seed(123)
    boot_res <- bootMer(model, boot_fun, nsim = 1000)
    cis <- apply(boot_res$t, 2, quantile, probs = c(0.025, 0.975))
    
    new_data$lower_ci <- cis[1, ]
    new_data$upper_ci <- cis[2, ]
    new_data$predicted <- predict(model, newdata = new_data, re.form = NA)
    
    ref_idx <- which(if (input$Anchor_type == "change") new_data$PGI_C == input$reference else {
      if (input$Comparison == "baseline") new_data$PGI_S_chg_from_bl == input$reference else new_data$PGI_S_chg == input$reference
    })
    ref_val <- new_data$predicted[ref_idx][1]
    
    plot_data <- if (input$Anchor_type == "change") {
      new_data %>% mutate(PGI = PGI_C)
    } else if (input$Anchor_type == "static" & input$Comparison == "baseline") {
      new_data %>% mutate(PGI = PGI_S_chg_from_bl)
    } else {
      new_data %>% mutate(PGI = PGI_S_chg)
    }
    
    plot_data <- plot_data %>%
      mutate(diff_from_ref = predicted - ref_val,
             lower_diff = lower_ci - ref_val,
             upper_diff = upper_ci - ref_val) %>%
      filter(as.numeric(PGI) != input$reference)
    
    ggplot(plot_data, aes(x = diff_from_ref, y = as.factor(PGI))) +
      geom_point(color = "grey10", size = 2) +
      geom_errorbarh(aes(xmin = lower_diff, xmax = upper_diff), height = 0.2, color = "grey30") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        x = "Mean Difference from Reference (Between-Group)",
        y = ifelse(input$Anchor_type == "change", "PGI-C", "PGI-S"),
        title = paste("Between-Group Forest Plot (", input$PROM_type, ")", sep = "")
      ) +
      theme_bw()
  })
  
  output$ci_table_within <- renderTable({
    req(model_fit())
    model <- model_fit()
    
    new_data <- if (input$PROM_type == "continuous" & input$Anchor_type == "change") {
      data.frame(PGI_C = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change") {
      data.frame(PGI_C = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = 1:4, SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = as.factor(1:4), SSID = NA)
    }
    
    boot_fun <- function(fit) predict(fit, newdata = new_data, re.form = NA)
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
        mutate(PGI_S = PGI_S_chg_from_bl, Expected_PROM = predicted) %>%
        arrange(desc(as.numeric(PGI_S))) %>%
        select(PGI_S, Expected_PROM, lower_ci, upper_ci)
    } else {
      new_data %>%
        mutate(PGI_S = PGI_S_chg, Expected_PROM = predicted) %>%
        arrange(desc(as.numeric(PGI_S))) %>%
        select(PGI_S, Expected_PROM, lower_ci, upper_ci)
    }
  }, striped = TRUE)
  
  output$ci_table_between <- renderTable({
    req(model_fit(), input$reference)
    model <- model_fit()
    
    new_data <- if (input$PROM_type == "continuous" & input$Anchor_type == "change") {
      data.frame(PGI_C = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "change") {
      data.frame(PGI_C = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = 1:4, SSID = NA)
    } else if (input$PROM_type == "continuous" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = 1:4, SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "last") {
      data.frame(PGI_S_chg = as.factor(1:4), SSID = NA)
    } else if (input$PROM_type == "categorical" & input$Anchor_type == "static" & input$Comparison == "baseline") {
      data.frame(PGI_S_chg_from_bl = as.factor(1:4), SSID = NA)
    }
    
    boot_fun <- function(fit) predict(fit, newdata = new_data, re.form = NA)
    set.seed(123)
    boot_res <- bootMer(model, boot_fun, nsim = 1000)
    cis <- apply(boot_res$t, 2, quantile, probs = c(0.025, 0.975))
    
    new_data$lower_ci <- cis[1, ]
    new_data$upper_ci <- cis[2, ]
    new_data$predicted <- predict(model, newdata = new_data, re.form = NA)
    
    ref_idx <- which(if (input$Anchor_type == "change") new_data$PGI_C == input$reference else {
      if (input$Comparison == "baseline") new_data$PGI_S_chg_from_bl == input$reference else new_data$PGI_S_chg == input$reference
    })
    ref_val <- new_data$predicted[ref_idx][1]
    
    plot_data <- if (input$Anchor_type == "change") {
      new_data %>% mutate(PGI = PGI_C)
    } else if (input$Anchor_type == "static" & input$Comparison == "baseline") {
      new_data %>% mutate(PGI = PGI_S_chg_from_bl)
    } else {
      new_data %>% mutate(PGI = PGI_S_chg)
    }
    
    plot_data %>%
      mutate(diff_from_ref = predicted - ref_val,
             lower_diff = lower_ci - ref_val,
             upper_diff = upper_ci - ref_val) %>%
      filter(as.numeric(PGI) != input$reference) %>%
      select(PGI, diff_from_ref, lower_diff, upper_diff)
  }, striped = TRUE)
}

shinyApp(ui, server)