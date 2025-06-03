# Load necessary packages
library(shiny)
library(DBI)
library(RMySQL)
library(DT)
library(dplyr)
library(openxlsx)

# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 8180)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-color: #2e2e2e;
        color: #f0f0f0;
      }
      .container {
        padding-top: 2em;
      }
      .sidebarPanel, .mainPanel {
        background-color: #3c3c3c;
        border: none;
        padding: 2em;
        border-radius: 8px;
      }
      h4 {
        color: #FFA500;
      }
      table.dataTable {
        background-color: #4a4a4a !important;
        color: #f0f0f0 !important;
      }
      .dataTables_wrapper .dataTables_filter input {
        background-color: #555 !important;
        color: #fff;
        border: 1px solid #888;
      }
      .dataTables_wrapper .dataTables_length select {
        background-color: #555 !important;
        color: #fff;
        border: 1px solid #888;
      }
      .stripe tbody tr.odd {
        background-color: #3c3c3c;
      }
      .stripe tbody tr.even {
        background-color: #4a4a4a;
      }
      .stripe tbody tr:hover {
        background-color: #FFA500 !important;
        color: #000 !important;
      }
      .dataTables_wrapper .dataTables_info {
        color: #ccc;
      }
      .dataTables_wrapper .dataTables_paginate .paginate_button {
        color: #FFA500 !important;
      }
      
      .btn, .form-control {
        background-color: #3a3a3a;
        color: #ffffff;
        border-color: #ff9900;
      }
      .btn:hover {
        background-color: #ff9900;
        color: #000;
      }
    "))
  ),
  sidebarLayout(
    sidebarPanel(
      # Add the logo at the top
      tags$img(src = "nsfw-logo.png", width = "100%", style = "margin-bottom: 1em;"),
      
      tags$h4("Examine variant effect:"),
      actionButton("showPerturb", "PerturbedDomain"),
      uiOutput("perturbInputs"),
      tags$hr(),
      uiOutput("returnButton"),
      tags$hr()
    ),
    mainPanel(
      DTOutput("mainTable"),
      tags$hr(),
      uiOutput("domainAnnotationUI"),
      
    )
  )
)

server <- function(input, output, session) {
  # 1. Connect to MySQL database when the server starts
  db_host <- Sys.getenv("DB_HOST", "127.0.0.1")  # fallback if env var missing
  db_user <- Sys.getenv("DB_USER", "root")
  db_password <- Sys.getenv("DB_PASSWORD", "")
  db_name <- Sys.getenv("DB_NAME", "InteractomeDB")
  db_port <- as.integer(Sys.getenv("DB_PORT", 3306))
  
  SFcon <- dbConnect(
    RMySQL::MySQL(),
    dbname = db_name,
    host = db_host,
    port = db_port,
    user = db_user,
    password = db_password
  )
  
  # 2. Load the local SURFME annotation file
  local_sf_data <- openxlsx::read.xlsx("meta/SURFME_v2023.xlsx")
  
  # 3. Query with JOIN from MySQL
  query <- "
    SELECT DISTINCT gi.gene_id, gi.gene_symbol, gi.chromosome, pi.uniprot_id, pi.protein_size
    FROM gene_info gi
    JOIN surfme_filter sf ON gi.gene_id = sf.gene_id
    JOIN protein_info pi ON gi.protein_id = pi.protein_id;
  "
  joined_data <- dbGetQuery(SFcon, query)
  
  # 4. Merge with local annotation by uniprot_id = Entry
  merged_data <- dplyr::left_join(joined_data, local_sf_data, 
                                  by = c("uniprot_id" = "Entry")) 
  print_data <- merged_data %>% 
    dplyr::select(uniprot_id, gene_id, gene_symbol, Protein.names1) %>% 
    dplyr::rename(
      UniprotID = uniprot_id,
      GeneID = gene_id,
      Symbol = gene_symbol,
      Annotation = Protein.names1)
  
  # Reactive values to track state
  rv <- reactiveValues(
    showPerturb = FALSE,
    perturbResult = NULL
  )
  # Show perturb input fields when PerturbedDomain is clicked
  observeEvent(input$showPerturb, {
    rv$showPerturb <- TRUE
  })
  
  # Reset view to default when Return is clicked
  observeEvent(input$returnToDefault, {
    rv$showPerturb <- FALSE
    rv$perturbResult <- NULL
  })
  # UI for GeneID and SNPposition
  output$perturbInputs <- renderUI({
    if (rv$showPerturb) {
      tagList(
        textInput("geneInput", "Gene ID", value = ""),
        numericInput("snpInput", "SNP Position", value = NA, min = 1),
        actionButton("runPerturb", "Run")
      )
    }
  })
  # Execute stored procedure on click
  observeEvent(input$runPerturb, {
    req(input$geneInput, input$snpInput)
    
    tryCatch({
      # Send the query
      result <- dbSendQuery(SFcon, paste0("CALL PerturbedDomains('", input$geneInput, "',", input$snpInput, ");"))
      
      # Fetch the first result set
      data <- dbFetch(result)
      
      # Clear this result set
      dbClearResult(result)
      
      # Force-read any remaining results (MySQL workaround)
      while (dbMoreResults(SFcon)) {
        suppressWarnings({
          next_result <- dbNextResult(SFcon)
          if (!is.null(next_result)) dbClearResult(next_result)
        })
      }
      
      # Save result
      rv$perturbResult <- data
    }, error = function(e) {
      rv$perturbResult <- data.frame(Error = e$message)
    })
  })
  # Render return button if perturb result is available
  output$returnButton <- renderUI({
    if (!is.null(rv$perturbResult)) {
      actionButton("returnToDefault", "Return")
    }
  })
  # Render main table or perturb result
  output$mainTable <- DT::renderDataTable({
    if (!is.null(rv$perturbResult)) {
      DT::datatable(
        rv$perturbResult,
        options = list(scrollX = TRUE, dom = 't'),
        class = "compact stripe hover",
        rownames = FALSE
      )
    } else {
      DT::datatable(
        print_data, 
        selection = "single",
        options = list(scrollX = TRUE, scrollY = "300px", paging = FALSE, dom = '<"top"lf>rt<"bottom"><"clear">'),
        style = "bootstrap",
        class = "compact stripe hover row-border order-column"
      )
    }
  })
  
  # 6. Render UI element under table upon row selection
  output$domainAnnotationUI <- renderUI({
    selected_row <- input$mainTable_rows_selected
    
    if (length(selected_row)) {
      selected_gene_id <- print_data[selected_row, "GeneID"]
      
      # Run the secondary query
      domain_query <- paste0("
        SELECT DISTINCT pdm.domain_id, pdm.domain_name, pdm.domain_start, pdm.domain_end, pdm.accuracy, pdm.E_value
        FROM gene_info gi
        JOIN protein_domain_map pdm ON gi.gene_id = pdm.gene_id
        WHERE gi.gene_id = '", selected_gene_id, "';")
      
      domain_data <- dbGetQuery(SFcon, domain_query, params = list(selected_gene_id))
      domain_data <- domain_data %>% 
        dplyr::group_by(domain_id) %>% 
        dplyr::arrange(E_value, .by_group = TRUE) %>%
        dplyr::slice(1)
      
      # Render as DT if any domains found
      if (nrow(domain_data) > 0) {
        # Save domain_data to reactive value
        output$domainTable <- DT::renderDataTable({
          DT::datatable(
            domain_data,
            rownames = FALSE,
            options = list(
              paging = FALSE,      # no pagination
              searching = FALSE,   # no search box
              info = FALSE,        # no 'Showing 1 to n of n entries'
              ordering = FALSE,     # no sorting
              dom = 't'            # only the table
            ),
            style = "bootstrap",
            class = "compact stripe hover",
            escape = FALSE
          )
        })
        
        tagList(
          tags$h4(paste("Domain annotations for Gene ID:", selected_gene_id)),
          DTOutput("domainTable")
        )
      } else {
        tags$p("No domain annotations found for the selected gene.")
      }
    } else {
      NULL
    }
  })
  
  
  # Disconnect from database when session ends
  session$onSessionEnded(function() {
    dbDisconnect(SFcon)
  })
}

shinyApp(ui = ui, server = server)