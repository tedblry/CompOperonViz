
library(shiny)
library(ggplot2)
library(DT)
library(plotly)
library(dplyr)
library(zoo)
library(stringr)
library(egg)
library(tidyr)
library(shinyWidgets)
library(colourpicker)


# Read Data
log2getmm_deseq_full_path <- "./Annotated_Differential_abundance_result.csv"
log2getmm_deseq_full <- read.csv(log2getmm_deseq_full_path, header = TRUE, stringsAsFactors = FALSE, row.names = "locus_tag")

operon_data <- read.csv("./operon_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Data Processing
# Fill down the Operon numbers for each gene
operon_data$operon <- zoo::na.locf(operon_data$operon)

# Filter out rows without position data and select relevant columns
operon_data_processed <- operon_data %>%
  filter(!is.na(start)) %>%
  select(start, end, operon)

# Add operon_unit column
operon_data_full <- operon_data_processed %>%
  group_by(operon) %>%
  mutate(operon_unit = ifelse(n() > 1, cur_group_id(), NA)) %>%
  ungroup() %>%
  mutate(operon_unit = ifelse(is.na(operon_unit), NA, 
                              paste0("unit_", match(operon, unique(operon[!is.na(operon_unit)]))))) 

# Function to combine Preferred_name and PFAMs
combine_preferred_pfam <- function(preferred, pfams) {
  if (is.na(preferred)) {
    return(as.character(pfams))
  } else {
    return(paste0("[", preferred, "]  ", pfams))
  }
}

# Combine DESeq2 Data with Operon Data
operon_lfc <- log2getmm_deseq_full %>%
  left_join(operon_data_full, by="start") %>%  
  # counting operon units (that is more than one gene assigned the same operon number)
  group_by(operon) %>%
  mutate(operon_unit = ifelse(n() > 1, cur_group_id(), NA)) %>%
  ungroup() %>%
  mutate(operon_unit = ifelse(is.na(operon_unit), NA, 
                              paste0("unit_", match(operon, unique(operon[!is.na(operon_unit)]))))) %>% 
  select(log2FoldChange, operon_unit, start, end.x, gene_length, padj, COG_category, Preferred_name, PFAMs) %>%
  mutate(COG_category = replace_na(COG_category, "S")) %>% 
  mutate(COG_category = ifelse(nchar(COG_category) > 1, "Joint COG", COG_category)) %>%
  mutate(log2FoldChange = pmin(pmax(log2FoldChange, -2), 2)) %>%
  mutate(pval_sig = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ NA_character_
  )) %>%
  tibble::rownames_to_column() %>% 
  dplyr::rename("end" = "end.x",
         "rel_pos" = "rowname") %>%
  mutate(rel_pos = as.numeric(rel_pos) - 1) %>%
  mutate(COG_category = as.factor(COG_category)) %>%
  mutate(gene_pfam = mapply(combine_preferred_pfam, Preferred_name, PFAMs)) %>%
  mutate(gene_pfam = str_replace_all(gene_pfam, ",", ", ")) %>%
  mutate(gene_name = str_c("- ", Preferred_name)) %>% 
  mutate(gene_pfam = str_c("- ", gene_pfam)) %>% 
  select(-c(start, end, gene_length, Preferred_name, PFAMs))

# Print Summary
cat("There are ", length(unique(operon_data$operon)), " operon units with more than 2 genes out of ", length(unique(operon_lfc$operon_unit)), " identified operons.\n")

# Create the Visualization
# Define COG_Operon_color
COG_Operon_color <- c(
  "M" = "#890000",
  "K" = "#CC79A7",
  "L" = "#673AB7",
  "H" = "#000079",
  "E" = "#008080",
  "G" = "#006000",
  "C" = "#76FF03", 
  "T" = "#FFEB3B",
  "P" = "#FF7000",
  "U" = "#F15C80",
  "J" = "#9C27B0", 
  "I" = "#03A9F4", 
  "V" = "#657D9B",
  "O" = "#C5FA30",
  "F" = "#FFB000",
  "D" = "#D8BFD8", 
  "Q" = "#795548",
  "S" = "#9E9E9E",
  "Joint COG" = "#ECECEC",
  "Operon_A" = "#d55e00",
  "Operon_B" = "#1f449c"
)

# Function to assign colors to operon_unit
assign_operon_colors <- function(operon_unit) {
  unique_units <- na.omit(unique(operon_unit))
  unit_colors <- c("Operon_A", "Operon_B")
  colors <- rep(unit_colors, length.out = length(unique_units))
  names(colors) <- unique_units
  return(colors[operon_unit])
}

# Pre-process the data
operon_lfc$ymin_operon <- ifelse(is.na(operon_lfc$operon_unit), 0, -2.1)
operon_lfc$ymax_operon <- ifelse(is.na(operon_lfc$operon_unit), 0, 2.1)

# Add operon_color column to operon_lfc
operon_lfc$operon_color <- ifelse(is.na(operon_lfc$operon_unit), NA, assign_operon_colors(operon_lfc$operon_unit))

ui <- fluidPage(
  # Custom CSS to adjust margins and ensure full width for sliderInput
  tags$head(
    tags$style(HTML("
      .custom-container {
        margin-left: 5%;
        margin-right: 5%;
        margin-bottom: 1%;
      }
      .shiny-input-container {
        width: 100% !important;
      }
    "))
  ),
  
  # Wrap titlePanel in a div with custom class
  div(class = "custom-container",
      titlePanel("Operon Data Visualization")
  ),
  div(class = "custom-container",
      tags$footer(
        HTML("<p>Please cite the following when using this tool:<br> 
           Byeongyeon Cho, Grace Moore, Loc-Duyen Pham et al. DNA characterization reveals potential operon-unit packaging of extracellular vesicle cargo from a gut bacterial symbiont, 04 December 2023, PREPRINT (Version 1) available at <a href='https://doi.org/10.21203/rs.3.rs-3689023/v1'>Research Square: doi.org/10.21203/rs.3.rs-3689023/v1</a></p>")
      )
  ),
  
  # User selection forms in the first row spanning the entire width
  fluidRow(
    div(class = "custom-container",
        sliderInput("posRange",
                    "Select Range of Relative Position:",
                    min = 0,
                    max = nrow(operon_lfc),
                    value = c(0, 150),  # Set default range from 0 to 150
                    step = 1)
    )
  ),
  
  # Additional controls in the second row
  fluidRow(
    div(class = "custom-container",
        checkboxInput("hideXLabels", "Hide x-axis labels", FALSE),
        actionButton("viewButton", "View Comparative Microbial Operon Plot")
    )
  ),
  
  # File upload inputs
  fluidRow(
    div(class = "custom-container",
        column(6, fileInput("log2getmm_deseq_full_file", "Upload CSV file containing the eggnNOG/prokka annotated comparative differential abundance DNA counts", accept = c(".csv"))),
        column(6, fileInput("operon_data_file", "Upload CSV file containing the operon-mapper v2 generated operon_list information", accept = c(".csv")))
    )
  ),
  
  
  # Tabs for plot and data table in the third row
  fluidRow(
    div(class = "custom-container",
        tabsetPanel(
          tabPanel("Plot", 
                   fluidRow(
                     column(2, colourInput("M_color", "M: Cell wall structure and biogenesis and outer membrane", "#890000")),
                     column(2, colourInput("K_color", "K: Transcription", "#CC79A7")),
                     column(2, colourInput("L_color", "L: Replication, recombination and repair", "#673AB7")),
                     column(2, colourInput("H_color", "H: Coenzyme metabolism", "#000079")),
                     column(2, colourInput("E_color", "E: Amino acid metabolism and transport", "#008080")),
                     column(2, colourInput("G_color", "G: Carbohydrate metabolism and transport", "#006000"))
                   ),
                   fluidRow(
                     column(2, colourInput("C_color", "C: Energy production and conversion", "#76FF03")),
                     column(2, colourInput("T_color", "T: Signal transduction", "#FFEB3B")),
                     column(2, colourInput("P_color", "P: Inorganic ion transport and metabolism", "#FF7000")),
                     column(2, colourInput("U_color", "U: Lipid metabolism", "#F15C80")),
                     column(2, colourInput("J_color", "J: Translation, including ribosome structure and biogenesis", "#9C27B0")),
                     column(2, colourInput("I_color", "I: Lipid metabolism", "#03A9F4"))
                   ),
                   fluidRow(
                     column(2, colourInput("V_color", "V: Defense mechanism", "#657D9B")),
                     column(2, colourInput("O_color", "O: Molecular chaperones and related functions", "#C5FA30")),
                     column(2, colourInput("F_color", "F: Nucleotide metabolism and transport", "#FFB000")),
                     column(2, colourInput("D_color", "D: Cell division and chromosome partitioning", "#D8BFD8")),
                     column(2, colourInput("Q_color", "Q: Secondary metabolites biosynthesis, transport and catabolism", "#795548")),
                     column(2, colourInput("S_color", "S: No functional prediction", "#9E9E9E"))
                   ),
                   plotlyOutput("operonPlot")),
          tabPanel("Data Table", DTOutput("dataTable"))
        )
    )
  )
)




# Global variables to store initial ranges
initial_x_range <- NULL
initial_y_range <- NULL

# Define Server
server <- function(input, output, session) {
  # Reactive values to store the data
  log2getmm_deseq_full <- reactiveVal()
  operon_data <- reactiveVal()
  
  observe({
    # When a new file is uploaded for log2getmm_deseq_full
    file1 <- input$log2getmm_deseq_full_file
    if (is.null(file1)) return()
    
    # Read the file and update the reactive value
    log2getmm_deseq_full(read.csv(file1$datapath, header = TRUE, stringsAsFactors = FALSE, row.names = "locus_tag"))
  })
  
  observe({
    # When a new file is uploaded for operon_data
    file2 <- input$operon_data_file
    if (is.null(file2)) return()
    
    # Read the file and update the reactive value
    operon_data(read.csv(file2$datapath, header = TRUE, stringsAsFactors = FALSE))
  })
    
  # Reactive value to store the initial range selected by the user
  initial_x_range <- reactive({
    input$posRange
  })
  
  filteredData <- reactive({
    req(input$viewButton)
    operon_lfc %>%
      filter(rel_pos >= input$posRange[1] & rel_pos <= input$posRange[2])
  })

    COG_Operon_color <- reactive({
      c(
        "M" = input$M_color,
        "K" = input$K_color,
        "L" = input$L_color,
        "H" = input$H_color,
        "E" = input$E_color,
        "G" = input$G_color,
        "C" = input$C_color,
        "T" = input$T_color,
        "P" = input$P_color,
        "U" = input$U_color,
        "J" = input$J_color,
        "I" = input$I_color,
        "V" = input$V_color,
        "O" = input$O_color,
        "F" = input$F_color,
        "D" = input$D_color,
        "Q" = input$Q_color,
        "S" = input$S_color,
        "Joint COG" = "#ECECEC",
        "Operon_A" = "#d55e00",
        "Operon_B" = "#1f449c"
      )
    })
    
  # Use COG_Operon_color() in the createPlot function
  createPlot <- function() {
    data_to_plot <- filteredData()
    text_size_base <- 12
    range <- input$posRange[2] - input$posRange[1]
    text_size <- text_size_base / sqrt(range/100)
    
    p <- ggplot(data_to_plot) + 
    
      geom_rect(data = data_to_plot, aes(xmin = rel_pos, xmax = rel_pos + 1, ymin = ymin_operon, ymax = ymax_operon, fill = operon_color), alpha = 0.1, color = NA) +
      geom_rect(data = data_to_plot, aes(xmin = rel_pos, xmax = rel_pos + 1, ymin = -14, ymax = -2.1, fill = NA), alpha = 0, color = NA) +
      ggtitle("Operon-unit masked log2fold normalized gene counts change, ordered by relative genomic position") +  # Add title here
      geom_rect(aes(xmin = rel_pos + 0.2, xmax = rel_pos + 0.8, ymin = 0, ymax = log2FoldChange, fill = COG_category)) +
      scale_fill_manual(values = COG_Operon_color()) +  # Use COG_Operon_color() here
      theme_minimal() +
      theme_void()+
      egg::theme_article() +
      scale_y_continuous(breaks = c(-2.5, 0, 2.5))  # Add this line to modify y-axis ticks
      
    
    p_plotly <- ggplotly(p) %>% 
      layout(autosize = T, height = 1200,  # Set autosize to TRUE and remove the width parameter
             margin = list(l = 50, r = 50, b = 200, t = 100, pad = 4))  # Increase the bottom margin
    
    # Add annotations with specified text size
    for(i in 1:nrow(data_to_plot)) {
      if (!is.na(data_to_plot$pval_sig[i])) {
        p_plotly <- p_plotly %>% 
          add_annotations(
            text = data_to_plot$pval_sig[i], 
            x = data_to_plot$rel_pos[i] + 0.5, 
            y = 2.1, 
            showarrow = FALSE, 
            textangle = 90,
            xanchor = "center",
            yanchor = "bottom",
            font = list(size = text_size)
          )
      }
      
      if (!is.na(data_to_plot$gene_pfam[i]) && !input$hideXLabels) {
        p_plotly <- p_plotly %>% 
          add_annotations(
            text = data_to_plot$gene_pfam[i], 
            x = data_to_plot$rel_pos[i] + 0.5, 
            y = -2.1, 
            showarrow = FALSE, 
            textangle = 90,
            xanchor = "center",
            yanchor = "top",
            font = list(size = text_size)
          )
      }
    }
    
    return(p_plotly)
  }
  
  
  # Render initial plot
  output$operonPlot <- renderPlotly({
    createPlot()
  })
  
  output$dataTable <- renderDT({
    datatable(filteredData(), options = list(searching = TRUE, pageLength = 20))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
