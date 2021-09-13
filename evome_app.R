# TODO: Create github repo for this, README, and associated files

# TODO:
# Rewrite docx help page as html (wait for Alex edits)
# Use includeHTML to display help page
# Change x and y range maxes to dynamic
# Include inf values in the table instead of dropping them
# Label genes above threshold in pt1 and pt2
# 

library(shiny)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggpointdensity)
library(ggpubr)

path <- "data/EVome_TPMs.csv"
df_evome_raw <- read_csv(path)
path2 <- "data/EVome_TPMs_part2.csv"
df_part2 <- read_csv(path2)
df_evome_raw <- full_join(df_evome_raw, df_part2)
male_genes_path <- "data/MaleSpecificCargoCandidates.csv"
df_male_genes <- read_csv(male_genes_path)

# Removing troublesome characters
names(df_evome_raw) <-
  gsub(x = names(df_evome_raw),
       pattern = "\\(|\\)",
       replacement = "")

organ_vec <- colnames(df_evome_raw)[4:10]
tissue_vec <- colnames(df_evome_raw)[11:36]
neuron_vec <- colnames(df_evome_raw)[37:length(colnames(df_evome_raw))]

server <- function(input, output, session) {
  
  df_evome <- reactive({
    if (input$exclude_male_genes == "yes") {
      df_evome_raw %>% 
        filter(., !(Gene_name %in% df_male_genes$Gene_name))
    } else {
      df_evome_raw
    }
  })
  
  nx_names <- reactive(c(
    input$nx1,
    input$nx2,
    input$nx3
  ))
  nx <- reactive(
    df_evome() %>%
      select(
        input$nx1,
        input$nx2,
        input$nx3
      ) %>%
      rowSums() / sum(length(input$nx1), length(input$nx2), length(input$nx3))
  )
  
  dx_names <- reactive(c(
    input$dx1,
    input$dx2,
    input$dx3
  ))
  dx <- reactive(
    df_evome() %>%
      select(
        input$dx1,
        input$dx2,
        input$dx3
      ) %>%
      rowSums() / sum(length(input$dx1), length(input$dx2), length(input$dx3))
  )
  
  ny_names <- reactive(c(
    input$ny1,
    input$ny2,
    input$ny3
  ))
  ny <- reactive(
    df_evome() %>%
      select(
        input$ny1,
        input$ny2,
        input$ny3
      ) %>%
      rowSums() / sum(length(input$ny1), length(input$ny2), length(input$ny3))
  )
  
  dy_names <- reactive(c(
    input$dy1,
    input$dy2,
    input$dy3
  ))
  dy <- reactive(
    df_evome() %>%
      select(
        input$dy1,
        input$dy2,
        input$dy3
      ) %>%
      rowSums() / sum(length(input$dy1), length(input$dy2), length(input$dy3))
  )
  
  # TODO: pivot_longer and filter features
  # TODO: note excluded rows from is.finite filter
  
  nxdx_nydy_df_table <- reactive({
      df_evome() %>%
        select(Peptide_count,Gene_name) %>%
        mutate(Sample_x=nx(), Background_x=dx(), Sample_y=ny(), Background_y=dy(),
               Sample_over_background_x=Sample_x/Background_x, Sample_over_background_y=Sample_y/Background_y,
               color = case_when(
                 Sample_over_background_x > input$x_thresh & Sample_over_background_y < input$y_thresh~ "onlyX",
                 Sample_over_background_y > input$y_thresh & Sample_over_background_x < input$x_thresh~ "onlyY",
                 (Sample_over_background_x > input$x_thresh & Sample_over_background_y > input$y_thresh) ~ "bothXY",
                 (Sample_over_background_x < input$x_thresh & Sample_over_background_y < input$y_thresh) ~ "none")) %>%
        rename(
          Sample_x = Sample_x,
          Sample_y = Sample_y,
          Background_x = Background_x,
          Background_y = Background_y,
          Sample_over_background_x = Sample_over_background_x,
          Sample_over_background_y = Sample_over_background_y
        )
  })
  
nxdx_nydy_df <- reactive({ 
  nxdx_nydy_df_table() %>%
    filter(is.finite(Sample_over_background_x) & is.finite(Sample_over_background_y))
  })
  
  pt1 <- reactive({
    ggplot(nxdx_nydy_df(), aes(x = Sample_x, y = Background_x)) +
      geom_pointdensity(size = 1) +
      # ggtitle("Sample Y vs. Background Y") +
      xlab("Avg. peptide abundance Background Y") +
      ylab("Avg. peptide abundance Sample Y") +
      theme_pubr() +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      # coord_equal() +
      # xlim(range(nxdx_nydy_df() %>% select(Sample_x,Background_x) %>% unlist())) +
      # ylim(range(nxdx_nydy_df() %>% select(Sample_x,Background_x) %>% unlist())) +
      scale_color_viridis_c(option = "B",
                            name = "Neighboring\npoints",
                            guide = "none")+
      theme(plot.title = element_text(hjust = 0.5),
            aspect.ratio = 1) +
      stat_cor(method = "spearman")
  })
  
  pt2 <-reactive({
    ggplot(nxdx_nydy_df(), aes(x = Sample_y, y = Background_y)) +
      geom_pointdensity(size = 1) +
      # ggtitle("Sample X vs. Background X") +
      xlab("Avg. peptide abundance Background X") +
      ylab("Avg. peptide abundance Sample X") +
      theme_pubr() +
      # xlim(range(nxdx_nydy_df() %>% select(Sample_y,Background_y) %>% unlist())) +
      # ylim(range(nxdx_nydy_df() %>% select(Sample_y,Background_y) %>% unlist())) +
      # # coord_equal() +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_color_viridis_c(option = "B",
                            name = "Neighboring\npoints",
                            guide = "none")+
      theme(plot.title = element_text(hjust = 0.5),
            aspect.ratio = 1) +
      stat_cor(method = "spearman")
  })
  
  # output$axes_analysis_plots <- renderPlot({pt1() + pt2()})
  output$pt1 <- renderPlot(pt1()) 
  output$pt2 <- renderPlot(pt2()) 
  
  output$rendered_table <- renderDataTable(nxdx_nydy_df_table())
  
  output$render_nx <- renderUI(HTML(paste(c("Sample Y", nx_names()), sep = '<br/>')))
  output$render_dx <- renderUI(HTML(paste(c("Background Y", dx_names()), sep = '<br/>')))
  output$render_ny <- renderUI(HTML(paste(c("Sample X", ny_names()), sep = '<br/>')))
  output$render_dy <- renderUI(HTML(paste(c("Background X", dy_names()), sep = '<br/>')))
  
  observe({
    print(dy_names())
    print(nxdx_nydy_df())
  })
  
  # Select/deslect all buttons 
  observe({
    # x-axis case organs
    if (input$nx1_select > 0) {
      if (input$nx1_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="nx1",
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="nx1",
          choices=organ_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # x-axis case tissues
    if (input$nx2_select > 0) {
      if (input$nx2_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="nx2",
          choices=tissue_vec,
          selected=tissue_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="nx2",
          choices=tissue_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # x-axis case neurons
    if (input$nx3_select > 0) {
      if (input$nx3_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="nx3",
          choices=neuron_vec,
          selected=neuron_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="nx3",
          choices=neuron_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # x-axis control organs
    if (input$dx1_select > 0) {
      if (input$dx1_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="dx1",
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dx1",
          choices=organ_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # x-axis control tissues
    if (input$dx2_select > 0) {
      if (input$dx2_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="dx2",
          choices=tissue_vec,
          selected=tissue_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dx2",
          choices=tissue_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # x-axis control neurons
    if (input$dx3_select > 0) {
      if (input$dx3_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="dx3",
          choices=neuron_vec,
          selected=neuron_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dx3",
          choices=neuron_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # y-axis case organs
    if (input$ny1_select > 0) {
      if (input$ny1_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="ny1",
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="ny1",
          choices=organ_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # y-axis case tissues
    if (input$ny2_select > 0) {
      if (input$ny2_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="ny2",
          choices=tissue_vec,
          selected=tissue_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="ny2",
          choices=tissue_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # y-axis case neurons
    if (input$ny3_select > 0) {
      if (input$ny3_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="ny3",
          choices=neuron_vec,
          selected=neuron_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="ny3",
          choices=neuron_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # y-axis control organs
    if (input$dy1_select > 0) {
      if (input$dy1_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="dy1",
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dy1",
          choices=organ_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # y-axis control tissues
    if (input$dy2_select > 0) {
      if (input$dy2_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="dy2",
          choices=tissue_vec,
          selected=tissue_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dy2",
          choices=tissue_vec,
          selected=NULL
        )
      }
    }
  })
  observe({
    # y-axis control neurons
    if (input$dy3_select > 0) {
      if (input$dy3_select %% 2 == 0) {
        updateCheckboxGroupInput(
          session=session,
          inputId="dy3",
          choices=neuron_vec,
          selected=neuron_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dy3",
          choices=neuron_vec,
          selected=NULL
        )
      }
    }
  })
  
  output$downloadPlot <- downloadHandler(
    filename = "evome_plot.png",
    content = function(file) {
      ggsave(file, plot_reactive())
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = "evome_data.csv",
    content = function(file) {
      write.csv(nxdx_nydy_df(), file)
    }
  )
  
  plot_reactive <- reactive(
    ggplot(nxdx_nydy_df(), aes(x = Sample_over_background_x, y = Sample_over_background_y, color = color)) +
      geom_point(alpha = 0.5, size = 0.8) +
      # # Covered by both thresholds
      # geom_point(
      #   data = nxdx_nydy_df() %>% filter((Sample_over_background_x > input$x_thresh) &
      #                                      (Sample_over_background_y > input$y_thresh)),
      #   aes(x = Sample_over_background_x, y = Sample_over_background_y, size = 1),
      #   color = "#601A4A",
      #   alpha = 0.5
      # ) +
      # geom_text_repel(data = nxdx_nydy_df() %>% filter((Sample_over_background_x > input$x_thresh) &
      #                                                    (Sample_over_background_y > input$y_thresh)),
      #                 aes(x = Sample_over_background_x, y = Sample_over_background_y, label = Gene_name),
      #                 size = 5) +
      # # Covered by x threshold
      # geom_point(
      #   data = nxdx_nydy_df() %>% filter((Sample_over_background_x > input$x_thresh) &
      #                                      (Sample_over_background_y < input$y_thresh)),
      #   aes(x = Sample_over_background_x, y = Sample_over_background_y, size = 1),
      #   color = "#EE442F",
      #   alpha = 0.5
      # ) +
      # geom_text_repel(data = nxdx_nydy_df() %>% filter((Sample_over_background_x > input$x_thresh) &
      #                                                    (Sample_over_background_y < input$y_thresh)),
      #                 aes(x = Sample_over_background_x, y = Sample_over_background_y, label = Gene_name),
      #                 size = 5) +
      # # Covered by y threshold
      # geom_point(
      #   data = nxdx_nydy_df() %>% filter((Sample_over_background_x < input$x_thresh) &
      #                                      (Sample_over_background_y > input$y_thresh)),
      #   aes(x = Sample_over_background_x, y = Sample_over_background_y, size = 1),
      #   color = "green",
      #   alpha = 0.5
      # ) +
    scale_color_brewer(palette = "Set1") +
      guides(color = "none") +
      geom_text_repel(data = nxdx_nydy_df() %>% filter(color != "none"),
                      aes(x = Sample_over_background_x, y = Sample_over_background_y, label = Gene_name),
                      size = 5) +
      theme_pubr() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(input$yrange) +
      xlim(input$xrange) +
      # scale_x_log10() +
      # scale_y_log10() +
      # ggtitle("Enrichment in case/Background Y vs. enrichment in case/Background X") +
      ylab(paste("Enrichment in case/Background X")) +
      xlab(paste("Enrichment in case/Background Y")) +
      coord_fixed() +
      theme(aspect.ratio = 1)
      # theme(
      #   axis.ticks.y.left = element_line(color = "black"),
      #   axis.ticks.x = element_line(color = "black"),
      #   axis.line.x = element_line(color = "black"),
      #   axis.line.y = element_line(color = "black"),
      #   axis.title.y = element_text(
      #     size = 10,
      #     color = "black",
      #     margin = margin(
      #       t = 0,
      #       r = 10,
      #       b = 0,
      #       l = 0
      #     )
      #   ),
      #   axis.title.x = element_text(
      #     size = 10,
      #     color = "black",
      #     margin = margin(
      #       t = 10,
      #       r = 0,
      #       b = 0,
      #       l = 0
      #     )
      #   ),
        # axis.text.y = element_text(size = 10, color = "black"),
        # axis.text.x = element_text(size = 10,
        #                            color = "black")
      )
  # )
  
  # Plot code modified from Inna
  output$nxdx_nydy_plot <- renderPlot({
    plot_reactive()
  },
  #ignoreNULL=FALSE
  )
  
  
  # output$instructions_text <- renderText({
  #   "    1a. In the Sample Y tab, select which data you want.\n
  #   1b. In the Background Y tab, select which data you want.\n
  #   1c. Check the main plot tab to see a average peptide abundances for Sample Y vs. Background Y.\n
  #   2a. In the Sample X tab, select which data you want.\n
  #   2b. In the Background X tab, select which data you want.\n
  #   2c. Check the main plot tab to see a average peptide abundances for Sample X vs. Background X.\n
  #   3a. Check the main plot tab to see enrichment in case/Background Y vs. enrichment in case/Background X.\n
  #   3b. If desired, adjust the thresholds for plot 3a with sliders beneath the display.\n
  #   3c. If desired, adjust the x and y ranges for plot 3a with sliders beneath the display.\n
  #   3d. If desired, exclude male-specific genes with a drop-down menu beneath the display.\n
  #   3e. Press 'download plot' to save plot 3a to your machine.\n
  #   4a. Check the table tab to see the data table producing plot 3a.\n
  #   4b. Press 'download table' to save this data table to your machine."
  # })
  
}

ui <- fluidPage(# App title
  headerPanel("MyEVome"),
  fluidPage(
    tabsetPanel(
      tabPanel(
        title="Main plot",
        # fluidRow(
          # column(6, align = "center", plotOutput('nxdx_nydy_plot'))
        # ),
#        plotOutput("axes_analysis_plots"),
       plotOutput("nxdx_nydy_plot"),
        fluidRow(
          column(3,
                 sliderInput(inputId="xrange", label="X range", max = 50, min = 0, value = c(0, 45))
          ),
          column(3,
                 sliderInput(inputId="yrange", label="Y range", max=50, min = 0, value = c(0, 45))
          ),
          column(3,
                 sliderInput(inputId="x_thresh", label="X threshold", min=0, max=50, value=5)
          ),
          column(3,
                 sliderInput(inputId="y_thresh", label="Y threshold", min=0, max=50, value=5)
          )),
        fluidRow(
        column(4, selectInput(inputId = "exclude_male_genes", 
                    label = "Exclude male-specific genes?",
                    choices = c("no", "yes"),
                    selected = "no")),
        column(4,
               headerPanel(""),
               downloadButton('downloadPlot', 'Download Plot'))
        ),
        fluidRow(
          column(3,
                 "Sample Y selection:",
                 wellPanel(htmlOutput("render_nx"))),
          column(3,
                 "Background Y selection:",
                 wellPanel(htmlOutput("render_dx"))),
          column(3,
                 "Sample X selection:",
                 wellPanel(htmlOutput("render_ny"))),
          column(3,
                 "Background X selection:",
                 wellPanel(htmlOutput("render_dy")))
        ),
        fluidRow(
          column(6, plotOutput('pt1')),
          column(6, plotOutput('pt2'))
        )
      ),
      tabPanel(
        title="Table",
        downloadButton('downloadData', 'Download Data'),
        dataTableOutput("rendered_table")
      ),
      tabPanel(
        title="Sample Y",
        fluidRow(
          column(4,
                 headerPanel(""),
                 actionButton("nx1_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="nx1",
                   label="Organs",
                   choices=organ_vec,
                   selected=NULL
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("nx2_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="nx2",
                   label="Tissues",
                   choices=tissue_vec,
                   selected=NULL
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("nx3_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="nx3",
                   label="Neurons",
                   choices=neuron_vec,
                   selected=c("Cholinergic_15")
                 ),
          ),
        ),
      ),
      tabPanel(
        title="Background Y",
        fluidRow(
          column(4,
                 headerPanel(""),
                 actionButton("dx1_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="dx1",
                   label="Organs",
                   choices=organ_vec,
                   selected=NULL
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("dx2_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="dx2",
                   label="Tissues",
                   choices=tissue_vec,
                   selected=NULL
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("dx3_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="dx3",
                   label="Neurons",
                   choices=neuron_vec,
                   selected=neuron_vec
                 ),
          ),
        ),  
      ),
      tabPanel(
        title="Sample X",
        fluidRow(
          column(4,
                 headerPanel(""),
                 actionButton("ny1_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="ny1",
                   label="Organs",
                   choices=organ_vec,
                   selected=NULL
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("ny2_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="ny2",
                   label="Tissues",
                   choices=tissue_vec,
                   selected=c("Ciliated_sensory_neurons")
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("ny3_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="ny3",
                   label="Neurons",
                   choices=neuron_vec,
                   selected=NULL
                 ),
          ),
        ),
      ),
      tabPanel(
        title="Background X",
        fluidRow(
          column(4,
                 headerPanel(""),
                 actionButton("dy1_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="dy1",
                   label="Organs",
                   choices=organ_vec,
                   selected=c("Intestine")
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("dy2_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="dy2",
                   label="Tissues",
                   choices=tissue_vec,
                   selected=tissue_vec
                 ),
          ),
          column(4,
                 headerPanel(""),
                 actionButton("dy3_select", label="Select/Deselect all"),
                 checkboxGroupInput(
                   inputId="dy3",
                   label="Neurons",
                   choices=neuron_vec,
                   selected=NULL
                 ),
          ),
        ),
      ),
      # tabPanel(
      #   title="Thresholds",
      #   tagList(
      #     sliderInput(inputId="x_thresh", label="x threshold", min=0, max=50, value=5),
      #     sliderInput(inputId="y_thresh", label="y threshold", min=0, max=50, value=5)
      #   )
      # ),
      tabPanel(
        title="Help",
        includeMarkdown("data/help.md")
      )
    )
  )
    
  )

# I cast spell, ~run shiny app~
shinyApp(ui, server)
