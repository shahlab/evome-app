# TODO: Create github repo for this, README, and associated files

# TODO:
# 1. Download button for plot/table (DONE)
# 2. Set default plot to Inna's paper plot (DONE)
# 3. List which points are excluded
# 4. Error message: "we recommend including Case in control"
# 5. Sortable table (DONE)
# 6. Fill out instructions tab
# 7. Replace old data frame with new one (DONE)
# 8. Select all option (DONE)
# 9. Layout, font sizes, descriptive nomenclature
# 10. Add two plots above enrichment, Case vs. control x, Case vs. control y (DONE)
# 11. Presets
# 12. Button to exclude male-specific vector (DONE)
# 13. Find out how to integrate with paper
# 14. Case A / case B (DONE)
# 15. Equal coordinate system
# 16. Buttons on top
# 17. Remove resizing of points by peptide count
# 18. Add dynamic point highlighting on axes_analyses_plots
# 19. Make plot download with white background (instead of transparent)
# 20. Increase downloaded plot size
# 

library(shiny)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggpointdensity)

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
      rowSums()
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
      rowSums()
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
  
  nxdx_nydy_df <- reactive({
      df_evome() %>%
        select(Peptide_count,Gene_name) %>%
        mutate(Case_x=nx(), Control_x=dx(), Case_y=ny(), Control_y=dy(),
               Case_over_control_x=Case_x/Control_x, Case_over_control_y=Case_y/Control_y,
               color = case_when(
                 Case_over_control_x > input$x_thresh & Case_over_control_y < input$y_thresh~ "onlyX",
                 Case_over_control_y > input$y_thresh & Case_over_control_x < input$x_thresh~ "onlyY",
                 (Case_over_control_x > input$x_thresh & Case_over_control_y > input$y_thresh) ~ "bothXY",
                 (Case_over_control_x < input$x_thresh & Case_over_control_y < input$y_thresh) ~ "none")) %>%
        filter(is.finite(Case_over_control_x) & is.finite(Case_over_control_y))
  })
  
  pt1 <- reactive({
    ggplot(nxdx_nydy_df(), aes(x = Case_x, y = Control_x)) +
      geom_pointdensity(size = 1) +
      ggtitle("Case A vs. control A") +
      xlab("Control A") +
      ylab("Case A") +
      theme_minimal() +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      coord_equal() +
      scale_color_viridis_c(option = "B",
                            name = "Neighboring\npoints",
                            guide = "none")+
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  pt2 <-reactive({
    ggplot(nxdx_nydy_df(), aes(x = Case_y, y = Control_y)) +
      geom_pointdensity(size = 1) +
      ggtitle("Case B vs. control B") +
      xlab("Control B") +
      ylab("Case B") +
      theme_minimal() +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      coord_equal() +
      scale_color_viridis_c(option = "B",
                            name = "Neighboring\npoints",
                            guide = "none")+
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$axes_analysis_plots <- renderPlot({pt1() + pt2()})
  
  output$rendered_table <- renderDataTable(nxdx_nydy_df())
  
  output$render_nx <- renderUI(HTML(paste(nx_names(), sep = '<br/>')))
  output$render_dx <- renderUI(HTML(paste(dx_names(), sep = '<br/>')))
  output$render_ny <- renderUI(HTML(paste(ny_names(), sep = '<br/>')))
  output$render_dy <- renderUI(HTML(paste(dy_names(), sep = '<br/>')))
  
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
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="nx3",
          choices=organ_vec,
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
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dx3",
          choices=organ_vec,
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
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="ny3",
          choices=organ_vec,
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
          choices=organ_vec,
          selected=organ_vec
        )
      }
      else {
        updateCheckboxGroupInput(
          session=session,
          inputId="dy3",
          choices=organ_vec,
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
    ggplot(nxdx_nydy_df(), aes(x = Case_over_control_x, y = Case_over_control_y, color = color)) +
      geom_point(alpha = 0.5, size = 0.8) +
      # # Covered by both thresholds
      # geom_point(
      #   data = nxdx_nydy_df() %>% filter((Case_over_control_x > input$x_thresh) &
      #                                      (Case_over_control_y > input$y_thresh)),
      #   aes(x = Case_over_control_x, y = Case_over_control_y, size = 1),
      #   color = "#601A4A",
      #   alpha = 0.5
      # ) +
      # geom_text_repel(data = nxdx_nydy_df() %>% filter((Case_over_control_x > input$x_thresh) &
      #                                                    (Case_over_control_y > input$y_thresh)),
      #                 aes(x = Case_over_control_x, y = Case_over_control_y, label = Gene_name),
      #                 size = 5) +
      # # Covered by x threshold
      # geom_point(
      #   data = nxdx_nydy_df() %>% filter((Case_over_control_x > input$x_thresh) &
      #                                      (Case_over_control_y < input$y_thresh)),
      #   aes(x = Case_over_control_x, y = Case_over_control_y, size = 1),
      #   color = "#EE442F",
      #   alpha = 0.5
      # ) +
      # geom_text_repel(data = nxdx_nydy_df() %>% filter((Case_over_control_x > input$x_thresh) &
      #                                                    (Case_over_control_y < input$y_thresh)),
      #                 aes(x = Case_over_control_x, y = Case_over_control_y, label = Gene_name),
      #                 size = 5) +
      # # Covered by y threshold
      # geom_point(
      #   data = nxdx_nydy_df() %>% filter((Case_over_control_x < input$x_thresh) &
      #                                      (Case_over_control_y > input$y_thresh)),
      #   aes(x = Case_over_control_x, y = Case_over_control_y, size = 1),
      #   color = "green",
      #   alpha = 0.5
      # ) +
    scale_color_brewer(palette = "Set1") +
      guides(color = "none") +
      geom_text_repel(data = nxdx_nydy_df() %>% filter(color != "none"),
                      aes(x = Case_over_control_x, y = Case_over_control_y, label = Gene_name),
                      size = 5) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(input$yrange) +
      xlim(input$xrange) +
      # scale_x_log10() +
      # scale_y_log10() +
      ggtitle("Enrichment in case/control A vs. enrichment in case/control B") +
      ylab(paste("Enrichment in case/control B")) +
      xlab(paste("Enrichment in case/control A")) +
      coord_equal()
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
  
  output$instructions_text <- renderText({
    "TODO"
  })
  
  # TODO: set default plot to the one from Inna's paper
  
}

ui <- fluidPage(# App title
  headerPanel("MyEVome"),
  fluidPage(
    tabsetPanel(
      tabPanel(
        title="Main plot",
        plotOutput("axes_analysis_plots"),
        plotOutput("nxdx_nydy_plot"),
        fluidRow(
          column(3,
                 sliderInput(inputId="yrange", label="Y Range", max=50, min = 0, value = c(0, 45))
                 ),
          column(3,
                  sliderInput(inputId="xrange", label="X Range", max = 50, min = 0, value = c(0, 45))
                 ),
          column(3,
                 sliderInput(inputId="x_thresh", label="x threshold", min=0, max=50, value=5)
          ),
          column(3,
                 sliderInput(inputId="y_thresh", label="y threshold", min=0, max=50, value=5)
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
                 htmlOutput("render_nx")),
          column(3,
                 htmlOutput("render_dx")),
          column(3,
                 htmlOutput("render_ny")),
          column(3,
                 htmlOutput("render_dy"))
        )
      ),
      tabPanel(
        title="Table",
        downloadButton('downloadData', 'Download Data'),
        dataTableOutput("rendered_table")
      ),
      tabPanel(
        title="Case A",
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
        title="Control A",
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
        title="Case B",
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
        title="Control B",
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
        textOutput("instructions_text")
      )
    )
  )
    
  )

# I cast spell, ~run shiny app~
shinyApp(ui, server)
