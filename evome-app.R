# TODO: Create github repo for this, README, and associated files

library(shiny)
library(tidyverse)
library(ggrepel)

wd <- "~/evome-app/"
setwd(wd)

path <- "data/EVome_TPMs.csv"
df_evome <- read_csv(path)

# Removing troublesome characters
names(df_evome) <-
  gsub(x = names(df_evome),
       pattern = "\\(|\\)",
       replacement = "")

organ_vec <- colnames(df_evome)[4:10]
tissue_vec <- colnames(df_evome)[11:36]
neuron_vec <- colnames(df_evome)[37:50]

server <- function(input, output, session) {
  
  nx <- reactive(
    df_evome %>%
      select(
        input$nx1,
        input$nx2,
        input$nx3
      ) %>%
      rowSums()
  )
  
  dx <- reactive(
    df_evome %>%
      select(
        input$dx1,
        input$dx2,
        input$dx3
      ) %>%
      rowSums()
  )
  
  ny <- reactive(
    df_evome %>%
      select(
        input$ny1,
        input$ny2,
        input$ny3
      ) %>%
      rowSums()
  )
  
  dy <- reactive(
    df_evome %>%
      select(
        input$dy1,
        input$dy2,
        input$dy3
      ) %>%
      rowSums()
  )
  
  # TODO: pivot_longer and filter features
  nxdx_nydy_df <- reactive(
    df_evome %>%
      select(
        Peptide_count,
        Gene_name
      ) %>%
      mutate(., nx=nx(), dx=dx(), ny=ny(), dy=dy()) %>%
      mutate(., nx_dx=nx/dx, ny_dy=ny/dy) %>%
      filter(is.finite(nx_dx) & is.finite(ny_dy))
  )
  
  observe({
    print(input$nx1)
    print(nxdx_nydy_df())
  })
  
  # Plot code modified from Inna
  output$nxdx_nydy_plot <- renderPlot({
    ggplot(nxdx_nydy_df(), aes(x = nx_dx, y = ny_dy), alpha = 0.5) +
      geom_point(color = "grey",
                 shape = 17,
                 size = 0.8) +
      # Covered by both thresholds
      geom_point(
        data = nxdx_nydy_df() %>% filter((nx_dx > input$x_thresh) &
                                (ny_dy > input$y_thresh)),
        aes(x = nx_dx, y = ny_dy, size = Peptide_count),
        color = "#601A4A",
        alpha = 0.5
      ) +
      geom_text_repel(data = nxdx_nydy_df() %>% filter((nx_dx > input$x_thresh) &
                                              (ny_dy > input$y_thresh)),
                      aes(x = nx_dx, y = ny_dy, label = Gene_name),
                      size = 2) +
      # Covered by x threshold
      geom_point(
        data = nxdx_nydy_df() %>% filter((nx_dx > input$x_thresh) &
                                (ny_dy < input$y_thresh)),
        aes(x = nx_dx, y = ny_dy, size = Peptide_count),
        color = "#EE442F",
        alpha = 0.5
      ) +
      geom_text_repel(data = nxdx_nydy_df() %>% filter((nx_dx > 20) &
                                            (ny_dy < input$y_thresh)),
                      aes(x = nx_dx, y = ny_dy, label = Gene_name),
                      size = 2) +
      # Covered by y threshold
      geom_point(
        data = nxdx_nydy_df() %>% filter((nx_dx < input$x_thresh) &
                                (ny_dy > input$y_thresh)),
        aes(x = nx_dx, y = ny_dy, size = Peptide_count),
        color = "green",
        alpha = 0.5
      ) +
      geom_text_repel(data = nxdx_nydy_df() %>% filter((nx_dx < input$x_thresh) &
                                              (ny_dy > input$y_thresh)),
                      aes(x = nx_dx, y = ny_dy, label = Gene_name),
                      size = 2) +
      theme_minimal() +
      ylim(0, 50) + # Replace 50 with input$ymax
      xlim(0, 50) + # Replace 50 with input$xmax
      ylab(paste("Enrichment in x cases over x controls")) +
      xlab(paste("Enrichment in y cases over y controls")) +
      theme(
        axis.ticks.y.left = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_text(
          size = 10,
          color = "black",
          margin = margin(
            t = 0,
            r = 10,
            b = 0,
            l = 0
          )
        ),
        axis.title.x = element_text(
          size = 10,
          color = "black",
          margin = margin(
            t = 10,
            r = 0,
            b = 0,
            l = 0
          )
        ),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10,
                                   color = "black")
      )
  })
  
  output$instructions_text <- renderText({
    "TODO"
  })
  
  output$nx_ui <- renderUI({
    fluidRow(
      column(4,
             checkboxGroupInput(
               inputId="nx1",
               label="Organs",
               choices=organ_vec,
               selected=list(organ_vec[1])
             ),
        ),
      column(4,
             checkboxGroupInput(
               inputId="nx2",
               label="Tissues",
               choices=tissue_vec,
               selected=NULL
             ),
      ),
      column(4,
             checkboxGroupInput(
               inputId="nx3",
               label="Neurons",
               choices=neuron_vec,
               selected=NULL
             ),
      ),
    )
  })
  
  output$dx_ui <- renderUI({
    fluidRow(
      column(4,
             checkboxGroupInput(
               inputId="dx1",
               label="Organs",
               choices=organ_vec,
               selected=organ_vec[2]
             ),
      ),
      column(4,
             checkboxGroupInput(
               inputId="dx2",
               label="Tissues",
               choices=tissue_vec,
               selected=NULL
             ),
      ),
      column(4,
             checkboxGroupInput(
               inputId="dx3",
               label="Neurons",
               choices=neuron_vec,
               selected=NULL
             ),
      ),
    )  
  })
  
  output$ny_ui <- renderUI({
    fluidRow(
      column(4,
             checkboxGroupInput(
               inputId="ny1",
               label="Organs",
               choices=organ_vec,
               selected=organ_vec[3]
             ),
      ),
      column(4,
             checkboxGroupInput(
               inputId="ny2",
               label="Tissues",
               choices=tissue_vec,
               selected=NULL
             ),
      ),
      column(4,
             checkboxGroupInput(
               inputId="ny3",
               label="Neurons",
               choices=neuron_vec,
               selected=NULL
             ),
      ),
    )
  })
  
  output$dy_ui <- renderUI({
    fluidRow(
      column(4,
             checkboxGroupInput(
               inputId="dy1",
               label="Organs",
               choices=organ_vec,
               selected=organ_vec[4]
             ),
      ),
      column(4,
             checkboxGroupInput(
               inputId="dy2",
               label="Tissues",
               choices=tissue_vec,
               selected=NULL
             ),
      ),
      column(4,
             checkboxGroupInput(
               inputId="dy3",
               label="Neurons",
               choices=neuron_vec,
               selected=NULL
             ),
      ),
    )
  })
  
  output$thresh_ui <- renderUI({
    tagList(
      sliderInput(inputId="x_thresh", label="x threshold", min=0, max=50, value=20),
      sliderInput(inputId="y_thresh", label="y threshold", min=0, max=50, value=20)
    )
  })
  
  output$plot_ui <- renderUI({
    tagList(
      textInput(inputId="xmax", label="Maximum x value", value=50),
      textInput(inputId="ymax", label="Maximum y value", value=50)
    )
  })
}

ui <- fluidPage(# App title
  headerPanel("EVome visualization (WORK IN PROGRESS)"),
  fluidPage(
    tabsetPanel(
      tabPanel(
        title="Plot",
        plotOutput("nxdx_nydy_plot"),
        uiOutput("plot_ui")
      ),
      tabPanel(
        title="Table",
        tableOutput("nxdx_nydy_df")
      ),
      tabPanel(
        title="x-axis Case",
        uiOutput("nx_ui")
      ),
      tabPanel(
        title="x-axis Control",
        uiOutput("dx_ui")
      ),
      tabPanel(
        title="y-axis Case",
        uiOutput("ny_ui")
      ),
      tabPanel(
        title="y-axis Control",
        uiOutput("dy_ui")
      ),
      tabPanel(
        title="Thresholds",
        uiOutput("thresh_ui")
      ),
      tabPanel(
        title="Help",
        textOutput("instructions_text")
      )
    )
  )
    
  )

# I cast spell, ~run shiny app~
shinyApp(ui, server)
