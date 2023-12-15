library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(ggplot2)
library(VariantAnnotation)
library(Gviz)
library(GenomicRanges)
library(DelayedMatrixStats)
library(sparseMatrixStats)

ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "VCF Explorer"),
  dashboardSidebar(),
  dashboardBody(
    fluidRow(
      box(
        fileInput("vcfFile", "Choose a VCF File"),
        selectInput("chromosomeInput", "Select Chromosome", choices = NULL,
                    width = "50%"),
        actionButton("submitBtn", "Submit"),
        div(
          style = "text-align: left; border: 1px solid #ccc; padding: 10px; border-radius: 5px; margin-top: 10px;",
          h4("Samples", style = "font-weight: bold; font-size: 14px;"),
          textOutput("sampleNames")
        )
      ),
      box(
        plotOutput("variantPlot", height = 300)
      ),
      box(
        plotOutput("plot3", height = 300), 
        width = 12
      )
    )
  )
)

server <- function(input, output, session) {
  # Create a reactive expression
  options(shiny.maxRequestSize=50*1024^2)
  vcfFile <- reactive({
    req(input$vcfFile)
    readVcf(input$vcfFile$datapath, "hg19")
  })
  
  observe({
    # Update the choices for the chromosome selection based on the VCF file
    if (!is.null(vcfFile())) {
      chromosomes <- unique(seqnames(vcfFile()))
      updateSelectInput(session, "chromosomeInput", choices = chromosomes)
    }
  })
  
  observeEvent(input$submitBtn, {
    # Count SNVs and Indels
    snvCount <- sum(isSNV(vcfFile()))
    indelCount <- sum(isIndel(vcfFile()))
    
    # Create a data frame for plotting
    variantCounts <- data.frame(
      Variant = c("SNV", "Indel"),
      Count = c(snvCount, indelCount)
    )
    
    # Plotting with updated colors
    output$variantPlot <- renderPlot({
      ggplot(variantCounts, aes(x = Variant, y = Count, fill = Variant)) +
        geom_bar(stat = "identity") +
        labs(title = "Variant Counts", x = "Variant Type", y = "Count") +
        theme_minimal() +
        scale_fill_manual(values = c("#06A3DA", "#34AD54"))
    })
    
    # Update the itrack based on the selected chromosome
    selectedChromosome <- input$chromosomeInput
    itrack <- IdeogramTrack(genome = "hg19", chromosome = selectedChromosome)
    vTrack <- AnnotationTrack(rowRanges(vcfFile()))
    
    output$plot3 <- renderPlot({
      plotTracks(list(itrack, vTrack))
    })
    
    # Display sample names below the select chromosome option
    output$sampleNames <- renderText({
      paste(samples(header(vcfFile())), collapse = ", ")
    })
  })
}

shinyApp(ui, server)
