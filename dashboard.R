library(shiny)
library(shinydashboard)
library(ggplot2)
library(VariantAnnotation)
library(Gviz)
library(GenomicRanges)

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
    # Subsetting VCF data based on the chromosome
    selectedChromosome <- input$chromosomeInput
    vcfSubset <- vcfFile()[which(seqnames(vcfFile()) == selectedChromosome)]
    
    # Counting SNVs and Indels for the selected chromosome
    snvCount <- sum(isSNV(vcfSubset))
    indelCount <- sum(isIndel(vcfSubset))
    
    #total
    totalSNVs <- sum(isSNV(vcfFile()))
    totalIndels <- sum(isIndel(vcfFile()))
    
    # Creating a data frame for plotting
    variantCounts <- data.frame(
      Variant = c("SNV", "Indel"),
      Count = c(snvCount, indelCount)
    )
    
    # Plotting 
    output$variantPlot <- renderPlot({
      ggplot(variantCounts, aes(x = Variant, y = Count, fill = Variant, label = Count)) +
        geom_bar(stat = "identity") +
        geom_text(position = position_stack(vjust = 0.5), color = "white", size = 5) +
        labs(title = paste("Variant Counts for Chromosome", selectedChromosome,
                           "\nTotal SNVs: ", totalSNVs, "\nTotal Indels: ", totalIndels),
             x = "Variant Type", y = "Count") +
        theme_minimal() +
        scale_fill_manual(values = c("#06A3DA", "#34AD54"))
    })
    
    
    # Updating the itrack based on the selected chromosome
    itrack <- IdeogramTrack(genome = "hg19", chromosome = selectedChromosome)
    vTrack <- AnnotationTrack(rowRanges(vcfSubset))
    
    output$plot3 <- renderPlot({
      plotTracks(list(itrack, vTrack))
    })
    
    # Displaying sample names 
    output$sampleNames <- renderText({
      paste(samples(header(vcfSubset)), collapse = ", ")
    })
  })
}




shinyApp(ui, server)