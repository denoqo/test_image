#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("R/global.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Gene Signatures Expressed Amplification"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("tcga_code",label = "Cancer Type", choices = tcga_codes, selected="chol_tcga"),
      selectInput("gene_scope", label = "Gene Scope", choices = c("Oncogene"="onco", "genome-wide" =  "gw"), selected="gw"),
      selectInput("mols_id", label = "Data type", multiple = T,
                  selected = c("gistic", "rppa_Zscores"), 
                  choices =c("gistic", "rna_seq_v2_mrna_median_Zscores", "rppa_Zscores")),
      selectInput("clin_selected", label = "Clinical Features",
                  selected = c(
                    "Overall Survival Status",
                    "Overall Survival (Months)",
                    "TMB (nonsynonymous)", 
                    "Fraction Genome Altered", 
                    "Neoplasm Histologic Grade",
                    "Sex"), 
                  multiple = TRUE, selectize = T,
                  choice = clin_colnames),
      selectInput("cna_cut", label = "CNA(levels) >=", choices = c(0,1,2), selected=1),
      selectInput("mrna_cut", label = "mRNA(z-score) >=", choices = c(0,1,2), selected=2),
      selectInput("rppa_cut", label = "Protein(z-score) >=", choices = c(0,1,2), selected=2),
    ),
    
    mainPanel(
      
      tabsetPanel(
        tabPanel(
          title = "Samples",
          uiOutput("gene_summary"),
          hr(),h1("Survival Curves"),
          plotOutput("survialPlot"),
          hr(),h1("Sample Groups"),
          fluidRow(
            DTOutput("tableSampleGroup") %>% column(width = 6),
            plotOutput("plotSampleGroup") %>% column(width = 6)
          )
          
          
        ),
        
        
        tabPanel(
          title = "Clinical",
          shinysky::busyIndicator(text = "Loading...", wait = 0),
          hr(),h1("Clinical Summary"),
          gt::gt_output("gt_asso"),
          hr(),h1("Clinical Visualization"),
          uiOutput("ui_boxbarPlot")
        ),
        
        
        tabPanel(
          title = "Signature",
          hr(),h1("Signature Visualization"),
          radioButtons("sortBand", "Sort by Chr?", 
                       selected = F, choices= c(F,T)),
          uiOutput("uiOncoSig"),
          hr(),h1("Signature list"),
          DTOutput("tableSig"),
          hr(),h1("Chromosome bands"),
          tableOutput("bandSig")
        )
      )
      
      
    )
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$ui_boxbarPlot <- renderUI({
    plotOutput(
      "boxbarPlot",
      height = paste0(ceiling(length(input$clin_selected) / 3) * 300, "px") )
  })
  
  tcga_code <- reactive({input$tcga_code})
  
  samples <- reactive({
    GetSample(input$tcga_code)
  })
  genes <- reactive({
    GetGene(input$gene_scope)
  })
  sig <- reactive({
    if(input$gene_scope == "onco") {
      sig <- GetGeneSignature(
        tcga_code = input$tcga_code, 
        gene_list = genes(),
        mols_id = input$mols_id,
        cna_cut = input$cna_cut,
        mrna_cut = input$mrna_cut,
        rppa_cut = input$rppa_cut
      )
    } else if (input$gene_scope == "gw") {
      sig <- GetGeneSignatureGW(
        tcga_code = input$tcga_code, 
        mae_list_gw = mae_list_gw, 
        mols_id = input$mols_id,
        cna_cut = input$cna_cut,
        mrna_cut = input$mrna_cut,
        rppa_cut = input$rppa_cut
      )
    }
    return(sig)
  })
  sample_group <- reactive({
    GetSampleGroup(sample_all = samples(), sig = sig())
  })
  clin_specific <- reactive({
    print("Filtering samples")
    clin_subset <- clin %>% filter(sample_id %in% samples()) 
    
    print("Binding group information")
    clin_group <- clin_subset %>%  right_join(sample_group())
    print("Finish Binding group information")
    
    return(clin_group)
  })
  group_clin <- reactive({
    # select clincial feature
    select(clin_specific(), group, sample_id, all_of(input$clin_selected))
  })
  
  
  output$gene_summary <- renderUI({
    ngenes <- length(unique(sig()[["gene"]]))
    nsampl <- length(unique(sig()[["sample_id"]]))
    totalSamples <- nrow(sample_group())
    
    return(
      fluidRow(align = "center",
               hr(),h1("Brief Suammry"),
               h4(glue::glue("There are {ngenes} genes with alterations in {nsampl} samples. ")),
               h4(glue::glue("The total sample number is {totalSamples} samples.")),
               
      )
    )
  })
  
  output$boxbarPlot <- renderPlot(res = 96, {
    PlotAssociation(group_clin())
  })
  output$survialPlot <- renderPlot(res = 96,{
    plot <- PerformSurvival(sample_group(), clinical = clin_cleanname,os_only = F)
    return(plot)
    print("Finish survialPlot")
  })
  output$gt_asso <- gt::render_gt({
    GetGTsummary(group_clin())
  })
  
  output$tableSig <- renderDT({
    sig()
  })
  output$oncoSig <- renderPlot(res=96,{
    ht = DrawOnco(sig(),samples(),annorow = T, sortBand = input$sortBand)
    return(ht)
  })
  output$uiOncoSig <- renderUI({
    lengene <- length(unique(sig()$gene))
    hei = ifelse(lengene < 100, lengene * 25, 2000)
    wid = sprintf("%spx", length(samples()) * 16+100)
    
    
    plotOutput("oncoSig", width = "100%",
               height = sprintf("%spx", hei))
  })
  output$bandSig <- renderTable({
    # sig
    siganno <- sig() %>% left_join(annodata)
    mat_band <- dplyr::count(siganno, sample_id, band) %>% 
      spread(band, n)
    return(mat_band)
  })
  
  
  output$tableSampleGroup <- renderDT({
    sample_group()
  })
  output$plotSampleGroup <- renderPlot(res = 96, {
    BarplotCount(dplyr::count(sample_group(), group))
  })
  
  
  
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
