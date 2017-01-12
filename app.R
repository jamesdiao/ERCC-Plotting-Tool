# Set working directory to file folder
outdir <- getSrcDirectory(function(dummy) {dummy})
setwd(outdir)

#rsconnect::deployApp('/Users/jamesdiao/Box Sync/Gerstein/ERCC-Plotting-Tool/')

# Install and load all required packages
pkg_list <- c("ggplot2","dplyr","shiny")
installed <- pkg_list %in% installed.packages()[,"Package"]
if (!all(installed))
  install.packages(pkg_list[!installed])
sapply(pkg_list, require, character.only = T)

### LOAD FILES
all_reads_pca <- readRDS("Dependencies/all_log_PCA_reduced.rds")
all_reads_tsne <- readRDS("Dependencies/all_log_tSNE_reduced.rds")
sample_map <- readRDS("Dependencies/sample_map.rds")
map <- readRDS("Dependencies/map.rds")

### IMPORT TSV
import_tsv <- function(data_file) {
  file.by.line <- suppressWarnings(readLines(data_file) %>% strsplit("\t"))
  data_summary <- do.call(rbind, file.by.line) %>% as.data.frame(stringsAsFactors = F)
  if (ncol(data_summary) == 19) {
    colnames(data_summary) <- c("Blank","Biosample_Name","Condition",
      "Anatomical_Location","Biofluid_Name","exRNA_Source","Cell_Culture_Source",
      "Profiling_Assay","RNA_Isolation_Kit","ERCC_QC_Meets_Standards",
      "ERCC_QC_Reference_Genome_Reads","ERCC_QC_Transcriptome_Genome_Ratio",
      "ERCC_QC_Transcriptome_Reads","Download_Data","Download_Advanced_Results",
      "Download_Metadata","RNA_Profile","External_References","Biosample_Metadata_Accession")
  } else {
    colnames(data_summary) <- c("Blank","Biosample_Name","Condition","Anatomical_Location",
      "Biofluid_Name","exRNA_Source","Cell_Culture_Source",
      "Profiling_Assay","RNA_Isolation_Kit","ERCC_QC_Meets_Standards",
      "ERCC_QC_Reference_Genome_Reads","ERCC_QC_Transcriptome_Genome_Ratio",
      "ERCC_QC_Transcriptome_Reads","Download_Data","Download_Metadata",
      "External_References","Biosample_Metadata_Accession")
  }
  data_summary <- data_summary[,apply(data_summary, 2, function(col) any(nzchar(col)))]
  data_summary <- suppressWarnings(data_summary %>% 
      mutate(ERCC_QC_Reference_Genome_Reads = gsub(",","",ERCC_QC_Reference_Genome_Reads) %>% as.integer) %>%
      mutate(ERCC_QC_Transcriptome_Reads = gsub(",","",ERCC_QC_Transcriptome_Reads) %>% as.integer) %>%
      mutate(ERCC_QC_Transcriptome_Genome_Ratio = ERCC_QC_Transcriptome_Genome_Ratio %>% as.numeric) %>% 
      mutate(Biofluid_Name = sub("Cerebrospinal fluid","CSF", 
             sub("Culture Media, Conditioned","Cultured_Media", Biofluid_Name))))
  biofluid_names <- table(data_summary$Biofluid_Name) %>% sort(decreasing = T) %>% names
  data_summary <- mutate(data_summary, Biofluid_Name = factor(Biofluid_Name, levels = biofluid_names))
}
data_summary <- import_tsv("Dependencies/Data_Summary_1567.tsv") %>% arrange(Biosample_Metadata_Accession)
#850, 1369, 1567

plottable <- "Dataset" %>% 
  c(colnames(data_summary[map,])[apply(data_summary[map,], 2, function(col) length(unique(col))) %>% between(2,20)])


### PRINCIPAL COMPONENTS ANALYSIS
pca_plot <- function(axis_x, axis_y, biofluid, color_by, keep, pca_object, smRNA, colorby) {
  data.frame(PCA_1 = pca_object[,axis_x], PCA_2 = pca_object[,axis_y], 
             Shape = biofluid, Color = color_by)[keep,] %>% 
    ggplot(aes(x = PCA_1, y = PCA_2, shape = Shape, color = Color)) + 
    geom_point(size = 2.5) + 
    ggtitle(sprintf("PCA Plot of %s Colored By %s", smRNA, colorby)) + 
    xlab(paste0("PCA_",axis_x)) + ylab(paste0("PCA_",axis_y)) + 
    theme(plot.title = element_text(size=22,face="bold"), 
          axis.title=element_text(size=16), 
          axis.text=element_text(size=11),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=11)
    )
}
### t-DISTRIBUTED STOCHASTIC NEIGHBOR EMBEDDING
tsne_plot <- function(biofluid, color_by, keep, tsne_object, smRNA, colorby) {
  temp <- data.frame(tsne_object, Shape = biofluid, Color = color_by)[keep,] %>%
    ggplot(aes(x = tSNE_1, y = tSNE_2, shape = Shape, color = Color)) + 
    geom_point(size = 2.5) + 
    ggtitle(sprintf("tSNE Plot of %s Colored By %s", smRNA, colorby)) + 
    theme(plot.title = element_text(size=22,face="bold"), 
          axis.title=element_text(size=16), 
          axis.text=element_text(size=11),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=11)
    )
}

data_opts <- unique(sample_map) %>% setNames(unique(sample_map)) %>% as.list
biofluid_opts <- levels(data_summary$Biofluid_Name)
biofluid_opts <- biofluid_opts %>% setNames(biofluid_opts) %>% as.list

ui <- shinyUI(fluidPage(
   
   titlePanel("Plotting Tool for 1075 Samples from the exRNA Atlas"),
   h4("James Diao, 12 January 2017"),
   h5("https://jamesdiao.shinyapps.io/ercc-plotting-tool"),
   fluidRow(
      column(4,
        h3("Control Panel"),
          wellPanel(
            radioButtons(inputId = "plotstyle", label = "Plotting Style", 
                         choices = c("tSNE", "PCA"), selected = "tSNE"),
          
            sliderInput(inputId = "pcs", label = "Principal Components (PCA Only)",
                        min = 1, max = 5, value = c(1,2), dragRange = T, ticks = F)
          ),
          wellPanel(
            radioButtons(inputId = "smRNA", label = "RNA Category", 
                         choices = c("miRNA", "piRNA", "tRNA", "snRNA")),
            radioButtons(inputId = "colorby", label = "Color By", 
                         choices = plottable)
          ),
          h3("Filtering"),
          wellPanel(
            checkboxGroupInput("checkdata", label = "Datasets", 
                               choices = data_opts,
                               selected = data_opts)
          ),
          wellPanel(
            checkboxGroupInput("checkfluid", label = "Biofluids", 
                               choices = biofluid_opts,
                               selected = biofluid_opts)
          ),
          tags$head(tags$style("#plot_out{height:80vh !important;}"))
      ),
      column(8, h3("Generate Plots"),
        wellPanel(
          actionButton(inputId = "run", label = "Make New Plot", 
                       icon = icon("bar-chart"), styleclass = "success"),
          
          downloadButton(outputId = "plot_down", label = "Download Plot")
        ),
         plotOutput("plot_out")
      )
   )
))


server <- shinyServer(function(input, output, session) {
  observeEvent(input$run, {
    keep_data <- sample_map %in% input$checkdata
    biofluid <- data_summary[map,]$Biofluid_Name
    keep_biofluid <- biofluid %in% input$checkfluid
    keep <- keep_data & keep_biofluid
    if (input$plotstyle == "PCA") {
      reads_pca <- all_reads_pca[[input$smRNA]]
      pcs <- input$pcs
      if (pcs[1] == pcs[2])
        pcs[2] <- pcs[1] + 1
      if (input$colorby == "Dataset"){
        plot_out <- pca_plot(pcs[1], pcs[2], biofluid, sample_map,
                             keep, reads_pca, input$smRNA, input$colorby)
      } else {
        plot_out <- pca_plot(pcs[1], pcs[2], biofluid, data_summary[,input$colorby][map], 
                             keep, reads_pca, input$smRNA, input$colorby)
      }

    } else {
      reads_tsne <- all_reads_tsne[[input$smRNA]]
      if (input$colorby == "Dataset"){
        plot_out <- tsne_plot(biofluid, sample_map, 
                              keep, reads_tsne, input$smRNA, input$colorby) 
      } else {
        plot_out <- tsne_plot(biofluid, data_summary[,input$colorby][map], 
                              keep, reads_tsne, input$smRNA, input$colorby)
      }
    }
    
    output$plot_out <- renderPlot({ plot_out })
    output$plot_down <- downloadHandler(
      filename = function() {
        sprintf("%s_Plot_%s_%s.pdf", input$plotstyle, input$smRNA, input$colorby)
      },
      content = function(file) {
        pdf(file, onefile = TRUE, width = 13, height = 10)
        if (exists("plot_out"))
          print(plot_out)
        dev.off()
      }
    )
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)




