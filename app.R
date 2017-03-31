#rsconnect::deployApp('/Users/jamesdiao/Documents/Gerstein/ERCC-Plotting-Tool/')
#setwd("/Users/jamesdiao/Documents/Gerstein/ERCC-Plotting-Tool")

require(plotly)
require(shinysky)
require(ggplot2)
require(dplyr)
require(shiny)
require(tsne)
options(warn = -1) 

path <- "Dependencies"

### LOAD FILES
log_rna_reads <- readRDS(sprintf("%s/all_log_rna_reads.rds", path))
all_reads_pca <- readRDS(sprintf("%s/all_log_PCA_reduced.rds", path))
all_reads_tsne_2d <- readRDS(sprintf("%s/all_log_tSNE_reduced.rds", path))
all_reads_tsne_3d <- readRDS(sprintf("%s/all_3D_log_tSNE_reduced.rds", path))
sample_map <- readRDS(sprintf("%s/sample_map.rds", path))
map <- readRDS(sprintf("%s/map.rds", path))

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
                                                                sub("Culture Media, Conditioned","Cultured Media", Biofluid_Name))))
  biofluid_names <- table(data_summary$Biofluid_Name) %>% sort(decreasing = T) %>% names
  data_summary <- mutate(data_summary, Biofluid_Name = factor(Biofluid_Name, levels = biofluid_names))
  locs <- apply(data_summary, 1, function(row) any(is.na(row)))
  data_summary[locs,] <- replace(data_summary[locs,], is.na(data_summary[locs,]), 0) 
  for (i in 1:ncol(data_summary)) {
    if (is.character(data_summary[,i])) {
      data_summary[,i] <- substr(data_summary[,i], 1, 29)
    }
  }
  return(data_summary)
}
data_summary <- import_tsv(sprintf("%s/Metadata_Summary.tsv", path)) %>% arrange(Biosample_Metadata_Accession)
#850, 1369, 1567

biofluid <- data_summary[map,]$Biofluid_Name
condition <- data_summary[map,]$Condition
anatomical <- data_summary[map,]$Anatomical_Location
exRNA_src <- data_summary[map,]$exRNA_Source
cell_src <- data_summary[map,]$Cell_Culture_Source
profiling <- data_summary[map,]$Profiling_Assay
rna_kit <- data_summary[map,]$RNA_Isolation_Kit
qc_std <- data_summary[map,]$ERCC_QC_Meets_Standards
gen_reads <- data_summary[map,]$ERCC_QC_Reference_Genome_Reads
tran_reads <- data_summary[map,]$ERCC_QC_Transcriptome_Reads
gt_ratio <- data_summary[map,]$ERCC_QC_Transcriptome_Genome_Ratio

hover_text <- sprintf('Biofluid: %s </br>Dataset: %s </br>Condition: %s </br>
                      Anatomy: %s </br>exRNA Source: %s </br>
                      Cell Source: %s </br>Profiling Assay: %s </br>
                      RNA Isolation Kit: %s </br> QC Standard: %s </br>Genome Reads: %s </br>
                      Trnscpt Reads: %s </br>Trnscpt/Genome Ratio: %s',
                      biofluid, 
                      abbreviate(sample_map,minlength = 20, method = 'both.sides'), 
                      condition, anatomical, exRNA_src, cell_src, profiling, 
                      rna_kit, qc_std, gen_reads, tran_reads, gt_ratio)
#hover_text <- sprintf('%s (%s)', 
#                      abbreviate(sample_map,minlength = 20, method = 'both.sides'), 
#                      biofluid)

plottable <- gsub("_"," ","Dataset" %>% 
                    c(colnames(data_summary[map,])[apply(data_summary[map,], 2, function(col) length(unique(col))) %>% between(2,20)]))

### PRINCIPAL COMPONENTS ANALYSIS
pca_plot <- function(axis_x, axis_y, biofluid, color_elements, keep, pca_object, smRNA, colorby, color_set) {
  data.frame(PCA_1 = pca_object[keep,axis_x], PCA_2 = pca_object[keep,axis_y], 
             Shape = biofluid[keep], Color = color_elements[keep]) %>% 
    ggplot(aes(x = PCA_1, y = PCA_2, shape = Shape, color = Color)) + 
    geom_point(size = 2.5) + 
    ggtitle(sprintf("PCA Plot of %s Colored by %s (%s Samples)", smRNA, gsub("_"," ",colorby), sum(keep))) + 
    xlab(paste0("PC ",axis_x)) + ylab(paste0("PC ",axis_y)) + 
    theme(plot.title = element_text(size=22,face="bold"), 
          axis.title=element_text(size=16), 
          axis.text=element_text(size=11),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=11)) + 
    scale_shape_manual(values = c(0, 16:17, 8:10, 2, 16, 1)) +
    scale_color_manual(values = color_set)
}

pca_plotly <- function(axis_x, axis_y, color_elements, keep, pca_object, smRNA, colorby) {
  data.frame(PCA_1 = pca_object[keep,axis_x], PCA_2 = pca_object[keep,axis_y], Color = color_elements[keep]) %>% 
  plot_ly(x = ~PCA_1, y = ~PCA_2, color = ~Color, 
               hoverinfo = 'text', text = hover_text[keep],
               marker = list(size = 6, symbol = 'square')) %>%
    add_markers() %>%
    layout(title = sprintf("PCA Plot of %s Colored by %s (%s Samples)", 
                           smRNA, gsub("_"," ",colorby), sum(keep)),
           xaxis = list(title = sprintf('PC %s', axis_x)),
           yaxis = list(title = sprintf('PC %s', axis_y)))
}

pca_plotly_3d <- function(axis_x, axis_y, axis_z, color_elements, keep, pca_object, smRNA, colorby) {
  data.frame(PCA_1 = pca_object[keep,axis_x], PCA_2 = pca_object[keep,axis_y], PCA_3 = pca_object[keep, axis_z], 
             Color = color_elements[keep]) %>%
  plot_ly(x = ~PCA_1, y = ~PCA_2, z = ~PCA_3, color = ~Color, 
          hoverinfo = 'text', text = hover_text[keep],
          marker = list(size = 5, symbol = 'square')) %>%
    add_markers() %>%
    layout(title = sprintf("PCA Plot of %s Colored by %s (%s Samples)", 
                           smRNA, gsub("_"," ",colorby), sum(keep)),
           scene = list(
             xaxis = list(title = sprintf('PC %s', axis_x)),
             yaxis = list(title = sprintf('PC %s', axis_y)),
             zaxis = list(title = sprintf('PC %s', axis_z))
           ))
}

### t-DISTRIBUTED STOCHASTIC NEIGHBOR EMBEDDING
tsne_plot <- function(biofluid, color_elements, keep, tsne_object, smRNA, colorby, color_set) {
  data.frame(tsne_object[keep,], Shape = as.factor(biofluid[keep]), Color = color_elements[keep]) %>%
    ggplot(aes(x = tSNE_1, y = tSNE_2, shape = Shape, color = Color)) + 
    geom_point(size = 2.5) + 
    ggtitle(sprintf("tSNE Plot of %s Colored by %s (%s Samples)", smRNA, gsub("_"," ",colorby), sum(keep))) + 
    xlab("tSNE 1") + ylab("tSNE 2") + 
    theme(plot.title = element_text(size=22,face="bold"), 
          axis.title=element_text(size=16), 
          axis.text=element_text(size=11),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=11)) + 
    scale_shape_manual(values = c(0, 16:17, 8:10, 2, 16, 1)) +
    scale_color_manual(values = color_set)
}

tsne_plotly <- function(color_elements, keep, tsne_object, smRNA, colorby) {
  data.frame(tSNE_1 = tsne_object[keep,1], tSNE_2 = tsne_object[keep,2], 
             Color = color_elements[keep]) %>% 
    plot_ly(x = ~tSNE_1, y = ~tSNE_2, color = ~Color,  
            hoverinfo = 'text', text = hover_text[keep],
            marker = list(size = 6, symbol = 'square')) %>%
    add_markers() %>%
    layout(title = sprintf("tSNE Plot of %s Colored by %s (%s Samples)", 
                           smRNA, gsub("_"," ",colorby), sum(keep)),
           xaxis = list(title = 'tSNE 1'),
           yaxis = list(title = 'tSNE 2'))
}

tsne_plotly_3d <- function(color_elements, keep, tsne_object, smRNA, colorby) {
  data.frame(tSNE_1 = tsne_object[keep,1], tSNE_2 = tsne_object[keep,2], tSNE_3 = tsne_object[keep,3], 
             Color = color_elements[keep]) %>% 
    plot_ly(x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, color = ~Color,  
            hoverinfo = 'text', text = hover_text[keep],
            marker = list(size = 5, symbol = 'square')) %>%
    add_markers() %>%
    layout(title = sprintf("tSNE Plot of %s Colored by %s (%s Samples)", 
                           smRNA, gsub("_"," ",colorby), sum(keep)),
           scene = list(
             xaxis = list(title = 'tSNE 1'),
             yaxis = list(title = 'tSNE 2'),
             zaxis = list(title = 'tSNE 3')
           ))
}

data_opts <- unique(sample_map) %>% setNames(unique(sample_map)) %>% as.list
  new_data_opts <- sprintf("%s (%s)", unlist(data_opts), table(sample_map)) %>% as.list
biofluid_opts <- levels(data_summary$Biofluid_Name)
  biofluid_opts <- biofluid_opts %>% setNames(biofluid_opts) %>% as.list
  new_biofluid_opts <- sprintf("%s (%s)", unlist(biofluid_opts), table(biofluid)) %>% as.list

ui <- shinyUI(fluidPage(
  
  titlePanel("Dimensionality Reduction Plotting Tool for the exRNA Atlas"),
  h4("James Diao, Version 1.0.3"),
  h5(a("https://github.com/jamesdiao/ERCC-Plotting-Tool", href="https://github.com/jamesdiao/ERCC-Plotting-Tool", target="_blank")),
  fluidRow(
    column(4,
           h3("Control Panel"),
           wellPanel(
             radioButtons(inputId = "plotstyle", label = "Plotting Style", 
                          inline = T, choices = c("ggplot2", "plotly"), selected = "ggplot2"),
             conditionalPanel(
               condition = "input.plotstyle == 'plotly'",
               radioButtons(inputId = "dim", label = "Dimension", 
                            inline = T, choices = c("2D", "3D"), selected = "2D")
             ),
             radioButtons(inputId = "embedding", label = "Embedding", 
                          inline = T, choices = c("PCA", "tSNE"), selected = "PCA"),
             conditionalPanel(
               condition = "input.embedding == 'PCA' & input.dim == '2D'",
               checkboxGroupInput(inputId = "pcs_2d", label = "Principal Components", 
                            inline = T, choices = sprintf("PC%s",1:5), selected = sprintf("PC%s",1:2))
               #selectizeInput(inputId = "pcs_2d", label = "Principal Components", 
              #                choices = sprintf("PC%s",1:10), selected = sprintf("PC%s",1:2), multiple = TRUE,
              #                options = list(maxItems = 2))
             ), 
             conditionalPanel(
               condition = "input.embedding == 'PCA' & input.dim == '3D'",
               checkboxGroupInput(inputId = "pcs_3d", label = "Principal Components", 
                           inline = T, choices = sprintf("PC%s",1:5), selected = sprintf("PC%s",1:3))
               #selectizeInput(inputId = "pcs_3d", label = "Principal Components", 
              #                choices = sprintf("PC%s",1:10), selected = sprintf("PC%s",1:3), multiple = TRUE,
              #                options = list(maxItems = 3))
             )
           ),
           wellPanel(
             radioButtons(inputId = "smRNA", label = "RNA Category", 
                          choices = c("miRNA", "piRNA", "tRNA", "snRNA")),
             radioButtons(inputId = "colorby", label = "Color By", 
                          choices = plottable)
           ),
           h3("Filtering"),
           #actionButton(inputId = 'recompute', label = 'Recompute Values'),
           #busyIndicator("In Progress: Please Wait", wait = 500),
           h4(),
           wellPanel(
             checkboxGroupInput("checkdata", label = "Datasets", 
                                choices = new_data_opts,
                                selected = new_data_opts)
           ),
           wellPanel(
             checkboxGroupInput("checkfluid", label = "Biofluids", 
                                choices = new_biofluid_opts,
                                selected = new_biofluid_opts)
           ),
           tags$head(tags$style("#plot_out{height:80vh !important;}")),
           tags$head(tags$style("#plotly_out{height:80vh !important;}"))
    ),
    column(8, h3("Generate Plots"),
           wellPanel(
             conditionalPanel(
               condition = "input.plotstyle == 'plotly'",
               actionButton(inputId = "run_plotly", label = "Make New Plot ", 
                            icon = icon("bar-chart"), styleclass = "success")
             ),
             conditionalPanel(
               condition = "input.plotstyle == 'ggplot2'",
               actionButton(inputId = "run_ggplot2", label = "Make New Plot ", 
                            icon = icon("bar-chart"), styleclass = "success"),
               downloadButton(outputId = "plot_down", label = "Download Plot")
             )
           ),
           shinyalert(id = 'alert', click.hide = TRUE, auto.close.after = 5),
           conditionalPanel(
             condition = "input.plotstyle == 'ggplot2'",
             plotOutput("plot_out")
           ),
           conditionalPanel(
             condition = "input.plotstyle == 'plotly'",
             plotlyOutput("plotly_out")
           )
    )
  )
))


server <- shinyServer(function(input, output, session) {
  
  coord <- reactiveValues()
  coord$pca <- all_reads_pca
  coord$tsne_2d <- all_reads_tsne_2d
  coord$tsne_3d <- all_reads_tsne_3d
  coord$keep <- rep(TRUE,length(sample_map))
  
  rv2 <- reactiveValues(lstval2=1:2,curval2=1:2)
  rv3 <- reactiveValues(lstval3=1:3,curval3=1:3)
  
  observeEvent(input$pcs_2d, {
    rv2$lstval <- rv2$curval
    rv2$curval <- input$pcs_2d %>% substr(3,3) %>% as.integer 
    if(length(rv2$curval) > 2) {
      updateCheckboxGroupInput(session, "pcs_2d", 
        selected = sprintf("PC%s", rv2$curval[!(rv2$curval %in% rv2$lstval)] %>% 
          c(sort(rv2$lstval)) %>% head(2) %>% sort))
    }
  })
  
  observeEvent(input$pcs_3d, {
    rv3$lstval <- rv3$curval
    rv3$curval <- input$pcs_3d %>% substr(3,3) %>% as.integer 
    if(length(rv3$curval) > 3) {
      updateCheckboxGroupInput(session, "pcs_3d", 
        selected = sprintf("PC%s", rv3$curval[!(rv3$curval %in% rv3$lstval)] %>% 
          c(sort(rv3$lstval)) %>% head(3) %>% sort))
    }
  })
  
  observeEvent(input$plotstyle, {
    if (input$plotstyle == 'ggplot2') {
      updateRadioButtons(session, "dim", selected = '2D')
    }
  })
  
  observeEvent(input$run_ggplot2 | input$run_plotly, {
    hideshinyalert(session, id = "alert")
  if (input$embedding == "PCA" & (
    (length(input$pcs_2d) < 2 & input$dim == '2D') | (length(input$pcs_3d) < 3 & input$dim == '3D')
    )) {
    showshinyalert(session, id = "alert", HTMLtext = 'Please select more principal components', 
                   styleclass = 'danger')
  } else {
    find_data <- function(original, new) {
      sapply(original, function(data) any(grepl(data,new)))
    }
    
    keep_data <- sample_map %>% find_data(input$checkdata)
    keep_biofluid <- biofluid %>% find_data(input$checkfluid)
    coord$keep <- keep_data & keep_biofluid
    
    new_data_opts <- unlist(data_opts) %>% setNames(NULL)
    new_data_opts <- sprintf("%s (%s)", new_data_opts, 
                             sapply(new_data_opts, function(data) sum(sample_map[coord$keep] == data)))
    data_remain <- unique(sample_map) %>% find_data(input$checkdata)
    updateCheckboxGroupInput(session, 'checkdata', choices = as.list(new_data_opts), 
                             selected = as.list(new_data_opts[data_remain]))
    
    new_biofluid_opts <- unlist(biofluid_opts) %>% setNames(NULL)
    new_biofluid_opts <- sprintf("%s (%s)", new_biofluid_opts, 
                             sapply(new_biofluid_opts, function(data) sum(biofluid[coord$keep] == data)))
    biofluid_remain <- levels(biofluid) %>% find_data(input$checkfluid)
    updateCheckboxGroupInput(session, 'checkfluid', choices = as.list(new_biofluid_opts), 
                             selected = as.list(new_biofluid_opts[biofluid_remain]))
    
    colorby <- gsub(" ","_",input$colorby)
    
    if (colorby == "Dataset"){
      color_elements <- sample_map
    } else {
      color_elements <- data_summary[,colorby][map]
    }
    
    set <- unique(color_elements)
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    color_set <- gg_color_hue(n = length(set)) %>% setNames(set)
    color_set <- color_set[set %in% unique(color_elements[coord$keep])]
    color_names <- names(color_set)
    if ( ("YES" %in% toupper(color_names)) & ("NO" %in% toupper(color_names)) )
      names(color_set) <- rev(color_names)
    
    if (input$embedding == "PCA") {
      reads_pca <- coord$pca[[input$smRNA]]
      pcs <- gsub("[^0-9]", "", unlist(input$pcs_2d) ) %>% as.numeric
      if (input$plotstyle == "plotly") {
        if (input$dim == "2D") {
          plotly_out <- pca_plotly(pcs[1], pcs[2], color_elements, coord$keep, reads_pca, input$smRNA, colorby)
        }
        if (input$dim == "3D") {
          pcs <- gsub("[^0-9]", "", unlist(input$pcs_3d) ) %>% as.numeric
          plotly_out <- pca_plotly_3d(pcs[1], pcs[2], pcs[3], color_elements, coord$keep, reads_pca, input$smRNA, colorby)
        }
        output$plotly_out <- renderPlotly({ plotly_out })
      }
      if (input$plotstyle == "ggplot2") {
        plot_out <- pca_plot(pcs[1], pcs[2] + (pcs[1] == pcs[2]), biofluid, color_elements, 
                             coord$keep, reads_pca, input$smRNA, colorby, color_set)
        output$plot_out <- renderPlot({ plot_out })
        
      }
    } 
    
    if (input$embedding == "tSNE") {
      reads_tsne_2d <- coord$tsne_2d[[input$smRNA]]
      reads_tsne_3d <- coord$tsne_3d[[input$smRNA]]
      if (input$plotstyle == "plotly") {
        if (input$dim == "2D") {
          plotly_out <- tsne_plotly(color_elements, coord$keep, reads_tsne_2d, input$smRNA, colorby)
        }
        if (input$dim == "3D") {
          plotly_out <- tsne_plotly_3d(color_elements, coord$keep, reads_tsne_3d, input$smRNA, colorby)
        }
        output$plotly_out <- renderPlotly({ plotly_out })
      }
      if (input$plotstyle == "ggplot2") {
        plot_out <- tsne_plot(biofluid, color_elements, coord$keep, 
                                  reads_tsne_2d, input$smRNA, colorby, color_set)
        output$plot_out <- renderPlot({ plot_out })
      }
    }
    
    output$plot_down <- downloadHandler(
      filename = function() {
        sprintf("%s_Plot_%s_%s.pdf", input$embedding, input$smRNA, colorby)
      },
      content = function(file) {
        pdf(file, onefile = TRUE, width = 13, height = 10)
        if (exists("plot_out"))
          print(plot_out)
        dev.off()
      }
    )
  }
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)






#observeEvent(input$recompute, {
#  
#  keep_data <- sample_map %in% input$checkdata
#  keep_biofluid <- biofluid %in% input$checkfluid
#  coord$keep <- keep_data & keep_biofluid

#  coord$pca <- lapply(log_rna_reads, function(reads) {
#    rna_pca <- prcomp(reads[coord$keep,], center = T, scale. = F)
#    reduced_cols <- seq(1,min(100,ncol(rna_pca$x)))
#    return(rna_pca$x[,reduced_cols])
#  })

#  coord$tsne <- lapply(log_rna_reads, function(reads) {
#    log_tsne_out <- tsne(reads[coord$keep,], k=2, initial_dims = 30, perplexity = 30, 
#                         max_iter = 500, epoch = 50)
#    return(data.frame("tSNE_1" = log_tsne_out[,1], "tSNE_2" = log_tsne_out[,2]))
#  })

#})







