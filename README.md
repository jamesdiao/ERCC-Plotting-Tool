# Dimensionality Reduction Plotting Tool for the exRNA Atlas

#### James Diao, 21 March 2017

Here are some brief notes on the usage of this visualization tool. <br />

-------------------------------------------------------------

### ONLINE APPLICATION

The tool is primarily hosted at: http://exrna-atlas.org/exat/precomputedJobs. It is also available via shinyapps.io at: https://jamesdiao.shinyapps.io/ercc-plotting-tool/.

-------------------------------------------------------------

### CLONE REPO

Clone this source repository, containing our scripts and data: <br />
- `$ git clone https://github.com/jamesdiao/ERCC-Plotting-Tool.git`

-------------------------------------------------------------

### STARTING THE APP

#### If you have RStudio:  
 - Open `app.R` using RStudio. 
 - If you haven't already, install the R package "shiny" (must be connected to the Internet).  
`> install.packages("shiny")`  
 - Click "Run App" on the top right of the main window. The app should open and run immediately. <br />

#### If you do not have RStudio:  
 - Open an R console.  
 - If you haven't already, install the R package "shiny" (must be connected to the Internet).  
`> install.packages("shiny")`  
 - Run the following command (with the correct directory substituted for ~):  
`> shiny::runApp("~/ERCC-Plotting-Tool")`  
 - If you directly open `app.R` in R, you should already be at the correct directory. In that case, you only need to run:  
`> shiny::runApp()`


-------------------------------------------------------------

### APP DETAILS

#### Control Panel  
- **Plotting Style**: choose between the two R packages that can be used to generate the plot. ggplot2 gives a static image, while plotly features a suite of interactive functions (zooming, filtering, rotating 3D plots, etc.)
- **Dimension**: choose between 2 and 3 dimensions. 3D is only available when plotting with the Plotly package. If 2D is selected and 3 PCs are given, the last one is ignored. If 3D is selected and 2 PCs are given, the panel will display an empty plot. 
- **Embedding**: choose between the two embeddings: t-distributed stochastic neighbor embedding (tSNE) and principal components analysis (PCA). 
- **Principal Components**: choose which principal components to view, from 1-100. Only available for PCA.  
- **RNA Category**: choose from miRNA, piRNA, tRNA, and snRNA.  
- **Color By**: choose from Dataset, Condition, Anatomical\_Location, Biofluid\_Name, exRNA\_Source, Cell\_Culture\_Source, RNA\_Isolation\_Kit, ERCC\_QC\_Meets\_Standards. All labels are drawn from the exRNA atlas gridview at http://exrna-atlas.org/  

#### Filtering  
- **Datsets**: uncheck boxes to exclude all data in that dataset. 
- **Biofluids**: uncheck boxes to exclude all data of that biofluid.  

#### Generate Plots  
- **Make New Plot**: plots data points with the specified embedding, filtering, and labels.  All biofluids are labeled by shape.  
- **Download Plot**: downloads a pdf of the shown plot, dimensions 13x10. 

-------------------------------------------------------------

### DESCRIPTION OF FILES/FOLDERS
1. `app.R` contains all of the code for running the ERCC Plotting Tool. 
2. `Dependencies/` contains additional data required by `app.R`
 - `all_log_rna_reads.rds` contains a list of the log-transformed data frames (sample x RNA) for miRNA, piRNA, tRNA, and snRNA reads. 
 - `all_log_PCA_reduced.rds` contains a list of the top 100 PCA axes for miRNA, piRNA, tRNA, and snRNA. PCA was performed using the prcomp command from the stats R package. The data (in reads-per-million) was log-transformed after a small constant was added (to remove 0s). Before PCA, the data was centered but not rescaled. 
 - `all_log_tSNE_reduced.rds` contains the 2-dimensional tSNE embeddings for miRNA, piRNA, tRNA, and snRNA. tSNE was performed using the tsne R package with default settings. The embedding was computed on log-transformed RPM data after a small constant was added (to remove 0s). 
  - `all_3D_log_tSNE_reduced.rds` contains the 3-dimensional tSNE embeddings, computed analogously to the above. 
  - `Data_Summary_1567.tsv` contains the metadata collected from the exRNA atlas gridview. Rows are mapped to samples and used to assign labels. 
  - `map.rds` contains the mapping from samples to sample names in `Data_Summary_1567.tsv`. 
  - `sample_map.rds` contains dataset identity of each sample.
3. `Make_Dependency_Files/` contains a `Make_Dependencies.Rmd`, which can be knitted to automatically generate all .rds dependencies. This folder also contains the original datasets from the ExRNA Atlas in `miRNA/`, `piRNA/`, `tRNA/`, `smallRNAQuants/`, and `QC_Results/`. 

-----------------------------------------------------------------

### CONTACT  

Please contact James Diao (james.diao@yale.edu) with any questions.

<br />
<br />
<br />


