################## WORKFLOW ########################
# 1. Create Table one, stratified by group (SC only, Bulk + SC, or Bulk only)
# 2. Read in sc data, build pseudo-bulk data, GT proportions
# 3. Build signature matrices, optimize marker selection
# 4. Deconvolution of SC+Bulk data, bias analysis
# 5. Deconvolution of real bulk data, visual analysis
# 7. Predictions, random forest optimization

# Later: SURVIVAL ANLAYSIS? TIME TO DEATH?

############################ BEGIN WORKFLOW ####################################

##### Install Libraries #####
# Create a vector with the package names
packages <- c("devtools", "Seurat", "tidyverse", "reshape2", "readxl", 
              "biotools", "knitr", "gtsummary", "curl", "Matrix", "flextable", 
              "purrr", "jsonlite", "officer", "xCell", "MCPcounter",
              "EPIC", "FARDEEP", "pROC", "glmnet", "randomForest",
              "conflicted")

# Loop through the packages and check if they are installed
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    if(package == "EPIC"){
      devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
    } else if("xCell"){
      devtools::install_github('dviraran/xCell')
    } else if("MCPcounter"){
      install_github("ebecht/MCPcounter",ref="master", subdir="Source")
    } else {
      # If the package is not installed, install it
      install.packages(package, dependencies = TRUE)
    }
  }
  # Load the package
  library(package, character.only = TRUE)
}

# Load Cibersort
source('CIBERSORT.R')

# Conflicting function preferences
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "base")
conflicts_prefer(base::setdiff)
conflicts_prefer(reshape2::melt)
conflicts_prefer(base::as.data.frame)
conflicts_prefer(readxl::read_xlsx)
conflicts_prefer(gtsummary::as_flextable)

##### File locations #######
# Working directory
setwd("/projects/b1036/davidhparkinson")

### On Quest ###
myfolderpath <- "/projects/b1036/davidhparkinson"
sc_datapath <- "/projects/b1036/zzreinstein/mcc_sc_tissue_oct22/all_mcc.rds"
bulk_datapath <- "/projects/b1036/davidhparkinson/Bulk+SC Info/bulk_data.rds"
sc_metadatapath <- "/projects/b1036/davidhparkinson/Bulk+SC Info/sc_metadata.xlsx"
bulk_metadatapath <- "/projects/b1036/davidhparkinson/Bulk+SC Info/bulk_metadata.xlsx"

### On Personal Computer ###
#myfolderpath <- "~/THESIS"
#scdatapath <- "~/Quest Files/all_mcc.rds"
#bulkdatapath <- "~/Quest Files/GRCh38+MCPyV_2pass_20230113.rds"
#metadatapath <- "~/Quest Files/sc_metadata.xlsx"
#bulkmetadatapath <- "~/Quest Files/bulk_metadata.xlsx"

###### READ IN DATA #############
bulkData <- readRDS(bulk_datapath) # BULK
scData <- readRDS(sc_datapath) # Single Cell
sc_metadata <- read_xlsx(sc_metadatapath) # SC metatdata
bulk_metadata <- read_xlsx(bulk_metadatapath) # Bulk metatdata


######## IMPORTANT FUNCTIONS USED ########

# Normalize a matrix to counts per million (CPM)
Normalize_col_CPM <- function(X){
  return(t(t(X) / colSums(X) * 1e6))
}

# Calculate GT proportions for SC only SC+Bulk
Get_GT_Props <- function(seurat_data, cells_to_remove = c("")){
  # Put cell types and sample IDs in a df
  Cell_labels <- data.frame(cell_type = seurat_data$type_B036, sample_id = seurat_data$orig.ident)
  
  # Create wide df
  cell_label_summary <- Cell_labels %>%
    group_by(sample_id, cell_type) %>%
    summarize(cell_count = n()) %>%
    spread(sample_id, cell_count)
  
  # Replace NAs w/ 0's
  cell_label_summary[is.na(cell_label_summary)] <- 0
  
  # Reshape data to long format
  df_long <- gather(cell_label_summary, key = "sample_id", value = "count", -cell_type)
  
  # Create boxplot - supplemental info for average cell-type counts per sample
  ggplot(df_long, aes(x = cell_type, y = log(count))) +
    geom_boxplot() +
    labs(x = "Cell Type", y = "Count") +
    theme(axis.text.x = element_text(angle = 90))
  
  
  # Decide on cell types to Remove/Keep
  toKeep <- as.character(cell_label_summary$cell_type[!(cell_label_summary$cell_type %in% cells_to_remove)])
  
  
  ### GENERATE Ground Truth data ###
  
  # Obtain summaries for no tumors and for cells to keep
  
  cell_label_summary_TO_USE <- cell_label_summary %>% 
    filter(!cell_type %in% cells_to_remove)
  
  # Divide each value in the data frame by the total number of cells for each column
  cell_proportions_GT <- data.frame(cell_label_summary$cell_type,
                                    apply(cell_label_summary %>% 
                                            select(-cell_type), 2, 
                                          function(x) x/sum(x)))
  
  rownames(cell_proportions_GT) <- cell_proportions_GT$cell_label_summary.cell_type
  cell_proportions_GT <- cell_proportions_GT[,-1]
  
  # THIS IS THE GT
  cell_proportions_GT <- as.matrix(cell_proportions_GT)
  return(cell_proportions_GT)
}

# MAIN DECON FUNCTION
Deconvolution <- function(bulkData_toUse, method, sigMatrix_toUse = "None"){
  if(method == "FARDEEP"){
    # Using FARDEEP
    est_proportions <- t(fardeep(sigMatrix_toUse,bulkData_toUse)$relative.beta)
  }
  else if(method == "EPIC"){
    # Using EPIC
    est_proportions <- t(EPIC(bulkData_toUse, reference = list(refProfiles = sigMatrix_toUse, sigGenes = rownames(sigMatrix_toUse)), withOtherCells = FALSE)$cellFractions)
  }
  else if(method == "xCell"){
    # Using xCell
    # HOW DO I INPUT SIG MATRIX?
    est_proportions <- xCellAnalysis(bulkData_toUse, signatures = , cell.types.use = NULL)
  }
  else if(method == "MCPcounter"){
    # Using MCPcounter
    # HOW DO I INPUT SIG MATRIX?
    est_proportions <- MCPcounter.estimate(bulkData_toUse)
  }
  else if(method == "CIBERSORT"){
    est_proportions <- CIBERSORT_david(sigMatrix_toUse, bulkData_toUse)
  }
  return(est_proportions)
}

# Calculate Deconvolution Accuracy. Returns average cc and rmse, and lists each of the values as a vector.
DeconAcc <- function(GT, Estimate, allResults = FALSE){
  # Define function to calculate RMSE and correlation coefficient
  calc_stats <- function(x, y) {
    rmse <- sqrt(mean((x - y)^2))
    corr <- cor(x, y)
    return(list(RMSE = rmse, Correlation = corr))
  }
  
  # Initialize empty vectors to store RMSE and correlation coefficients
  rmse_vec <- numeric()
  corr_vec <- numeric()
  est_reordered <- Estimate[match(rownames(GT), rownames(Estimate)), ]
  
  # Loop over each sample id
  for (sample_id in colnames(GT)) {
    # Get the data for this sample from both matrices
    sc_GT_props_sample <- GT[,sample_id]
    est_sample <- est_reordered[,sample_id]
    
    # Calculate RMSE and correlation coefficient
    stats <- calc_stats(sc_GT_props_sample, est_sample)
    
    # Add results to the vectors
    rmse_vec <- c(rmse_vec, stats$RMSE)
    corr_vec <- c(corr_vec, stats$Correlation)
  }
  # Calculate and return average RMSE and correlation coefficient
  avg_rmse <- mean(rmse_vec)
  avg_corr <- mean(corr_vec)
  
  if(allResults){
    final <- list(avg_rmse = avg_rmse,
                  avg_corr = avg_corr,
                  rmse_vec = rmse_vec,
                  corr_vec = corr_vec)
  } else{
    final <- list(avg_rmse = avg_rmse,
                  avg_corr = avg_corr)
  }
  
  return(final)
}

# Create a function to create a plot for each deconvolution method
create_plot <- function(data, method) {
  ggplot(data[data$method_type == method,], aes(x = logfc, y = minpct, size = 1/RMSE, color = CC)) +
    geom_point() +
    labs(title = method, x = "logfc", y = "minpct") +
    color_scale +
    scale_size_continuous(limits = c(floor(1/max(accuracy_df$RMSE)), ceiling(1/min(accuracy_df$RMSE))), range = c(1,10)) +
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5), labels = c(0.25, 0.5, 0.75, 1, 1.25, 1.5)) +
    theme_minimal()
}

# PCA Plot given Cell Type prop estimates
PlotPCA <- function(data, toRemove){
  # Run PCA
  pca <- prcomp(t(data[!(rownames(data) %in% toRemove),]))
  print(round(head(pca$sdev^2 / sum(pca$sdev^2),2)*100,1))
  
  # Extract first two components
  pca_scores <- data.frame(Sample = colnames(data), pca$x[, 1:2])
  
  # Create scatter plot
  ggplot(pca_scores, aes(x = PC1, y = PC2)) +
    geom_point() +
    labs(x = "PC1", y = "PC2", title = "PCA of Cell Type Proportions")
  
  
  pca_scores_color <- merge(pca_scores, bulk_only_metadata_filtered, by.x = "Sample", by.y = "Tumor_Sample_Barcode") %>% subset(Response_Binary != "NA")
  ggplot(pca_scores_color, aes(x = PC1, y = PC2, color = factor(Response_Binary))) +
    geom_point() +
    labs(x = "PC1", y = "PC2", title = "PCA of Cell Type Proportions") +
    scale_color_manual(values = c("#00BFC4", "#F8766D", "#c400b4")) +
    guides(color = guide_legend(title = "Drug Efficacy"))
}

# Run predictions and return accuracy stats

Predictions <- function(data, vars, cell_types, train_per, cutoff, pred_method = "RF"){
  # Check input arguments
  if (!is.data.frame(data)) stop("data must be a dataframe")
  if (!all(vars %in% colnames(data))) stop("vars must be valid column names in data")
  if (!all(cell_types %in% colnames(data))) stop("cell_types must be valid column names in data")
  if (train_per <= 0 || train_per >= 1) stop("train_per must be between 0 and 1")
  if (cutoff < 0 || cutoff > 1) stop("cutoff must be between 0 and 1")
  
  vars_to_keep <- c(vars, cell_types)
  # Prepare the data for Lasso regression
  X <- data[, vars_to_keep]
  y <- ifelse(data$Response_Binary == "R", 1,0)
  
  # Split the data into training and testing sets
  set.seed(123)
  train_idx <- sample(nrow(data), train_per * nrow(data))
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Standardize the predictors
  X_train_scaled <- scale(X_train)
  X_test_scaled <- scale(X_test)
  
  if(pred_method == "LASSO"){
    # Fit a Lasso model using the penalized package
    lasso_model <- cv.glmnet(x = as.matrix(X_train_scaled), y = y_train, alpha = 1)
    
    # Get the predicted probabilities and predicted classes for the test set
    predicted_prob <- predict(lasso_model, newx = as.matrix(X_test_scaled), type = "response")
    predicted_class <- ifelse(predicted_prob > cutoff, 1, 0)
  }
  
  if(pred_method == "RF"){
    # Fit a Random Forest model
    rf_model <- randomForest(x = X_train_scaled, y = y_train)
    
    # Get the predicted probabilities and predicted classes for the test set
    predicted_prob <- predict(rf_model, X_test_scaled)
    predicted_class <- ifelse(predicted_prob > cutoff, 1, 0)
  }
  
  # Calculate the AUC and prediction accuracy
  auc <- roc(y_test, predicted_prob)$auc
  accuracy <- mean(predicted_class == y_test)
  return(list(Accuracy = accuracy, AUC = auc))
}


############ PART 1: Create Table 1 ###################

# Current bulk vars of interest
voi <- c("Tumor_Sample_Barcode","Patient_name","Viral_Status","Age","Sex","Tissue","Response_post_bx","Response_Binary",
              "Primary_vs_Met","OS_days","IO_OS_days",
              "Timepoint","IO_Naive")

# Select variables of interest and filter for Timepoint = 1 (pre-IO)
## Also remove patient 68 and 82 LN, since there are 2 and we only have skin for the matched bulk samples. 
sc_bulk_metadata_filtered <- sc_metadata %>% subset(Patient_name %in% intersect(sc_metadata$Patient_name, bulk_metadata$Patient_name)) %>%
  subset(Tumor_Sample_Barcode != "MO_MCC_068LN" & Tumor_Sample_Barcode != "MO_MCC_082LN") %>% 
  select(voi) %>% subset(Timepoint == 1) %>% subset(is.na(IO_Naive)|IO_Naive == "yes")

sc_only_metadata_filtered <- sc_metadata %>% subset(!(Patient_name %in% intersect(sc_metadata$Patient_name, bulk_metadata$Patient_name))) %>%
  subset(Tumor_Sample_Barcode != "MO_MCC_068LN" & Tumor_Sample_Barcode != "MO_MCC_082LN") %>% 
  select(voi) %>% subset(Timepoint == 1) %>% subset(is.na(IO_Naive)|IO_Naive == "yes")

bulk_only_metadata_filtered <- bulk_metadata %>% subset(!(Patient_name %in% intersect(sc_metadata$Patient_name, bulk_metadata$Patient_name))) %>%
  select(voi) %>% subset(Timepoint == 1) %>% subset(IO_Naive == TRUE)

bulk_sc_metadata_filtered <- bulk_metadata %>% subset((Patient_name %in% intersect(sc_metadata$Patient_name, bulk_metadata$Patient_name))) %>%
  select(voi) %>% subset(Timepoint == 1) %>% subset(IO_Naive == TRUE) %>% subset(Patient_name %in% sc_bulk_metadata_filtered$Patient_name)

# Grab bulk barcode labels for deconvolution later
bulk_only_labels <- bulk_only_metadata_filtered$Tumor_Sample_Barcode
sc_only_labels <- sc_only_metadata_filtered$Tumor_Sample_Barcode

sc_bulk_labels <- sc_bulk_metadata_filtered$Tumor_Sample_Barcode
bulk_sc_labels <- bulk_sc_metadata_filtered$Tumor_Sample_Barcode

sc_IO_labels <- (sc_metadata %>% subset(Tumor_Sample_Barcode != "MO_MCC_068LN" & Tumor_Sample_Barcode != "MO_MCC_082LN") %>% 
                   subset(Timepoint == 1) %>% subset(is.na(IO_Naive)|IO_Naive == "yes"))$Tumor_Sample_Barcode

# Select variables to keep prior to merging
toDrop <- c("Tumor_Sample_Barcode","Timepoint","IO_Naive")

# Select these variables and ensure numeric variables
scToMerge <- sc_only_metadata_filtered %>% select(-toDrop) %>%
  mutate(Age = as.numeric(Age),
         OS_days = as.numeric(OS_days),
         IO_OS_days = as.numeric(ifelse(IO_OS_days == "NA", NA, IO_OS_days)))

sc_bulk_ToMerge <- sc_bulk_metadata_filtered %>% select(-toDrop) %>%
  mutate(Age = as.numeric(Age),
         OS_days = as.numeric(OS_days),
         IO_OS_days = as.numeric(ifelse(IO_OS_days == "NA", NA, IO_OS_days)))

bulkToMerge <- bulk_only_metadata_filtered %>% select(-toDrop) %>% 
  mutate(Age = as.numeric(Age),
         OS_days = as.numeric(OS_days),
         IO_OS_days = as.numeric(ifelse(IO_OS_days == "NA", NA, IO_OS_days)))

# Merge Data
merged_data <- bind_rows(
  scToMerge %>% mutate(Type = "SC Only"),
  sc_bulk_ToMerge %>% mutate(Type = "SC+Bulk"),
  bulkToMerge %>% mutate(Type = "Bulk Only")
)

merged_data_forTab <- merged_data %>% mutate(Viral_Status = ifelse(Viral_Status == "VP","MCPyV+",
                                                               ifelse(Viral_Status == "VN","MCPyV-",NA)),
                                             Sex = ifelse(Sex == "M", "Male", "Female"),
                                             Tissue = ifelse(Tissue == "LN","Lymph Node",
                                                             ifelse(Tissue == "Skin","Skin","Other")),
                                             Response_post_bx = case_when(
                                               Response_post_bx == "CR" ~ "Complete Response",
                                               Response_post_bx == "NA" ~ "No IO treatment",
                                               Response_post_bx == "PD" ~ "Disease Progression",
                                               Response_post_bx == "PR" ~ "Partial Response",
                                               Response_post_bx == "SD" ~ "Stable Disease",
                                               TRUE ~ NA_character_),
                                             Response_Binary = case_when(
                                               Response_Binary == "R" ~ "Response",
                                               Response_Binary == "NA" ~ "No IO treatment",
                                               Response_Binary == "NR" ~ "No Response",
                                               TRUE ~ NA_character_)) %>%
                                        mutate(Response_post_bx = factor(Response_post_bx,
                                                                         levels = c("Disease Progression",
                                                                                    "Stable Disease",
                                                                                    "Partial Response",
                                                                                    "Complete Response")),
                                               Response_Binary = factor(Response_Binary,
                                                                         levels = c("No Response",
                                                                                    "Response")),
                                               Tissue = factor(Tissue,
                                                                        levels = c("Skin",
                                                                                   "Lymph Node",
                                                                                   "Other")))

#### GENERATE TABLE 1 #####

vars_for_tab <- c("Viral_Status","Age","Sex",
                  "Tissue","Response_post_bx","Response_Binary","Primary_vs_Met",
                  "OS_days","IO_OS_days","Type")
stratifier <- "Type"
normalVars <- "Age"
nonNormalVars <- c("OS_days", "IO_OS_days")
varsToFactor <- vars_for_tab[!vars_for_tab %in% c(stratifier,normalVars,nonNormalVars)]

# Modify the order of variables in the table
merged_data_forTab <- merged_data_forTab %>%
  select(Age, Sex, Tissue, Viral_Status, Primary_vs_Met, Response_post_bx, Response_Binary, OS_days, IO_OS_days, Type)

tbl_summary_1 <-
  merged_data_forTab %>% 
    tbl_summary(
      by = stratifier,
    label = list(Age ~ "Age, years", Tissue ~ "Location of Biopsy", Viral_Status ~ "Viral Status", Primary_vs_Met ~ "Metastacicity",
                   Response_post_bx ~ "Response to IO", Response_Binary ~ "Binary Response to IO",
                   OS_days ~ "Survival Time following biopsy, days", IO_OS_days ~ "Survival Time following IO, days"),
      statistic = list(normalVars ~ "{mean} ({sd})",
                       nonNormalVars ~ "{median} ({p25}, {p75})"),
                       #varsToFactor ~ "{n} ({p}%)"),
      digits = list(normalVars ~ c(0,1),
                    nonNormalVars ~ c(0, 0)),
      missing = "ifany",
    missing_text = "No IO treatment"
    ) %>%
    add_overall() %>% bold_labels()


# Display Table 1
tbl_summary_1

# Export the table to a Word document
doc <- read_docx()
doc <- body_add_flextable(doc, as_flex_table(tbl_summary_1), align = "left")
print(doc, target = "Manuscript Figures/Table1.docx")




############ PART 2: Generate pseudo-bulk data, GT props, make UMAPs ###################

# Identify cells to remove
levels_to_remove <- c("Mast", "Melano")

# Identify rows to remove
scData$removeCells <- scData$type_B036 %in% levels_to_remove
# Remove rows
scData <- subset(scData, removeCells == FALSE)
# Remove levels
scData$type_B036 <- factor(scData$type_B036, levels = setdiff(levels(scData$type_B036), levels_to_remove))

# Filter scData for SC Only and SC+Bulk
sc_only_Data_filt <- subset(scData, subset = orig.ident %in% sc_only_labels)
sc_bulk_Data_filt <- subset(scData, subset = orig.ident %in% sc_bulk_labels)

# Make UMAPs for those who underwent IO treatment
sc_IO_Data_filt <- subset(scData, orig.ident %in% sc_IO_labels) %>% 
  subset(Responder_post_bx == "yes"|Responder_post_bx == "no")

# Delete Massive SC dataset
rm(scData)

sc_IO_Data_filt <- NormalizeData(sc_IO_Data_filt)
sc_IO_Data_filt <- FindVariableFeatures(sc_IO_Data_filt)
sc_IO_Data_filt <- ScaleData(sc_IO_Data_filt)
sc_IO_Data_filt <- RunPCA(sc_IO_Data_filt)
#ElbowPlot(sc_IO_Data_filt)
sc_IO_Data_filt <- RunUMAP(sc_IO_Data_filt, dims = 1:20)
unique(sc_IO_Data_filt$orig.ident)

# Create a vector of sample labels for display
sample_labels <- paste0("Sample ", 1:21)

# Convert orig.ident to a factor and relabel the levels for display
sc_IO_Data_filt$orig.ident <- factor(sc_IO_Data_filt$orig.ident)
levels(sc_IO_Data_filt$orig.ident) <- sample_labels

ID_UMAP <- DimPlot(sc_IO_Data_filt, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample ID")
Cell_Type_UMAP <- DimPlot(sc_IO_Data_filt, reduction = "umap", group.by = "type_B036") +
  ggtitle("Cell Types")
EFF_UMAP <- DimPlot(sc_IO_Data_filt, reduction = "umap", group.by = "Responder_post_bx") +
  ggtitle("IO efficacy") +
  guides(color = guide_legend(direction = "horizontal", override.aes = list(size = 3)))
#DimPlot(sc_IO_Data_filt, reduction = "umap", group.by = "Response_post_bx")

# Arrange the plots side by side - FIGURE 1
gridExtra::grid.arrange(ID_UMAP, Cell_Type_UMAP, EFF_UMAP, ncol = 3)

### MAKE PSUEDO BULK FOR TESTING ###
sc_only_filtAgg <- AggregateExpression(sc_only_Data_filt, group.by = "orig.ident", slot = "count")$RNA
sc_bulk_filtAgg <- AggregateExpression(sc_bulk_Data_filt, group.by = "orig.ident", slot = "count")$RNA

# THIS IS THE PSEUDO BULK DATA for each group
pseudo_bulk_sc_only <- Normalize_col_CPM(sc_only_filtAgg)
pseudo_bulk_sc_bulk <- Normalize_col_CPM(sc_bulk_filtAgg)

# Get GT Props - USE THESE FOR TESTING!!
sc_only_GT_props <- Get_GT_Props(sc_only_Data_filt)
sc_bulk_GT_props <- Get_GT_Props(sc_bulk_Data_filt)

# Save the Seurat object as an RDS file
saveRDS(scData_filt, file = "/projects/b1036/davidhparkinson/Bulk+SC Info/scData_filtered.rds")

####### PART 3 - BUILDING SIGNATURE MATRIX AND OPTIMIZING GENE MARKER SELECTION ################

# Aggregate sc data by cell-type within sample ID
sig_info_AGG <- AggregateExpression(sc_only_Data_filt, group.by = c("type_B036", "orig.ident"), slot = "counts")$RNA

# Change column names
new_colnames <- gsub("_.*", "", colnames(sig_info_AGG))
colnames(sig_info_AGG) <- new_colnames

# Average values by column name
sig_info_AGG_df <- as.data.frame(sig_info_AGG)
result_df <- data.frame(matrix(nrow = nrow(sig_info_AGG_df), ncol = 0))

# Loop through unique column names
for (col_name in unique(colnames(sig_info_AGG_df))) {
  # Calculate the mean for columns with the same name
  mean_col <- rowMeans(sig_info_AGG_df[, colnames(sig_info_AGG_df) == col_name, drop = FALSE])
  
  # Add the mean column to the result data.frame
  result_df[[col_name]] <- mean_col
}

###### Convert data.frame back to matrix - this is the basis for the sig matrix!!!!!!!!!
averaged_matrix <- as.matrix(result_df)
rownames(averaged_matrix) <- rownames(sig_info_AGG)


######### FIND MARKER GENES #############
# Define the desired values for logfc.threshold and min.pct
logfc.thresholds <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5)
min.pcts <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# Initialize an empty list to store the marker gene matrices
marker_genes_list <- list()
marker_genes_info_list <- list()
# Loop over the combinations of logfc.threshold and min.pct values
for (logfc.threshold in logfc.thresholds) {
  for (min.pct in min.pcts) {
    # Call Seurat::FindAllMarkers() with the specified parameters
    marker_genes <- Seurat::FindAllMarkers(sc_only_Data_filt, slot = "data",
                                           logfc.threshold = logfc.threshold,
                                           min.pct = min.pct)
    # Add the resulting df to the list
    marker_genes_info_list[[paste0("logfc.", logfc.threshold, ".minpct.", min.pct)]] <- marker_genes
    marker_genes_list[[paste0("logfc.", logfc.threshold, ".minpct.", min.pct)]] <- unique(marker_genes$gene)
  }
}

# Save list as JSON file
write_json(marker_genes_list, "marker_genes.json")


####### PART 4 - BUILDING SIGNATURE MATRIX AND OPTIMIZING GENE MARKER SELECTION ################

##### GET ACCURACY RESULTS FOR ALL METHODS + MARKER GENE ALGORITHMS #####
decon_methods <- c("EPIC","CIBERSORT","FARDEEP")
accuracy_info <- list()
i <- 0
for(marker_genes in marker_genes_list){
  i <- i + 1
  sig_mat <- averaged_matrix[rownames(averaged_matrix) %in% marker_genes,]
  for(method in decon_methods){
    est <- Deconvolution(pseudo_bulk_sc_only, method, sig_mat)
    accuracy_info[[paste0(method, "_", names(marker_genes_list)[i])]] <- DeconAcc(sc_only_GT_props,est, TRUE)
    print(paste0(names(marker_genes_list)[i]," version of the ", method, " method is COMPLETE."))
  }
}

# Unlist the nested lists and create a data frame
accuracy_df <- do.call(rbind, lapply(accuracy_info, function(x) data.frame(CC = x$avg_corr, RMSE = x$avg_rmse)))

# Extract the method, logfc, and minpct information
accuracy_df$method_type <- sapply(rownames(accuracy_df), function(x) gsub("_(.*)", "", x))
accuracy_df$logfc <- as.numeric(sapply(rownames(accuracy_df), function(x) gsub("^.*logfc\\.([0-9\\.]+)\\.minpct.*$", "\\1", x)))
accuracy_df$minpct <- as.numeric(gsub(".*\\.minpct\\.[0-9]+", "\\1", rownames(accuracy_df)))

# Rank each method by CC and RMSE
accuracy_df$Rank <- rank(-accuracy_df$CC + accuracy_df$RMSE, ties.method = "min")

# Save accuracy info to be saved!
write_csv(accuracy_df, "accuracy_info.csv")

###### BEST GENE LIST #######
# Best method details
best_option <- list(logfc = (accuracy_df %>% subset(Rank == 1))$logfc,
                    minpct = (accuracy_df %>% subset(Rank == 1))$minpct,
                    CC = (accuracy_df %>% subset(Rank == 1))$CC,
                    RMSE = (accuracy_df %>% subset(Rank == 1))$RMSE)

print("Best options stats:")
print(best_option)
best_marker_genes
# Marker Gene List
best_marker_genes <- marker_genes_list[[paste0("logfc.",best_option$logfc,".minpct.",best_option$minpct)]]

# SIGNATURE MATRIX TO USE!
SIG_MATRIX_FINAL <- averaged_matrix[rownames(averaged_matrix) %in% best_marker_genes,]

## FIGURE 3 ##
heatmap(SIG_MATRIX_FINAL,
        #xlab = "Cell Type",
        ylab = "Gene Expression Profile")









# Read in the marker genes list in case you need it!
marker_genes_list <- fromJSON("marker_genes.json")

# Read in the accuracy info in case you need it!
accuracy_df <- read_csv("accuracy_info.csv")

####### MAKE PLOTS ##########

# Define the color scale based on CC values
color_scale <- scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1))

# Create the three plots
plot_epic <- create_plot(accuracy_df, "EPIC")
plot_cibersort <- create_plot(accuracy_df, "CIBERSORT")
plot_fardeep <- create_plot(accuracy_df, "FARDEEP")

##### Arrange the plots side by side - FIGURE 2A #####
gridExtra::grid.arrange(plot_epic, plot_cibersort, plot_fardeep, ncol = 3)



####### THIS CAN BE FIXED FOR A POTENTIAL FIGURE 2B LATER ########## ##############
# Filter lists with names starting with "EPIC_logfc.0.75.minpct."
filtered_accuracy_info <- accuracy_df %>% subset(method_type == "EPIC", logfc = 0.75)

# Function to extract minpct from name
get_minpct <- function(name) {
  as.numeric(str_extract(name, "(?<=EPIC_logfc.0.75.minpct.)\\d+(\\.\\d+)?"))
}

###### Create dataframe
#result_df <- filtered_accuracy_info %>%
#  enframe(name = "minpct") %>%
#  mutate(minpct = map_dbl(minpct, get_minpct),
#         value = map(value, ~tibble(RMSE = .x$rmse_vec, CC = .x$corr_vec))) %>%
#  unnest(value)


###### Create density plot - FIGURE 2B
#ggplot(result_df, aes(x = CC, y = 1/RMSE, color = as.factor(minpct))) +
#  geom_point() +
#  geom_density_2d() +
#  labs(x = "CC", y = "1/RMSE", color = "minpct") +
#  theme_minimal()

################### END FIX ###############################################



####### PART 4 - DECONVOLUTE SC+BULK DATA, BIAS ANALYSIS ################

# Clean bulk data - fix gene IDs
bulk_filtered <- (bulkData[6:nrow(bulkData),])
rownames(bulk_filtered) <- bulk_filtered$ENSEMBL_GeneID
bulk_filtered <- as.matrix(bulk_filtered %>% select(-1))

### Rename the gene IDs ###
gene_ann <- read_csv("/projects/b1036/Jasen/code/experiment.code/ce.0096.updateGeneAnno.R/gene_anno.csv")

new_row_names <- ifelse(rownames(bulk_filtered) %in% gene_ann$gene_id,
                        gene_ann$symbol[match(rownames(bulk_filtered), gene_ann$gene_id)],
                        rownames(bulk_filtered))

rownames(bulk_filtered) <- new_row_names

# Separate bulk data into SC+Bulk and Bulk only
bulk_sc_data <- bulk_filtered[,colnames(bulk_filtered) %in% bulk_sc_labels]
bulk_only_data <- bulk_filtered[,colnames(bulk_filtered) %in% bulk_only_labels]

# Normalize to CPM prior to deconvolution
bulk_sc_forDecon <- Normalize_col_CPM(bulk_sc_data)
bulk_only_forDecon <- Normalize_col_CPM(bulk_only_data)

# Deconvolute SC+Bulk and Bulk only using the optimal marker genes
bulk_sc_est <- Deconvolution(bulk_sc_forDecon, "EPIC", SIG_MATRIX_FINAL)
bulk_only_est <- Deconvolution(bulk_only_forDecon, "EPIC", SIG_MATRIX_FINAL)


###### SC+BULK Bias #########
# Adjust colnames for agreement
colnames(bulk_sc_est) <- colnames(sc_bulk_GT_props)

# Calculate Error for each cell-type
df_results <- data.frame(row.names = colnames(bulk_sc_est))
for(i in 1:nrow(bulk_sc_est)){
  err <- (bulk_sc_est[i,]-sc_bulk_GT_props[i,])/sc_bulk_GT_props[i,]
  df_results[[rownames(bulk_sc_est)[i]]] <- err
}

# Create box plots for each column using ggplot - FIGURE 4
ggplot(gather((df_results)), aes(x = key, y = value)) + 
  geom_boxplot() + 
  labs(x = "", y = "Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  # -1 means the estimate was 0. Inf means GT was 0.
  ylim(c(-1,1000)) +
  ylab("Percent change from SC proportions") +
  xlab("Cell Type")

# Average Percent errors
print(paste0("Endo error: ",mean(df_results$Endo[is.finite(df_results$Endo)])))
print(paste0("Fibro error: ",mean(df_results$Fibro[is.finite(df_results$Fibro)])))
print(paste0("Plasma error: ",mean(df_results$Plasma[is.finite(df_results$Plasma)])))
print(paste0("Tumor cycling error: ",mean(df_results$`Tumor cycling`[is.finite(df_results$`Tumor cycling`)])))


####### PART 4 - DECONVOLUTE BULK ONLY DATA, VISUAL ANALYSIS ################

# PCA plot of all bulk data with IO treatment
PlotPCA(bulk_only_est, c("Fibro", "Endo", "Plasma", "Tumor cycling", "Tumor"))

# Convert bulk_only estimates matrix to long format
proportions_long <- reshape2::melt(bulk_only_est, varnames = c("Cell_Type", "Sample"))

# Join with metadata and create boxplot with color
proportions_color <- merge(proportions_long, bulk_only_metadata_filtered, by.x = "Sample", by.y = "Tumor_Sample_Barcode")
ggplot(proportions_color %>% subset(Response_Binary != "NA") %>% 
         mutate(`IO Efficacy` = 
                  ifelse(Response_Binary == "NR","No","Yes"))
       , aes(x = Cell_Type, y = value, fill = `IO Efficacy`)) +
  geom_boxplot(color = "black", alpha = 0.8, outlier.color = "red", outlier.alpha = 0.8) +
  facet_wrap(~ `IO Efficacy`, ncol = 1) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "#c400b4")) +
  labs(x = "Cell Type", y = "Proportion", title = "Boxplot of Cell Type Proportions by Drug Efficacy") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Run t-tests for each cell-type
for(i in unique(proportions_color$Cell_Type)){
  t_res <- wilcox.test((proportions_color %>% subset(Response_Binary == "R") %>% subset(Cell_Type == i))$value,
                  (proportions_color %>% subset(Response_Binary == "NR") %>% subset(Cell_Type == i))$value)
  print(paste0(i,": ", round(t_res$p.value,4)))
}

# Fit a logistic regression model, check for significance, print p-values
for(i in rownames(bulk_only_est)){
  mod <- paste0("Response ~ Sex + Age + Primary_vs_Met + `",i,"`")
  model <- glm(mod, data = bulk_only_forLASSO %>% mutate(Response = ifelse(Response_Binary == "R",1,0)), family = "binomial")
  print(paste0(i,": ",summary(model)$coefficients[5,4]))
}

# Plot the log of Cell Type proportions stratified by IO efficacy
ggplot(proportions_color %>% subset(Response_Binary != "NA") %>% 
         mutate(`IO Efficacy` = 
                  ifelse(Response_Binary == "NR","No","Yes")), 
       aes(x = Cell_Type, y = log(value), color = `IO Efficacy`)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Log of the proportion") +
  xlab("Cell Type")


# RUN PREDICTIONS WITH CLINICAL VARIABLES, THEN ADD IN CELL TYPE PROPS #

# Transpose the bulk_only_est dataframe and convert it to a tibble, log normalize
bulk_only_est_t <- as_tibble(t(log(bulk_only_est)))

# Rename the column names to match the sample IDs. Merge with metadata
bulk_only_est_t$Tumor_Sample_Barcode <- colnames(bulk_only_est)
bulk_only_toPred <- merge(bulk_only_metadata_filtered, bulk_only_est_t, by = "Tumor_Sample_Barcode")

# Prep for LASSO and Random Forest, by converting vars to numeric
bulk_only_forPred_FINAL <- bulk_only_toPred %>% subset(Response_Binary != "NA") %>% 
  dplyr::mutate(Primary_vs_Met = ifelse(Primary_vs_Met == "Metastasis",1,0),
                Viral_Status = ifelse(Viral_Status == "VP",1,0),
                Sex = ifelse(Sex == "M",1,0),
                Age = as.numeric(Age))


# Insert biased cell types
bad_cells <- c("Fibro", "Endo","Plasma", "Tumor cycling","Tumor")
#bad_cells <- c("Fibro")
#good_cells <- rownames(bulk_only_est)[!(rownames(bulk_only_est) %in% bad_cells)]

# Initialize an empty list to store the results
comb_list <- list()

# Loop through all possible vector lengths from 1 to 5
for (i in 1:length(bad_cells)) {
  # Generate all possible combinations of length i
  combs <- combn(bad_cells, i)
  for(j in 1:ncol(combs)){
    comb_list[[paste(combs[,j], collapse = ", ")]] <- Predictions(bulk_only_forLASSO, c("Viral_Status", "Age", "Sex"), rownames(bulk_only_est)[!(rownames(bulk_only_est) %in% combs[,j])], 0.75, 0.75, "RF")
  }
}

pred_acc_results <- data.frame(Removed_cells = names(comb_list),
                 AUC = unlist(lapply(comb_list, function(x) x$AUC)),
                 accuracy = unlist(lapply(comb_list, function(x) x$Accuracy)))
write_csv(pred_acc_results,"/projects/b1036/davidhparkinson/Manuscript Figures/Table2.csv")

bad_immune_cells <- c("Fibro", "Endo","Plasma")
# Without cells
Predictions(bulk_only_forLASSO, c("Viral_Status", "Age", "Sex"), c(), 0.75, 0.75, "RF")

# Immune cells ONLY
Predictions(bulk_only_forLASSO, c("Viral_Status", "Age", "Sex"), rownames(bulk_only_est)[!(rownames(bulk_only_est) %in% c("Tumor cycling","Tumor"))], 0.75, 0.75, "RF")

# With ALL cells
Predictions(bulk_only_forLASSO, c("Viral_Status", "Age", "Sex"), rownames(bulk_only_est), 0.75, 0.75, "RF")

# Without any biased cells 
Predictions(bulk_only_forLASSO, c("Viral_Status", "Age", "Sex"), rownames(bulk_only_est)[!(rownames(bulk_only_est) %in% c(bad_immune_cells,c("Tumor cycling","Tumor")))], 0.75, 0.75, "RF")

# Immune cells ONLY without biased cells 
Predictions(bulk_only_forLASSO, c("Viral_Status", "Age", "Sex"), rownames(bulk_only_est)[!(rownames(bulk_only_est) %in% bad_immune_cells)], 0.75, 0.75, "RF")

########## DONE ############