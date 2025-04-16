# Title: Microarray Data Analysis with GEOquery and limma
pacman::p_load(
    GEOquery, tidyverse, ggrepel, limma, oligo, DT, pheatmap, tidyplots, affy, oligoClasses, testit
)
# function
isLog2Transformed <- function(data) {
    qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    shouldBeLogged <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
    return(!shouldBeLogged)
}

# First obtain metadata
id <- "GSE19080"
# Check if the GEO ID is valid
# creating a temp folder
if (!dir.exists(paste0(".temp/", id))) {
    dir.create(paste0(".temp/", id))
} else {
    message("Folder already exists!")
}

# Obtaining Metadata
meta <- getGEO(id, GSEMatrix = TRUE, destdir = ".temp")
meta <- meta[[1]]
fData(meta)
varLabels(meta)
head(pData(meta))
colnames(exprs(meta))
dim(meta)

# GPL4133
pd
pd <- pData(meta) |>
    tibble::rownames_to_column("ID") |>
    dplyr::select("characteristics_ch1", "characteristics_ch1.3", ID, supplementary_file) |>
    dplyr::mutate(file = str_split(supplementary_file, "/") |> map_chr(tail, 1)) |>
    dplyr::rename(atl_subtype = "characteristics_ch1", gender = "characteristics_ch1.3") |>
    dplyr::mutate(across(
        c(atl_subtype, gender),
        ~ stringr::str_remove(., "^.*: ")
    )) |>
    mutate(atl_subtype2 = case_when(
        str_detect(atl_subtype, "ATL") ~ "ATL",
        str_detect(atl_subtype, "HAM") ~ "HAMTSP",
        str_detect(atl_subtype, "AC") ~ "AC",
        str_detect(atl_subtype, "Healthy") ~ "HD",
        .default = NA
    ))

## A loop to download all the files in the .temp of this environment
for (i in 1:length(pd7$supplementary_file)) {
    url <- pd7$supplementary_file[i]
    destfile <- file.path(paste0(".temp/", id, "/", pd7$file[i]))
    # Download the file
    tryCatch(
        {
            download.file(url, destfile, mode = "wb", method = "curl")
        },
        error = function(e) {
            # Fallback to default method if curl fails
            download.file(url, destfile, mode = "wb")
        }
    )
    # Optional: Extract if it's a tar file
    if (grepl("\\.tar(\\.gz)?$", destfile)) {
        untar(destfile, exdir = paste0(".temp/", id))
    }
}

# exploring one file
con <- gzfile(file.path("/Users/denriquez/Documents/GitHub/microarray_tutorial/.temp/GSE19080/GSM472372_HISH0553.txt.gz"))
file_lines <- readLines(con, n=1000)
close(con)

### Re-read the raw files
pd7 = pd |>
    filter(ID%in%c(paste0("GSM4723", c(56:73, 82:93))))

agilent_data1 <- read.maimages(
    files = file.path(".temp", id, pd7$file),
    source = "quantarray",
    green.only = FALSE,
    names = pd7$ID,
    other.columns = list(
        Flag = "Ignore Filter"))

head(agilent_data1)

#################
pd7_2 = pd |>
    filter(ID%in%c(paste0("GSM4723", 74:81)))
pd7_2
agilent_data2 <- read.maimages(
    files = file.path(".temp", id, pd7_2$file),
    source = "genepix",
    green.only = FALSE,
    names = pd7_2$ID)

################

#####Convert targets to tibble for easy joining
agilent_data1$targets
agilent_data1$targets <- agilent_data1$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(pd7, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")

#### gene Annotation
gpl2 <- getGEO("GPL9686")
dim(fData(meta)[,c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")])

annot <- Table(gpl2)[, c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")] 

# Add gene symbols to agilent_data1
agilent_data1$genes$GENE_SYMBOL <- annot$SYMBOL[match(agilent_data1$genes$Name, annot$GB_ACC)]

agilent_data1$genes$is_control<-is.na(agilent_data1$genes$SYMBOL)
table((agilent_data1$genes$is_control))
head(agilent_data1$genes)
class(agilent_data1)

# (A) Use negative controls for background correction
agilent_data1_bg <- limma::backgroundCorrect(agilent_data1, method="normexp", offset=20, negctrl=neg_controls) 
 # Normalization
agilent_data1_bg <- limma::normalizeWithinArrays(agilent_data1_bg, method="loess")
agilent_data1_bg <- limma::normalizeBetweenArrays(agilent_data1_bg, method="quantile")
head(agilent_data1_bg) 

#####
plot(agilent_data1_bg$R, agilent_data1_bg$G, 
     xlab="Cy5 (R)", ylab="Cy3 (G)", 
     main="Raw Negative Controls")
boxplot(agilent_data1_bg$M[agilent_data1_bg$genes$is_control == TRUE, ], 
        main="Normalized Negative Controls")

###########################################
# Step 2: Average gene probes (excluding controls)

#gene_avg <- avereps(gene_data, ID=gene_data$genes$GENE_SYMBOL)


# Collapse technical replicates (same probe ID)
agilent_avg <- limma::avereps(agilent_data1_bg, ID=agilent_data1_bg$genes$Name)
gene_data <- agilent_avg[!agilent_avg$genes$is_control, ]
dim(gene_data)
# For single-channel EList (it is not obligatory)
eset <- ExpressionSet(
    assayData = agilent_data$E,
    phenoData = AnnotatedDataFrame(agilent_data$targets),
    featureData = AnnotatedDataFrame(agilent_data$genes)
)

agilent_two_color_qc <- function(agilent_data,
                                min_valid_probes = 0.2,
                                max_probe_na = 0.1,
                                max_sample_na = 0.1,
                                intensity_cutoff_quantile = 0.2,
                                min_samples_above_cutoff = 0.2,
                                verbose = TRUE) {
    # 1. Input Validation --------------------------------------------------------
    required_components <- c("R", "G", "genes", "targets")
    missing_comps <- setdiff(required_components, names(agilent_data))
    if(length(missing_comps)) {
        stop("Missing required components: ", paste(missing_comps, collapse=", "))
    }
    
    # 2. Initial Probe Validation -----------------------------------------------
    if(verbose) message("\n=== INITIAL PROBE VALIDATION ===")
    
    # 2. A - Calculate valid probes (non-NA in both channels)
    is_valid <- !is.na(agilent_data$R) & !is.na(agilent_data$G)
    valid_prop <- rowMeans(is_valid)
    
    # 2. B - Identify probes to keep
    keep_probes <- valid_prop >= min_valid_probes
    
    if(verbose) {
        message("Probe validity summary:")
        message("- Total probes: ", length(keep_probes))
        message("- Valid probes (≥", min_valid_probes*100, "% good values): ", sum(keep_probes))
        message("- Invalid probes removed: ", sum(!keep_probes))
        
        # Show distribution of valid values per probe
        print(quantile(valid_prop, probs = seq(0, 1, 0.1)))
    }

    # 2. C - Filter all components
    agilent_data <- agilent_data[keep_probes, ]
    agilent_data$genes <- agilent_data$genes[keep_probes, ]

    # 3. Robust Data Conversion ---------------------------------------------------
    clean_channel <- function(x, channel_name) {
        # Convert to matrix if needed
        if(!is.matrix(x)) x <- as.matrix(x)
        
        # Identify non-numeric values
        problematic <- !matrix(as.numeric(grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x)), 
                             nrow=nrow(x), ncol=ncol(x))
        
        if(any(problematic)) {
            if(verbose) {
                message("Found ", sum(problematic), " non-numeric values in ", channel_name)
                message("Sample problematic values: ", 
                       paste(unique(as.vector(x)[problematic][1:5]), collapse=", "), "...")
            }
            
            # Convert to numeric, forcing NAs for non-convertible values
            x <- matrix(as.numeric(x), nrow=nrow(x), ncol=ncol(x))
            
            # Replace NAs with channel median
            x[is.na(x)] <- median(x, na.rm=TRUE)
        } else {
            storage.mode(x) <- "numeric"
        }
        x
    }
    
    agilent_data$R <- clean_channel(agilent_data$R, "Red channel")
    agilent_data$G <- clean_channel(agilent_data$G, "Green channel")
    
    # Handle background channels if they exist
    if("Rb" %in% names(agilent_data)) {
        agilent_data$Rb <- clean_channel(agilent_data$Rb, "Red background")
    }
    if("Gb" %in% names(agilent_data)) {
        agilent_data$Gb <- clean_channel(agilent_data$Gb, "Green background")
    }
    
    # 3. Control Probe Handling -------------------------------------------------
    if(verbose) message("\nIdentifying control probes...")
    control_types <- c("^SSC", "^Arabidopsis")
    is_control <- Reduce(`|`, lapply(control_types, function(p) grepl(p, agilent_data$genes$Name)))
    
    if(verbose) {
        message("Control probe summary:")
        print(table(gsub("_.*", "", agilent_data$genes$Name[is_control])))
    }
    
    # 4. Dimensional Checks ----------------------------------------------------
    if(verbose) message("\nVerifying dimensions...")
    dim_checks <- list(
        list(nrow(agilent_data$R), nrow(agilent_data$G), "Channel row mismatch"),
        list(ncol(agilent_data$R), ncol(agilent_data$G), "Channel column mismatch"),
        list(nrow(agilent_data$R), nrow(agilent_data$genes), "Gene annotation row mismatch"),
        list(ncol(agilent_data$R), nrow(agilent_data$targets), "Sample annotation mismatch")
    )
    
    for(chk in dim_checks) {
        if(chk[[1]] != chk[[2]]) stop(chk[[3]])
    }
    
    # 5. Missing Value Analysis ------------------------------------------------
    if(verbose) message("\nAnalyzing missing values...")
    na_stats <- list(
        R = list(
            probes = rowMeans(is.na(agilent_data$R)),
            samples = colMeans(is.na(agilent_data$R))
        ),
        G = list(
            probes = rowMeans(is.na(agilent_data$G)),
            samples = colMeans(is.na(agilent_data$G))
        )
    )
    
    # Identify completely failed probes/samples
    complete_failures <- list(
        probes = which(rowSums(is.na(agilent_data$R)) == ncol(agilent_data$R) | 
                       rowSums(is.na(agilent_data$G)) == ncol(agilent_data$G)),
        samples = which(colSums(is.na(agilent_data$R)) == nrow(agilent_data$R) &
                       colSums(is.na(agilent_data$G)) == nrow(agilent_data$G))
    )
    
    # 6. Data Cleaning ---------------------------------------------------------
    if(length(complete_failures$probes) > 0) {
        if(verbose) message("\nRemoving ", length(complete_failures$probes), " completely failed probes")
        agilent_data <- agilent_data[-complete_failures$probes, ]
    }
    
    # 7. Quality Control Processing --------------------------------------------
    if(verbose) message("\nRunning QC pipeline...")
    
    # Background correction
    agilent_data <- limma::backgroundCorrect(agilent_data, method="normexp", offset=1)
    
    # Normalization
    agilent_data <- limma::normalizeWithinArrays(agilent_data, method="loess")
    agilent_data <- limma::normalizeBetweenArrays(agilent_data, method="quantile")
    
    # 2. Calculate NA proportions safely
    calculate_na_proportion <- function(x) {
        if(ncol(x) == 1) {
            as.numeric(is.na(x))
        } else {
            rowSums(is.na(x)) / ncol(x)
        }
    }
    
    na_proportion_R <- calculate_na_proportion(agilent_data$R)
    na_proportion_G <- calculate_na_proportion(agilent_data$G)
    
    # 3. Create filter mask
    keep <- !is_control & 
            (na_proportion_R <= max_probe_na) & 
            (na_proportion_G <= max_probe_na)
    
    # 4. Apply filtering
    agilent_data <- list(
        R = agilent_data$R[keep, , drop = FALSE],
        G = agilent_data$G[keep, , drop = FALSE],
        genes = agilent_data$genes[keep, ],
        targets = agilent_data$targets
    )

    # Filter low expression
    A_values <- (agilent_data$R + agilent_data$G)/2
    if(!is.matrix(A_values)) A_values <- as.matrix(A_values)
    
    cutoff <- quantile(A_values, intensity_cutoff_quantile, na.rm=TRUE)
    min_samples <- max(1, ceiling(ncol(A_values) * min_samples_above_cutoff))
    keep_expr <- rowSums(A_values > cutoff, na.rm=TRUE) >= min_samples
    
    filtered_data <- agilent_data[keep_expr, ]
    filtered_data$genes <- droplevels(filtered_data$genes)
    
    # 8. Return Results --------------------------------------------------------
    structure(
        list(
            data = filtered_data,
            qc_metrics = list(
                initial_probes = nrow(agilent_data$R) + length(complete_failures$probes),
                final_probes = nrow(filtered_data$R),
                controls_removed = sum(is_control),
                na_probes_removed = sum(!keep),
                low_expr_removed = sum(keep) - sum(keep_expr),
                intensity_cutoff = cutoff,
                min_samples = min_samples,
                missing_values = na_stats
            )
        ),
        class = "agilent_qc_results"
    )
}

# Enhanced print method
print.agilent_qc_results <- function(x, ...) {
    cat("Agilent Two-Color QC Results\n")
    cat("----------------------------\n")
    cat("Final dimensions:", x$qc_metrics$final_probes, "genes ×", 
       ncol(x$data$R), "samples\n\n")
    cat("Filtering metrics:\n")
    cat("- Controls removed:", x$qc_metrics$controls_removed, "\n")
    cat("- High NA probes removed:", x$qc_metrics$na_probes_removed, "\n")
    cat("- Low expression removed:", x$qc_metrics$low_expr_removed, "\n")
    cat("- Final intensity cutoff:", round(x$qc_metrics$intensity_cutoff, 2), "\n")
    cat("- Required samples above cutoff:", x$qc_metrics$min_samples, "\n")
    invisible(x)
}

gse = agilent_two_color_qc(agilent_data1)
print.agilent_qc_results(gse)


# First perform differential expression analysis to find significant probes
# Create design matrix based on final_subtype
pd7
design <- model.matrix(~ 0 + factor(pd7$atl_subtype2))
rownames(design) <- pd7$ID
colnames(design) <- levels(factor(pd7$atl_subtype2))
head(design, 20)
gene_data
# Make contrasts (adjust based on your comparisons)
levels(factor(pd7$atl_subtype2))
contrast.matrix <- makeContrasts(
    ATL_HAMTSP = ATL - HAMTSP,
     ATL_AC = ATL - AC,
     HAM_TSP_AC = HAMTSP - AC,
    levels = design
)
fit2 <- lmFit(gene_data, design) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()

# Get significant probes (FDR < 0.01)
top_probes <- topTable(fit2, number = Inf, adjust.method = "BH",  p.value = 0.01)

top_probes |>
    arrange(adj.P.Val)

topprobes <- top_probes |>
    as.data.frame() |>
    tibble::rownames_to_column("probeid") 
    
topprobes |>
    arrange(adj.P.Val) |>
    head(50)

topprobes |>
    filter(GENE_SYMBOL == "ZCCHC12")

topprobes |>
    mutate(sign = ifelse(adj.P.Val < 0.01 & abs(logFC) > 1.2, "Sign", "No")) |>
    dplyr::mutate(adp = -log10(adj.P.Val)) |>
    tidyplots::tidyplot(x = logFC, y = adp, color = sign) |>
    tidyplots::add_data_points(alpha = .5) |>
    tidyplots::add_data_labels_repel(label = GENE_SYMBOL, data = filter_rows(sign == "Sign"), color = "black") |>
    tidyplots::adjust_x_axis_title("Log2(Fold Change)") |>
    tidyplots::adjust_y_axis_title("-Log10(Adjusted Pvalue)") |>
    tidyplots::remove_legend() |>
    tidyplots::save_plot("ATLvATLc.png",
        bg = "transparent"


top_probes


#### heatmap

degenes = top_probes |>
    tibble::rownames_to_column("probeid") |>
    pull(probeid) 

genes = gse_final$E[degenes,] |>
  as.data.frame() |>
  tibble::rownames_to_column("ENSEMBL_ID") |>
  left_join(annot[,c("ENSEMBL_ID", "GENE_SYMBOL")], by=c("ENSEMBL_ID")) |>
  tibble::column_to_rownames("GENE_SYMBOL") |>
  dplyr::select(-ENSEMBL_ID) |>
  as.matrix()
genes

pd7m = pd7 |>
  tibble::column_to_rownames("ID")
pd7m

# 1. Prepare annotation data - ensure it's properly formatted
annotation_data <- data.frame(
  agec = factor(pd7m$agec),  # Convert to factor explicitly
  atl = factor(pd7m$atl_subtype2),
  row.names = rownames(pd7m)       # Ensure this matches your column names in 'genes'
)

# 2. Verify the annotation data
print(head(annotation_data))
print(table(annotation_data$atl))

# 3. Create color mapping - ensure colors match actual factor levels
categories <- levels(annotation_data$atl)  # Use levels() instead of unique()
categories_2 = levels(annotation_data$agec)
n_categories <- length(categories)
n2_categories = length(categories)

# Check if we have enough colors
if(n_categories != 3) {
  warning(paste("You have", n_categories, "categories but provided 3 colors"))
  # Generate enough colors if needed
  colors <- colorRampPalette(c("#66C2A5", "#FC8D62", "#8DA0CB"))(n_categories)
} else {
  colors <- c("#66C2A5", "#FC8D62", "#8DA0CB")
}

annotation_colors <- list(
  ATLType = setNames(colors, categories)
)

# 4. Verify matrix and annotation dimensions
stopifnot(
  identical(colnames(genes), rownames(annotation_data)),
  ncol(genes) == nrow(annotation_data)
)

# 5. Create the heatmap with error handling

pheatmap(
    mat = genes,
    scale = "row",
    annotation_col = annotation_data,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellheight=15,
    main = "Differential gene expression ",
  filename = "2publication_heatmap.png",
  width = 10,  # Nature standard single-column
  height = 7,
  units = "in",
  family = "Arial"  # Embed font
)
?pheatmap

annotation_data <- data.frame(
  ATLType = pd7$atltype,  # Replace 'atltype' with your actual column name
  row.names = rownames(pd7m)  # Must match your sample names
)
annotation_data 

categories <- unique(annotation_data$ATLType)

annotation_colors <- list(
  ATLType = c("#66C2A5", "#FC8D62", "#8DA0CB")  # Must match factor levels
)
names(annotation_colors$ATLType) <- levels(annotation_data$ATLType)

pheatmap(
  mat = genes,
  scale = "row",
  annotation_col = annotation_data,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,  # Often better for large datasets
  show_colnames = FALSE,
  main = "Gene Expression by ATL Type"
)
```
