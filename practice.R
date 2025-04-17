#--------------------------------------------------------
# Title: Microarray Data Analysis with GEOquery and limma

#--------------------------------------------------------
# Part 1: loading packages and functions
#--------------------------------------------------------
pacman::p_load(
    GEOquery, tidyverse, ggrepel, limma, oligo, DT, pheatmap, tidyplots, affy, oligoClasses, testit
)

# function is log2transformed
isLog2Transformed <- function(data) {
    qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    shouldBeLogged <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
    return(!shouldBeLogged)
}

#--------------------------------------------------------
# Part 2: Obtaining metadata and raw data
#--------------------------------------------------------

# First obtain metadata
id <- "GSE19080" # GPL4133

# creating a temp folder
if (!dir.exists(paste0(".temp/", id))) {
    dir.create(paste0(".temp/", id))
} else {
    message("Folder already exists!")
}

# Obtaining Metadata
meta <- getGEO(id, GSEMatrix = TRUE, destdir = ".temp")
meta <- meta[[1]]
head(fData(meta))
head(pData(meta))
colnames(exprs(meta))
dim(meta)

# Metadata wrangling
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

## Download all the files in the .temp of this environment
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

#--------------------------------------------------------
# Part 3: Capturing targets
#--------------------------------------------------------
# Exploring the structure of one file
con <- gzfile(file.path("/Users/denriquez/Documents/GitHub/microarray_tutorial/.temp/GSE19080/GSM472372_HISH0553.txt.gz"))
file_lines <- readLines(con, n=1000)
close(con)

# Subsettting files (here we have two different platforms) 
pd7 = pd |>
    filter(ID%in%c(paste0("GSM4723", c(56:73, 82:93))))

# Re-read the raw files, this was a quantarray (described in metadata)
agilent_data <- read.maimages(
    files = file.path(".temp", id, pd7$file),
    source = "quantarray",
    green.only = FALSE,
    names = pd7$ID,
    other.columns = list(
        Flag = "Ignore Filter"))


#--------------------------------------------------------
# Part 4: Manipulating targets
#--------------------------------------------------------

# Convert targets to tibble for easy joining

agilent_data$targets = agilent_data1$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(pd7, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")

#### gene Annotation
gpl2 <- getGEO("GPL9686")

dim(fData(meta)[,c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")])

annot <- Table(gpl2)[, c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")] 

# Add gene symbols to agilent_data1
agilent_data$genes$GENE_SYMBOL <- annot$SYMBOL[match(agilent_data$genes$Name, annot$GB_ACC)]


agilent_data$genes=agilent_data$genes |>
    mutate(is_control=factor(ifelse(is.na(GENE_SYMBOL), "control", "gene")))

table((agilent_data$genes$is_control))
head(agilent_data$genes)

#--------------------------------------------------------
# Part 5: QC
#--------------------------------------------------------
# (A) Use negative controls for background correction

agilent_two_color_qc <- function(agilent_data,
                                verbose = TRUE) {
    # 1. Input Validation --------------------------------------------------------
    required_components <- c("R", "G", "genes", "targets")
    missing_comps <- setdiff(required_components, names(agilent_data))
    if(length(missing_comps)) {
        stop("Missing required components: ", paste(missing_comps, collapse=", "))
    }
    # 2. Dimensional Checks ----------------------------------------------------
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

    # 3. Missing Value Analysis ------------------------------------------------
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
    
    # 4. Background correction
    agilent_data <- limma::backgroundCorrect(agilent_data, method="normexp", offset=20) 
    # 5. Normalization
    agilent_data <- limma::normalizeWithinArrays(agilent_data, method="loess")
    agilent_data <- limma::normalizeBetweenArrays(agilent_data, method="quantile")

    # 6. technical replicates reduction by probename
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$Name)
    agilent_data <- agilent_data[agilent_data$genes$is_control=="gene", ]
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$GENE_SYMBOL)
    agilent_data$genes <- droplevels(agilent_data$genes)
    return(agilent_data)
}
colnames(agilent_data$genes)

gse1 = agilent_two_color_qc(agilent_data)

##### PCA

dim(gse1)
boxplot(gse1$M, main="Normalized Negative Controls")
pd7
pca <- prcomp(t(gse1$M))
rownames(pca$x)
pca = pca$x |>
    as.data.frame() 
pca |>
    tibble::rownames_to_column("ID") |>
    dplyr::left_join(pd7, by="ID") |> 
    tidyplots::tidyplot(x=PC1, y=PC2, color=atl_subtype) |>
    tidyplots::add_data_points()

#--------------------------------------------------------
# Part 6: DE
#--------------------------------------------------------
# Create design matrix based on final_subtype
design <- model.matrix(~ 0 + factor(pd7$atl_subtype2))
rownames(design) <- pd7$ID
colnames(design) <- levels(factor(pd7$atl_subtype2))
head(design, 20)

# Make contrasts (adjust based on your comparisons)
levels(factor(pd7$atl_subtype2))
contrast.matrix <- makeContrasts(
    ATL_HAMTSP = ATL - HAMTSP,
    ATL_AC = ATL - AC,
    HAM_TSP_AC = HAMTSP - AC,
    levels = design
)
fit2 <- lmFit(gse1, design) %>%
    contrasts.fit(contrast.matrix) %>%
    eBayes()

# Get significant probes (FDR < 0.01)
top_probes <- topTable(fit2, number = Inf, adjust.method = "BH") #,  p.value = 0.01)
top_probes

#--------------------------------------------------------
# Part 7: Manipulation
#--------------------------------------------------------
colnames(top_probes)

"""
Proximal TCR signaling: CTLA-4, PD-1, SHP-1, SHP-2, LYP, Cbl-b, GRAIL, SIT, PAG, Dok family members, Drak2, and CD5. 
These regulators often target kinases like Lck and ZAP-70, adaptor proteins like LAT and SLP-76, or the CD3ζ chain itself.   
Calcium signaling: Some negative regulators, such as CTLA-4 and Dok family members, can influence calcium mobilization and downstream signaling events.   
MAPK pathway: Negative regulation of the MAPK pathway, which is critical for gene expression and effector functions, is mediated by CTLA-4, PD-1, SHP-2, Dok family members, and RASA2.   
NF-κB pathway: Several negative regulators, including CTLA-4, PD-1, Cbl-b, Itch, A20, and Peli1, target the NF-κB signaling pathway, which controls the expression of genes involved in inflammation, apoptosis, and immune responses.   
NFAT pathway: CTLA-4, PD-1, and MDM2 are among the negative regulators that can modulate the activity of NFAT transcription factors, which are critical for cytokine production.   
PI3K/Akt/mTOR pathway: The PI3K/Akt/mTOR pathway, involved in cell survival, growth, and metabolism, is targeted by negative regulators such as PD-1, Cbl-b, Peli1, and TSC1/TSC2.
"""
probes = top_probes |>
    dplyr::select(-c(1:5, 7)) |>
    as_tibble() |>
    arrange((ATL_AC)) #|>

probes |>
    dplyr::filter(GENE_SYMBOL=="NT5E") # |>
    dplyr::filter(ATL_AC > ATL_HAMTSP, HAM_TSP_AC>0, adj.P.Val < 0.05)  |>
    pull(GENE_SYMBOL) 

early_tcr=c("LCK", "CD3D", "CD3E", "CD3G", "LAT", "LCP2") # ZAP70, SLP76, ITK

neg_tcr=c("CD5", "PTPN6", "PTPN22","SOCS1", "CBLB", "DUSP14", "CD6", "PDCD1LG2", "PTPN12", "DOK1", "MDM2", "JUNB", "SHC1", "UBD", "ATM", "BTK", "IL1RL1", "IFNGR2", "IFNGR1", "GNAI3", "SOCS2") #DRAK2 , PTPN11


probess

probess = probes |>
    dplyr::left_join(tcr_regulation_df, by=c("GENE_SYMBOL"="elements")) |>
    pull("GENE_SYMBOL")

tcr_regulation_df <- data.frame(
  elements = c(early_tcr, neg_tcr),
  type = c(rep("early_tcr", length(early_tcr)), 
           rep("neg_tcr", length(neg_tcr)))
)
neg_tcr

annotationdf =data.frame(
    row.names = pd7$ID,
    Subtype=pd7$atl_subtype2)

rownames(annotationdf)
heatmapgenes = gse1$A[neg_tcr,]
heatmapgenes
gse1$A[c("NFATC4"),]

probess 
pheatmap(
    mat = heatmapgenes,
    scale = "row",
    annotation_col = annotationdf,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellheight=15,
    main = "Differential gene expression ") #,
  filename = "2publication_heatmap.png",
  width = 10,  # Nature standard single-column
  height = 7,
  units = "in",
  family = "Arial"  # Embed font
)
library(DT)
annot |>
    dplyr::filter(str_detect(GENE_NAME, "adenosin") | str_detect(SYMBOL, "NT5E")) #|>
    pull(SYMBOL)

library(org.Hs.eg.db)
library(AnnotationDbi)

IF
colnames(annot)
topprobes <- top_probes |>
    as.data.frame() |>
    tibble::rownames_to_column("probeid") 
    
topprobes |>
    arrange(adj.P.Val) |>
    head(50)

top_probes |>
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
    scale = "column",
    #annotation_col = annotation_data,
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

#####################
pd7_2 = pd |>
    filter(ID%in%c(paste0("GSM4723", 74:81)))
pd7_2
agilent_data2 <- read.maimages(
    files = file.path(".temp", id, pd7_2$file),
    source = "genepix",
    green.only = FALSE,
    names = pd7_2$ID)

################



#--------------------------------------------------------
# Part 3B: HD - Capturing targets
#--------------------------------------------------------
# Subsettting files (here we have two different platforms) 
hd_pd7 = pd |>
    filter(ID%in%c(paste0("GSM4723", 74:81)))

# Re-read the raw files, this was a quantarray (described in metadata)
agilent_data <- read.maimages(
    files = file.path(".temp", id, hd_pd7$file),
    source = "genepix",
    green.only = FALSE,
    names = hd_pd7$ID)

head(agilent_data)
#--------------------------------------------------------
# Part 4: Manipulating targets
#--------------------------------------------------------

# Convert targets to tibble for easy joining

agilent_data$targets = agilent_data$targets |>
    as_tibble(rownames = "sampleName") |>
    left_join(hd_pd7, by = c("sampleName" = "ID")) |>
    as.data.frame() |>
    tibble::column_to_rownames("sampleName")

agilent_data$targets
#### gene Annotation

gpl2 <- getGEO("GPL9686")

dim(fData(meta)[,c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")])

annot <- Table(gpl2)[, c("ID", "SYMBOL", "GENE_NAME", "GB_ACC")] 

head(agilent_data$genes)
# Add gene symbols to agilent_data1
agilent_data$genes <- agilent_data$genes |>
    tidyr::separate(
    Name, 
    into = c("GENE_SYMBOL", "Annot", "probename"), 
    sep = ":"  # Adjust delimiter
  ) |>
  dplyr::left_join(annot[, c("SYMBOL", "GB_ACC")], by=c("probename"="GB_ACC")) 

table(agilent_data$genes$is_control)


agilent_data$genes=agilent_data$genes |>
    mutate(is_control=factor(ifelse(is.na(SYMBOL), "control", "gene")))

table((agilent_data$genes$is_control))
table((agilent_data$genes$is_control))

#--------------------------------------------------------
# Part 5: QC
#--------------------------------------------------------
# (A) Use negative controls for background correction

agilent_two_color_qc <- function(agilent_data,
                                verbose = TRUE) {
    # 1. Input Validation --------------------------------------------------------
    required_components <- c("R", "G", "genes", "targets")
    missing_comps <- setdiff(required_components, names(agilent_data))
    if(length(missing_comps)) {
        stop("Missing required components: ", paste(missing_comps, collapse=", "))
    }
    # 2. Dimensional Checks ----------------------------------------------------
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

    # 3. Missing Value Analysis ------------------------------------------------
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
    
    # 4. Background correction
    agilent_data <- limma::backgroundCorrect(agilent_data, method="normexp", offset=20) 
    # 5. Normalization
    agilent_data <- limma::normalizeWithinArrays(agilent_data, method="loess")
    agilent_data <- limma::normalizeBetweenArrays(agilent_data, method="quantile")

    # 6. technical replicates reduction by probename
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$Name)
    agilent_data <- agilent_data[agilent_data$genes$is_control=="gene", ]
    agilent_data <- limma::avereps(agilent_data, ID=agilent_data$genes$GENE_SYMBOL)
    agilent_data$genes <- droplevels(agilent_data$genes)
    return(agilent_data)
}


gse = agilent_two_color_qc(agilent_data)
dim(gse1)
dim(gse)
class(gse)
boxplot(gse$M, main="Normalized Negative Controls")
pca <- prcomp(t(gse$M))
rownames(pca$x)
pca = pca$x |>
    as.data.frame() 
pca |>
    tibble::rownames_to_column("ID") |>
    dplyr::left_join(pd7, by="ID") |> 
    tidyplots::tidyplot(x=PC1, y=PC2, color=atl_subtype) |>
    tidyplots::add_data_points()