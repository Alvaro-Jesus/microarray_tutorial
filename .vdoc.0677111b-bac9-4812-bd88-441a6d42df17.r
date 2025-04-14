#
#
pacman::p_load(
    GEOquery, tidyverse, ggrepel, limma, oligo, DT, pheatmap, tidyplots,  affy, oligoClasses, testit
)
# function
isLog2Transformed <- function(data){
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  shouldBeLogged <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
  return(!shouldBeLogged)
  }


# First obtain metadata
id="GSE33615"
meta = getGEO(id, GSEMatrix=TRUE, destdir=".temp" ) 
meta=meta[[1]]
meta
varLabels(meta)
head(pData(meta))
# GPL4133
pd
pd <- pData(meta) |>
  tibble::rownames_to_column("ID") |>
  dplyr::select("characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3", ID, supplementary_file) |>
  dplyr::mutate(file = str_split(supplementary_file, "/") |> map_chr(tail, 1)) |>
  dplyr::rename(atl_subtype="characteristics_ch1.1", age="characteristics_ch1.2", gender="characteristics_ch1.3") |>
  dplyr::mutate(across(
    c(atl_subtype, age, gender),
    ~ stringr::str_remove(., "^.*: ")  
  )) |>
  mutate(atl_subtype2 = dplyr::recode_factor(atl_subtype, 
    "Acute"="ATLa",
    "Lymphoma"="ATLa",
    "Chronic"="ATLc",
    "Smoldering"="ATLc",
    "HTLV-1 uninfected"="HD",
  )) |>
  dplyr::filter(atl_subtype2!="Unkown")


for (i in 1:length(pd7$supplementary_file)) {
  url <- pd7$supplementary_file[i]
  destfile <- file.path(paste0(".temp/", id, "/", pd7$file[i]))
  
  # Download the file
  tryCatch({
    download.file(url, destfile, mode = "wb", method = "curl")
  }, error = function(e) {
    # Fallback to default method if curl fails
    download.file(url, destfile, mode = "wb")
  })
  
  # Optional: Extract if it's a tar file
  if (grepl("\\.tar(\\.gz)?$", destfile)) {
    untar(destfile, exdir = paste0(".temp/", id))
  }
}

### Re-read the raw files
agilent_data <- read.maimages(
  files = paste0(".temp/", id, "/", pd7$file),
  source = "agilent",
  green.only = TRUE,
  names = pd7$ID,
  columns = list(
    E = "gMedianSignal",       # Foreground signal
    Eb = "gBGMedianSignal"     # Background signal
  ),
  other.columns = "gIsWellAboveBG"  # Quality flag (if column exists)
)
# Convert targets to tibble for easy joining
agilent_data$targets = agilent_data$targets |>
as_tibble(rownames = "sampleName") |>
  left_join(pd7, by = c("sampleName"="ID")) |>
  as.data.frame() |>
  tibble::column_to_rownames("sampleName")

#### gene Annotation
gpl2 = getGEO("GPL4133")
annot <- Table(gpl2)[, c("ID", "GENE_SYMBOL", "ENSEMBL_ID", "SPOT_ID", "CONTROL_TYPE")]
agilent_data$genes <- annot[match(agilent_data$genes$ProbeName, annot$SPOT_ID), ]



# For single-channel EList (it is not obligatory)
eset <- ExpressionSet(assayData = agilent_data$E,
                     phenoData = AnnotatedDataFrame(agilent_data$targets),
                     featureData = AnnotatedDataFrame(agilent_data$genes))


agilent_qc = function(
  agilent_data,
  max_row_na = 0.1 * ncol(agilent_data$E),  # Allow 10% 
  max_col_na = 0.1 * nrow(agilent_data$E)
  ){
    # Calculate NA statistics
    row_na_counts <- rowSums(is.na(agilent_data$E))
    col_na_counts <- colSums(is.na(agilent_data$E))
    # Print summary tables
    message("\nMissing values per ROW (probe/gene):")
    print(table(row_na_counts))
    message("\nMissing values per COLUMN (sample):")
    print(table(col_na_counts))
    # Check for critical issues
    problematic_rows <- sum(row_na_counts > max_row_na)
    problematic_cols <- sum(col_na_counts > max_col_na)
    
    if (problematic_rows > 0 || problematic_cols > 0){
    stop(
      sprintf(
        "\nSANITY CHECK FAILED:\n- %d probes/genes exceed %d missing values\n- %d samples exceed %d missing values\n",
        problematic_rows, max_row_na,
        problematic_cols, max_col_na
        ),
        call. = FALSE)} 
    else {
    message("\nSANITY CHECK PASSED: No excessive missing values detected")}
    # "3 tables"
    # check if # of rows in pdata is same and cols in edata
    # and rows of edata match rows of fdata
    assert(dim(agilent_data$E)[2]==dim(agilent_data$targets)[1])
    assert(dim(agilent_data$E)[1]==dim(agilent_data$genes)[1])
    message("\nSANITY CHECK PASSED: All data match between rows and columns")
    # Quality control
    message("\nQC: background correction")
    agilent_data_qc  <- limma::backgroundCorrect(agilent_data , method="normexp", offset=20) 
    #agilent_data = normalizeWithinArrays(agilent_data ,method="loess") we only have one color
    message("\nQC: normalize between arrays")
    agilent_data_nor <- limma::normalizeBetweenArrays(agilent_data_qc, method="quantile")
    # remove those without a symbol
    message("\nQC: removing those without a symbol")
    agilent_data_rem <- agilent_data_nor[!is.na(agilent_data_nor$genes$ENSEMBL_ID) & 
                                    agilent_data_nor$genes$ENSEMBL_ID != "", ]
    # todo: averaging or take one with best FC?
    print(dim(agilent_data_rem))
    message("\nQC: average for repetitives")
    agilent_data_avg <- limma::avereps(agilent_data_rem, ID=agilent_data_rem$genes$ENSEMBL_ID) 
    print(dim(agilent_data_avg)) # 33492   202
    message("\nQC: final dimensions")
    # is it normalized?
    message("\nIs it normalized:")
    print(isLog2Transformed(agilent_data_avg$E))
    message("\nDensities plot ->:")
    plotDensities(agilent_data_avg, main = "", legend = FALSE)
    message("\nRemoving low expressing genes across...")
    # remove low expressing geneS
    cutoff <- median(agilent_data_avg$E)
    is_expressed <- agilent_data_avg$E > cutoff
    keep <- rowSums(is_expressed) > 2
    ## check how many genes are removed / retained.
    ## subset to just those expressed genes
    gse <- agilent_data_avg[keep,]
    # QC Checks
    message("\nRunning quality checks...")
    print(paste("Original probes:", nrow(agilent_data$E)))
    print(paste("After filtering:", nrow(agilent_data_rem$E)))
    print(paste("After averaging:", nrow(agilent_data_avg$E)))
    print(paste("After averaging:", nrow(gse$E)))
    return(gse)
    }

gse_final = agilent_qc(agilent_data)
dim(agilent_data_avg)

# First perform differential expression analysis to find significant probes
# Create design matrix based on final_subtype
pd7
design <- model.matrix(~0 + factor(pd7$atl_subtype2))
design 
rownames(design) = pd7$ID
colnames(design) <- levels(factor(pd7$atl_subtype2))
head(design, 20)

agilent_data_avg$E |>
    as.data.frame() |>
 tibble::rownames_to_column("id") |>
 group_by(id) |>
 distinct(id)

# Fit linear model
fit <- lmFit(gse, design)
# Make contrasts (adjust based on your comparisons)
levels(factor(pd7$atl_subtype2))
contrast.matrix <- makeContrasts(
  ATLa_ATLc= ATLa - ATLc,
  ATLa_HD= ATLa - HD,
  ATLc_HD= ATLc - HD,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2
# Get significant probes (FDR < 0.01)
top_probes <- topTable(fit2, number = Inf,  adjust.method = "BH", sort.by = "B", p.value = 0.01)
colnames(top_probes)
head(top_probes, 20)
top_probes |>
    tibble::rownames_to_column("probeid") |>
    left_join(annot, by=c("probeid"="ENSEMBL_ID")) |>
    arrange(desc(ATLa_ATLc)) |>
    head(50)

#
#
#
