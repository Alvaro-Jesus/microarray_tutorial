#
#
#
#
pacman::p_load(
    GEOquery, tidyverse, ggrepel, limma, oligo, DT, marray, pheatmap, tidyplots,  affy, oligoClasses, Biobase
)
pacman::p_load(
    hgu133plus2.db,
    org.Hs.eg.db, reshape2, survminer, arrayQualityMetrics, testit, biomaRt, AgiMicroRna,
)
#
#
#
#
# First obtain metadata
GSE55851_meta = getGEO("GSE55851", GSEMatrix=TRUE, destdir=".temp" ) 
GSE55851_meta=GSE55851_meta[[1]]

# Explore columns
class(GSE55851_meta)
varLabels(GSE55851_meta)
head(pData(GSE55851_meta))

# Data wrangling
pd <- pData(GSE55851_meta) |>
  tibble::rownames_to_column("ID") |>
  dplyr::select("gender:ch1", "atl subtype:ch1", "age:ch1", ID, supplementary_file) |>
  dplyr::mutate(file = str_split(supplementary_file, "/") |> map_chr(tail, 1)) |>
  dplyr::rename(gender="gender:ch1", atl_sub="atl subtype:ch1", age= "age:ch1") |>
  mutate(atltype = dplyr::recode_factor(atl_sub, 
    "ATL-acute"="ATLa",
    "ATL-chronic"="ATLc",
    "ATL-smoldering"="ATLc",
    "HTLV-1 asymptomatic carrier"="AC",
    "Normal"="HD"
  ), cd7=case_when(
      str_detect(file, "-D.") ~ "C7D", 
      str_detect(file, "-P") ~"CD7P",
      str_detect(file, "N") ~"CD7N",
      .default="U"))

### for CD7+CADM1+ (infected CD4 is homogenous) across all subtypes
### for CD7+CADM1- non-infected CD4 is homogeneous across all subtypes
### for CD7-CADM1+ (neoplastic cells) are homogeneous across all subtypes
## Exlude HD to evaluate DE genes in CD7+CADM1- vs CD7+CADM1+

pd7 = pd |>
    dplyr::filter(atltype!="HD", cd7!="C7D")

pd7
### Re-read the raw files
agilent_data <- read.maimages(
  files = paste0(".temp/GSE55851/", pd7$file),
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
  left_join(pd, by = c("sampleName"="ID")) |>
  as.data.frame() |>
  tibble::column_to_rownames("sampleName")


#### gene Annotation
gpl2 = getGEO("GPL10332")
annot <- Table(gpl2)[, c("ID", "GENE_SYMBOL", "ENSEMBL_ID", "SPOT_ID", "CONTROL_TYPE")]
agilent_data$genes <- annot[match(agilent_data$genes$ProbeName, annot$SPOT_ID), ]


#
#
#
#
# For single-channel EList (it is not obligatory)
eset <- ExpressionSet(assayData = agilent_data$E,
                     phenoData = AnnotatedDataFrame(agilent_data$targets),
                     featureData = AnnotatedDataFrame(agilent_data$genes))


# Sanity check
agilent_data
table(rowSums(is.na(agilent_data$E)))
table(colSums(is.na(agilent_data$E)))

# "3 tables"
# check if # of rows in pdata is same and cols in edata
# and rows of edata match rows of fdata
assert(dim(agilent_data$E)[2]==dim(agilent_data$targets)[1])
assert(dim(agilent_data$E)[1]==dim(agilent_data$genes)[1])

# Quality control
agilent_data  <- limma::backgroundCorrect(agilent_data , method="normexp", offset=20) 
#agilent_data = normalizeWithinArrays(agilent_data ,method="loess") we only have one color
agilent_data <- limma::normalizeBetweenArrays(agilent_data, method="quantile")
agilent_data <- agilent_data[!(is.na(agilent_data$genes$ENSEMBL_ID)), ]
#
#
#
