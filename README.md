# DDX41 Curation Package

This is a R package containing methods to support classification of pathogenicity of DDX41 variants according to ACMG/AMP criteria.


## Installation
```r
# install devtools if you don't have it
# install.packages("devtools")
library(devtools)
install_github("blomberylab/ddx41")
```

## Using the package
To use the package, users need to provide gnomAD, ToMMo and KOVA population datasets (these datasets need to have some specific column names which can be seen in the example datasets provided), and a dataset with known variant curations (to be used for PS1 and PM5 curations).
DDX41 variant counts from existing literature, In silico values such as REVEL score, Grantham score, spliceAI values and AlphaMissense class are also required.
In addition, variant information such as hgvsc, hgvsp, cds start position , hg38 start position, exon and consequence are also needed.

### Included Datasets
```r
library(ddx41)

#Load datasets, you can use your own with the same column names and data types
bulk_curations_file <- system.file("extdata", "BulkCurations.csv", package = "ddx41")
gnomad4_file <- system.file("extdata", "gnomAD_v4.1.0_ENSG00000183258_2024_09_13_20_49_51.csv", package = "ddx41")
tommo_kova_file <- system.file("extdata", "tommo_KOVA.ddx41.csv", package = "ddx41")
bulk_curations <- read.csv(bulk_curations_file)
gnomad4 <- read.csv(gnomad4_file)
tommo_kova <- read.csv(tommo_kova_file)
```

### Example
Given below is an example, with information needed for variant curations.
``` r
#variant info, user to provide
variant_info <- list()
variant_info$hgvsc <- "c.415_418dup"
variant_info$hgvsp <- "p.Asp140GlyfsTer2"
variant_info$cds_start <- 414
variant_info$exon <- "5/17"  
variant_info$consequence_VEP <- "frameshift_variant"
variant_info$revel_score <- NA
variant_info$alphamissense_class <- NA
variant_info$grantham_score <- NA
variant_info$spliceAI_score <- c("DS_AL"=0.03, "DS_DL"=0.02, "DS_AG"=0.04, "DS_DG"=0)
variant_info$literature_gcount <- 160
variant_info$literature_cases <- 34164
variant_info$literature_somatic <- c("p.Arg525His", "p.Gly530Asp", "p.Pro321Leu", "p.Cys264Tyr", "p.Thr232Ala", "p.Asp344Glu",
                                     "p.Glu345Asp", "p.Lys494Glu", "p.Ala376Val", "p.Gly402Leu", "p.Gly530Ala", "p.Gly530Arg",
                                     "p.Gly589Arg", "p.Leu126Ser", "p.Leu269Met", "p.Leu548His", "p.Leu598Ile", "p.Leu87His",
                                     "p.Met378Ile", "p.Pro265Leu", "p.Pro379Gln", "p.Val408Phe")  #Somatic HGVsps
variant_info$total_pp4 <- 235
variant_info$single_recurrent <- 108
variant_info$single_non_recurrent <- 19
variant_info$multi_somatic <- 4

#Modify config_values if you want to change the default settings and configurations
print(config_values)

result <- ddx41::curateVariant(variant_info, gnomad4, tommo_kova, bulk_curations, config_values)
```

### Curation Results
Output of curateVariant function is a named list.
```r
#Final curation output
result$decision
#Individual rule decisions
result$rule_curations
```

### Example input dataset
We have also provided a dataset containing values required as input for some variants.
```r
variant_info_file <- system.file("extdata", "VariantInfo.csv", package = "ddx41")
data <- read.csv(variant_info_file, stringsAsFactors = F)
```

## Online Curation Tool
Check our online DDX41 curation tool at https://blombery-lab.shinyapps.io/ddx41/ if you are interested, which offers single and bulk modes for DDX41 variant curations.

## Contact Us
If you have any questions, require assistance, or would like to provide feedback, please do not hesitate to get in touch with us.
You may contact us via email at: molecular.haematology@petermac.org
