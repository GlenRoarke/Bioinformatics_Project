Create_TSS_promoters_sites
================
Glen Roarke
2023-06-20

Purpose - This Rmarkdown file provides code to create a TSS for use in
Griffin using the TSS configuration.

# Import data

## griffin sites

This data can be used to run on a small subset of TSS for transcription
factors present in the current griffin sites.

``` r
#sets wd to location of rmd file R studio only not console
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import yaml sites selected from existing sites
data <- yaml.load_file("sites.yaml")

# Extract keys from the data
keys <- names(data$site_lists)

# Create a dataframe with the keys
griffin_sites <- data.frame(Key = keys)

# Mutate the "Key" column by removing everything after the full stop
griffin_sites$Gene_id <- str_remove(griffin_sites$Key, "\\..*")
```

## EPD promoters & chromosome annotation

Some manipulation is required to get the EPD promoter sites into a
format that griffin can use.

``` r
###########

# import promoters, a current site example and human chr annotation
promoter_tss <- vroom("Hs_EPDnew_006_hg38.sga", col_names = FALSE)
```

    ## Rows: 29598 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): X1, X2, X4, X6
    ## dbl (2): X3, X5
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# import chromosome annotation
human_ref <- vroom("GCF_000001405.39_GRCh38.p13_assembly_report.txt")
```

    ## Rows: 640 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (9): # Sequence-Name, Sequence-Role, Assigned-Molecule, Assigned-Molecul...
    ## dbl (1): Sequence-Length
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# import gene drivers from a publication
crc_drivers <- vroom("Driver_genes_41586_2022_5202_MOESM8_ESM.csv", delim = "\t")
```

    ## Rows: 568 Columns: 1
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): IntOGen
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#rename the column used for the join
crc_drivers <- crc_drivers %>% rename(Gene_id = IntOGen)


# generate the UCSC style to accession id mappings
chroms <- human_ref %>% filter(str_detect(`RefSeq-Accn`,"^NC")) %>% 
  select(`RefSeq-Accn`, `UCSC-style-name`) %>% 
  rename(Chroms = `RefSeq-Accn`)


# rename promoter columns to griffin inputs
promoter_tss <- promoter_tss %>% rename(Chroms = X1, TSS = X3, Strand = X4, Gene = X6)


# list of genes and chromosomes
chrom_p <- promoter_tss %>% distinct(Chroms) %>% arrange(Chroms)
genes_p <- promoter_tss %>% distinct(Gene)

# map chromosome description to accession ID

# inner join TSS sites to accession number to retrieve the UCSC style name.
tss <- promoter_tss %>% inner_join(chroms) %>% 
  select(`UCSC-style-name`, TSS, Strand, Gene) %>% rename(Chrom = `UCSC-style-name`) %>% 
  mutate(Gene_id = str_remove(Gene, "_.*"))
```

    ## Joining, by = "Chroms"

An overview of the new sites file structure.

``` r
head(tss)
```

    ## # A tibble: 6 × 5
    ##   Chrom     TSS Strand Gene      Gene_id
    ##   <chr>   <dbl> <chr>  <chr>     <chr>  
    ## 1 chr1   959256 -      NOC2L_1   NOC2L  
    ## 2 chr1   960633 +      KLHL17_1  KLHL17 
    ## 3 chr1   966482 +      PLEKHN1_1 PLEKHN1
    ## 4 chr1   976681 -      PERM1_1   PERM1  
    ## 5 chr1  1000097 -      HES4_1    HES4   
    ## 6 chr1  1000511 +      ISG15_2   ISG15

# create TSS CRC list

There are a large amount of TSS in the EPD file, so here a subset of CRC
drivers is used based on a recent publication.

This create one file for every row in the data frame.

``` r
# match tss sites by gene drivers
tss_crc_drivers <- tss %>% inner_join(crc_drivers)
```

    ## Joining, by = "Gene_id"

``` r
# review any TSS site that are not linked correctly.
tss_missing <- crc_drivers %>% anti_join(tss)
```

    ## Joining, by = "Gene_id"

``` r
tss_missing
```

    ## # A tibble: 29 × 1
    ##    Gene_id 
    ##    <chr>   
    ##  1 ACKR3   
    ##  2 ADGRB1  
    ##  3 AFF3    
    ##  4 CEBPA   
    ##  5 CNOT9   
    ##  6 DCAF12L2
    ##  7 FAM186A 
    ##  8 FANCF   
    ##  9 FAT3    
    ## 10 FOXD4L1 
    ## # … with 19 more rows

``` r
# Specify the folder path
folder_path <- "C:/Users/fh22528/OneDrive - University of Bristol/Project/scripts/griffin_analysis/TFBs_sites_analysis/tss_v2_crc"

# Create the folder if it doesn't exist
dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)

# Loop through each row and write to a separate file in the folder
for (i in 1:nrow(tss_crc_drivers)) {
  # Extract the current row
  row <- tss_crc_drivers[i, ]
  
  # Create a filename based on the row's ID and folder path
  filename <- file.path(folder_path, paste0(row$Gene, ".h38.TSS.txt"))
  
  # Write the row to a file
  write.table(row, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}
```

# All EPD sites

The below code chunk creates a folder of all 30,000 TSS sites available
in the EPD database. The site list can easily be modified by changing
the folder location and data frame object.

``` r
# Specify the folder path
folder_path <- "C:/Users/fh22528/OneDrive - University of Bristol/Project/scripts/griffin_analysis/TFBs_sites_analysis/tss_v1"

# Create the folder if it doesn't exist
dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)

# Loop through each row and write to a separate file in the folder
for (i in 1:nrow(tss)) {
  # Extract the current row
  row <- tss[i, ]
  
  # Create a filename based on the row's ID and folder path
  filename <- file.path(folder_path, paste0(row$Gene, ".h38.TSS.txt"))
  
  # Write the row to a file
  write.table(row, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}
```

``` r
save.image(file = "06_create_TSS_griffin_sites.RData")
```
