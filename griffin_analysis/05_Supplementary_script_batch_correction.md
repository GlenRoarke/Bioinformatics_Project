Supplementary 5 - Batch correction
================
Glen Roarke
2023-08-09

``` r
install.packages("vroom")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("stringr")
install.packages("patchwork")
```

``` r
knitr::opts_chunk$set(echo = TRUE)

library(vroom)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(mixOmics)
library(sva)
library(patchwork)
library(ggbiplot)
```

``` r
load("inputs/HQ_central_coverage_workspace.RData")
```

``` r
b1_preproc_ex <- b1_preproc %>% filter(Mapped_read_pairs >= 0.8)
b2_preproc_ex <- b2_preproc %>% filter(Mapped_read_pairs >= 0.8)
b3_preproc_ex <- b3_preproc %>% filter(Mapped_read_pairs >= 0.8)

# add in batch information as a factor
b1_preproc_ex$batch <- as.factor(1)
b2_preproc_ex$batch <- as.factor(2)
b3_preproc_ex$batch <- as.factor(3)

preproc_all.tm <- bind_rows(b1_preproc_ex, b2_preproc_ex, b3_preproc_ex)
```

``` r
# create clinical data filtered by read mapping rate
clinical_df_all.tm <- preproc_all.tm %>% mutate(patient_id = gsub("^[0]+|[A-Za-z]", "", sample_id))%>%
  inner_join(patient_data, by = "patient_id") %>% 
  dplyr::rename(CAP_TRG = `CAP-TRG`) %>% 
  dplyr::select(patient_id, sample_id, Age, Gender, CAP_TRG, pT, pN,batch) %>%
  dplyr::mutate(Gender = as.character(Gender), CAP_TRG = as.character(CAP_TRG)) %>% 
  dplyr::mutate(CAP_TRG = if_else(is.na(CAP_TRG), "Unknown", CAP_TRG)) %>% 
  dplyr::mutate(Gender = if_else(Gender == 1, "Female", "Male"))


# create a list to filer out of coverage data. 

b1_filter.tm <- b1_preproc_ex$sample_id
b2_filter.tm <- b2_preproc_ex$sample_id
b3_filter.tm <- b3_preproc_ex$sample_id
```

``` r
#create coverage matrices

# create matrices 
#batch 1 
b1_mat_ex.tm <- 
  b1_cov %>% 
  dplyr::filter(sample_id %in% b1_filter.tm) %>% 
  dplyr::select(site_name, sample_id, central_coverage) %>% 
  dplyr::arrange(sample_id) %>% 
  pivot_wider(names_from = sample_id, values_from = central_coverage) %>% 
  column_to_rownames(var = "site_name") %>% 
  as.matrix()

dim(b1_mat_ex.tm)
```

    ## [1] 270  25

``` r
b2_mat_ex.tm <- 
  b2_cov %>% 
  dplyr::filter(sample_id %in% b2_filter.tm) %>% 
  dplyr::select(site_name, sample_id, central_coverage) %>% 
  dplyr::arrange(sample_id) %>% 
  pivot_wider(names_from = sample_id, values_from = central_coverage) %>% 
  column_to_rownames(var = "site_name") %>% 
  as.matrix()

dim(b2_mat_ex.tm)
```

    ## [1] 270  21

``` r
b3_mat_ex.tm <- 
  b3_cov %>% 
  filter(sample_id %in% b3_filter.tm) %>% 
  dplyr::select(site_name, sample_id, central_coverage) %>% 
  arrange(sample_id) %>% 
  pivot_wider(names_from = sample_id, values_from = central_coverage) %>% 
  column_to_rownames(var = "site_name") %>% 
  as.matrix()

dim(b3_mat_ex.tm)
```

    ## [1] 270  25

``` r
# merge all batches into one matrix - read mapping above 0.8
all_cov_mat.tm <- t(cbind(b1_mat_ex.tm, b2_mat_ex.tm, b3_mat_ex.tm))

dim(all_cov_mat.tm)
```

    ## [1]  71 270

``` r
#create a data frame with samples and batches 
coverage.batch <- preproc_all.tm %>% dplyr::select(sample_id, batch)

# only join high quality samples
cap_trg.trt <- preproc_all.tm %>% 
  inner_join(clinical_df_all.tm, by = "sample_id") %>% 
  dplyr::select(sample_id, CAP_TRG) %>%  mutate(CAP_TRG = as.factor(CAP_TRG))
```

``` r
coverage.batch$sample_id
```

    ##  [1] "10A"   "11A"   "12B"   "13B"   "14A"   "16B"   "17A"   "17B"   "18A"  
    ## [10] "18B"   "19A"   "80A"   "80B"   "81A"   "81B"   "82A"   "82B"   "83A"  
    ## [19] "84A"   "85A"   "85B"   "86A"   "87A"   "87B"   "89A"   "10B"   "11B"  
    ## [28] "12A"   "13A"   "15A"   "15B"   "47A"   "47B"   "48A"   "48B"   "49A"  
    ## [37] "49B"   "50Aa"  "63A"   "63B"   "83B"   "84B"   "86B"   "88A"   "88B"  
    ## [46] "89B"   "027A"  "027B"  "027C"  "040D"  "052A"  "052D"  "056A"  "056E" 
    ## [55] "057A"  "057D"  "061B"  "095A"  "101A"  "101D"  "S001A" "S001B" "S001C"
    ## [64] "S002A" "S002B" "S002C" "S003A" "S003B" "S003C" "S004A" "S004B"

``` r
coverage.batch$batch
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2
    ## [39] 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
    ## Levels: 1 2 3

``` r
cap_trg.trt$sample_id
```

    ##  [1] "10A"   "11A"   "12B"   "13B"   "14A"   "16B"   "17A"   "17B"   "18A"  
    ## [10] "18B"   "19A"   "80A"   "80B"   "81A"   "81B"   "82A"   "82B"   "83A"  
    ## [19] "84A"   "85A"   "85B"   "86A"   "87A"   "87B"   "89A"   "10B"   "11B"  
    ## [28] "12A"   "13A"   "15A"   "15B"   "47A"   "47B"   "48A"   "48B"   "49A"  
    ## [37] "49B"   "50Aa"  "63A"   "63B"   "83B"   "84B"   "86B"   "88A"   "88B"  
    ## [46] "89B"   "027A"  "027B"  "027C"  "040D"  "052A"  "052D"  "056A"  "056E" 
    ## [55] "057A"  "057D"  "061B"  "095A"  "101A"  "101D"  "S001A" "S001B" "S001C"
    ## [64] "S002A" "S002B" "S002C" "S003A" "S003B" "S003C" "S004A" "S004B"

``` r
cap_trg.trt$CAP_TRG
```

    ##  [1] Unknown 2       1       0       0       2       3       3       1      
    ## [10] 1       Unknown 0       0       3       3       0       0       3      
    ## [19] 2       Unknown Unknown 2       2       2       2       Unknown 2      
    ## [28] 1       0       3       3       3       3       2       2       0      
    ## [37] 0       3       2       2       3       2       2       3       3      
    ## [46] 2       0       0       0       Unknown Unknown Unknown 0       0      
    ## [55] 1       1       3       Unknown 0       0       Unknown Unknown Unknown
    ## [64] Unknown Unknown Unknown Unknown Unknown Unknown Unknown Unknown
    ## Levels: 0 1 2 3 Unknown

``` r
#apply Combat batch correction parametric
cov.mat.combat = t(ComBat(dat=t(all_cov_mat.tm), batch= coverage.batch$batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE))
```

    ## Found3batches

    ## Adjusting for0covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding parametric adjustments

    ## Adjusting the Data

``` r
# non parametric 
cov.mat.combat.2 = t(ComBat(dat=t(all_cov_mat.tm), batch= coverage.batch$batch, mod=NULL, par.prior=FALSE, mean.only=TRUE))
```

    ## Using the 'mean only' version of ComBat

    ## Found3batches

    ## Adjusting for0covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding nonparametric adjustments

    ## Adjusting the Data

``` r
# mixOmics PCA - may have to revert to precomp here
cov.pca.before <- pca(all_cov_mat.tm, ncomp = 15, scale = TRUE, center = TRUE)
cov.pca.combat <- pca(cov.mat.combat, ncomp = 15, scale = TRUE, center = TRUE)

cov.pca.combat.2 <- pca(cov.mat.combat.2, ncomp = 15, scale = TRUE, center = TRUE)

plot(cov.pca.before, main = "PCA before correction")
```

![](05_Supplementary_script_batch_correction_files/figure-gfm/Batchcorrection%20and%20PCA-1.png)<!-- -->

``` r
plot(cov.pca.combat, main = "PCA after ComBat correction")
```

![](05_Supplementary_script_batch_correction_files/figure-gfm/Batchcorrection%20and%20PCA-2.png)<!-- -->

``` r
# plot the samples

#before
pca.before <- plotIndiv(cov.pca.before, comp = c(1, 2), 
          group = coverage.batch$batch,
          #pch = cap_trg.trt$CAP_TRG,
          ellipse = TRUE,
          ellipse.level = 0.95,
          legend = T,
          legend.title = "Batch",
          legend.title.pch = 'CAP_TRG',
          title = "Before correction")
```

![](05_Supplementary_script_batch_correction_files/figure-gfm/Batchcorrection%20and%20PCA-3.png)<!-- -->

``` r
ggsave("final_figures/Figure_6.1_PCA_before_correction.pdf", height = 4, width = 6)

# after parametric
plotIndiv(cov.pca.combat, comp = c(1, 2), 
          group = coverage.batch$batch,
          #pch = cap_trg.trt$CAP_TRG,
          ellipse = TRUE,
          ellipse.level = 0.95,
          legend = T,
          legend.title = "Batch",
          legend.title.pch = 'CAP_TRG',
          title = "after Combat correction- parametric")
```

![](05_Supplementary_script_batch_correction_files/figure-gfm/Batchcorrection%20and%20PCA-4.png)<!-- -->

``` r
ggsave("final_figures/Figure_6.1_PCA_after_correction_para.pdf", height = 4, width = 6)


#after non-parametric
plotIndiv(cov.pca.combat.2, comp = c(1, 2), 
          group = coverage.batch$batch,
          #pch = cap_trg.trt$CAP_TRG,
          ellipse = TRUE,
          ellipse.level = 0.95,
          legend = T,
          legend.title = "Batch",
          legend.title.pch = 'CAP_TRG',
          title = "Combat correction - non-parametric")
```

![](05_Supplementary_script_batch_correction_files/figure-gfm/Batchcorrection%20and%20PCA-5.png)<!-- -->

``` r
ggsave("final_figures/Figure_6.1_PCA_after_correction_non-para.pdf", height = 4, width = 6)


# Open pdf file
pdf(file= "final_figures/Figure_6.2_Scree_plot_batch_effect.pdf", height = 4, width = 8)

# create a 1X2 grid
par( mfrow= c(1,3) )

plot(cov.pca.before, main = "PCA before correction")
plot(cov.pca.combat, main = "PCA after ComBat correction")
plot(cov.pca.combat.2, main = "PCA ComBat non-parametric")

# Close the PDF device and save the plots to the file
dev.off()
```

    ## png 
    ##   2

# Save data

``` r
save.image(file = "inputs/05_PCA_Correction.RData")
```

# SessionInfo

``` r
#sessionInfo()
installed.packages()[names(sessionInfo()$otherPkgs), "Version"]
```

    ##     ggbiplot       scales         plyr          sva BiocParallel   genefilter 
    ##       "0.55"      "1.2.1"      "1.8.8"     "3.46.0"     "1.32.6"     "1.80.3" 
    ##         mgcv         nlme     mixOmics      lattice         MASS    patchwork 
    ##     "1.8-40"    "3.1-157"     "6.22.0"    "0.20-45"     "7.3-57"      "1.1.2" 
    ##      forcats      stringr        dplyr        purrr        readr        tidyr 
    ##      "0.5.2"      "1.4.1"     "1.0.10"      "0.3.5"      "2.1.3"      "1.2.1" 
    ##       tibble    tidyverse      ggplot2        vroom 
    ##      "3.1.8"      "1.3.2"      "3.4.0" "1.6.0.9000"
