library(dplyr)
library(tidyr)
library(adegenet)
library(radiator)
library(DescTools)
library(factoextra)


#' convert a long-format dataframe into a genind for pca plotting
#' 
#' @param long_df_strata the dataframe with genotypes formatted with populations or species as strata
#' @param outfile name of file for outputting PDF
#' @details Genotypes should be coded with 3 integers for each alleles. 
#' 6 integers in total for the genotypes. e.g. 001002 or 111333 
#' (for heterozygote individual). 6 integers WITH separator: e.g. 001/002 or 111/333 
#' (for heterozygote individual). The separator can be any of these: 
#' "/", ":", "_", "-", ".", and will be removed.

convert_to_genind <- function(long_df_strata) {
  
  # reformat genotypes for genind conversion
  # create 3 digit integers from the genotypes
  long_df_strata$GT3 <- Format(long_df_strata$GT, ldigits = 3, digits = 0)
  
  # fix NAs
  long_df0s <- long_df_strata %>%
    mutate(GT3 = ifelse(is.na(GT3), "000", GT3))
  
  # Now combine the GT3 column per indiv/marker:
  # make the genos characters and then try pasting them as strings
  long_df0s$GT3 <- as.character(long_df0s$GT3)
  
  long_df3digit <- long_df0s %>%
    group_by(INDIVIDUALS, MARKERS) %>% 
    arrange(GT3, .by_group = TRUE) %>% 
    summarise(GENOTYPE = toString(GT3))
  
  # paste strings together
  long_df3digit$GENOTYPE <- gsub(", ","",long_df3digit$GENOTYPE)
  
  # add back on species identity as strata
  df_for_conversion <- long_df0s %>% 
    select(-GT, -GT3) %>%
    left_join(., long_df3digit) %>%
    unique() %>%
    rename(GT = GENOTYPE) %>%
    mutate(GT = ifelse(GT == "000000", NA, GT))
  
  df_for_conversion$STRATA <- as.factor(df_for_conversion$STRATA)
  
  # use the radiator package for this conversion
  df_for_conversion
  
  print(df_for_conversion)
  
}


# convert to genind object and then plot PCA into PDF

convert_to_genind_pca <- function(long_df_strata, outfile) {
  
      # reformat genotypes for genind conversion
      # create 3 digit integers from the genotypes
      long_df_strata$GT3 <- Format(long_df_strata$GT, ldigits = 3, digits = 0)
      
      # fix NAs
      long_df0s <- long_df_strata %>%
        mutate(GT3 = ifelse(is.na(GT3), "000", GT3))
      
      # Now combine the GT3 column per indiv/marker:
      # make the genos characters and then try pasting them as strings
      long_df0s$GT3 <- as.character(long_df0s$GT3)
      
      long_df3digit <- long_df0s %>%
        group_by(INDIVIDUALS, MARKERS) %>% 
        arrange(GT3, .by_group = TRUE) %>% 
        summarise(GENOTYPE = toString(GT3))
      
      # paste strings together
      long_df3digit$GENOTYPE <- gsub(", ","",long_df3digit$GENOTYPE)
      
      # add back on species identity as strata
      df_for_conversion <- long_df0s %>% 
        select(-GT, -GT3) %>%
        left_join(., long_df3digit) %>%
        unique() %>%
        rename(GT = GENOTYPE) %>%
        mutate(GT = ifelse(GT == "000000", NA, GT))
      
      df_for_conversion$STRATA <- as.factor(df_for_conversion$STRATA)
      
      # use the radiator package for this conversion
      genind_df <- write_genind(df_for_conversion)
      
      
      # Now that the data is a genind object, go ahead and run the PCA.
      # Make PCA
      
      # Allele presence absence data are extracted and NAs replaced using tab:
      datasetX <- tab(genind_df, NA.method="mean") # double check that is this the appropriate method.
      
      # make PCA
      dataset_pca1 <- dudi.pca(datasetX, center = TRUE, scannf = FALSE, scale=FALSE, nf = 10)
      
      # colors
      mycol <- c("darkcyan", "deepskyblue", "salmon", "darkblue", "sienna", "slateblue", "darkorange","firebrick", "purple4", "darkmagenta", "goldenrod")
      
      # visualize the scree plot
      fviz_eig(dataset_pca1)
      
      # plot
      # fviz_pca_ind(dataset_pca1,
      #              col.ind = fac(pop(genind_df)),
      #              palette = mycol,
      #              addEllipses = FALSE,
      #              repel = TRUE)
      
      # # plot with factor labels
      pdf(paste0("pdf_outputs/",outfile,".pdf"), width = 10, height = 10)
      s.class(dataset_pca1$li, fac=pop(genind_df), wt = rep(1, length(pop(genind_df))), clabel = .8, grid = FALSE, axesell=FALSE,
              xax=1, yax=2, col=transp(mycol,.8),
               cstar=0, cpoint=0.9)
      dev.off()


}