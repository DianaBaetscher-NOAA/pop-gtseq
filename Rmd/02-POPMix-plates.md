02-POPMix-plates
================
Diana Baetscher
2023-03-10

Data from 8 March 2023. POP mixed plate genotyped with POP gtseq panel
combined (merged) with POP gtseq3 data originally processed in
`01-initial-analysis.Rmd`.

Using the rds file output from microhaplot: 1. read in rds files 2.
apply read depth filters 3. apply allele balance filter

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.4.1      v purrr   0.3.4 
    ## v tibble  3.1.2      v dplyr   1.0.10
    ## v tidyr   1.2.0      v stringr 1.4.0 
    ## v readr   1.4.0      v forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(readxl)
library(stringr)
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

``` r
source("../R/microhaplot-genos-funcs.R")
source("../R/pca-funcs.R")
```

    ## Loading required package: ade4

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## 
    ##    /// adegenet 2.1.3 is loaded ////////////
    ## 
    ##    > overview: '?adegenet'
    ##    > tutorials/doc/questions: 'adegenetWeb()' 
    ##    > bug reports/feature requests: adegenetIssues()

    ## Warning: package 'DescTools' was built under R version 4.1.3

    ## Registered S3 method overwritten by 'DescTools':
    ##   method         from 
    ##   reorder.factor gdata

``` r
#### Call genos from the microhaplot rds files ####

# this next step uses a little data frame that has the names and GTSeq run numbers
# I got that file like this:
# dsb:rds_files dianabaetscher$ pwd
# /Users/dianabaetscher/Desktop/NOAA_grad/git-repos/rockfish-species-id/new_baseline_data/feather_files
# dsb:feather_files dianabaetscher$ ls -l | awk 'BEGIN {print "gtseq_run", "file"} NR >1 && !/posinfo/ {num = $NF; gsub(/gtseq/, "", num); gsub(/[._].*$/, "", num);  print num, $NF}' > ../rds-file-list.txt 


# get the names of the files
fdf <- read.table("../data/rds-file-list.txt", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df()
```

    ## Warning: `tbl_df()` was deprecated in dplyr 1.0.0.

    ## Warning: Please use `tibble::as_tibble()` instead.

``` r
dir <- "../data/rds_files/"


# cycle over them, read them and add the gtseq_run column on each.
# at the end, bind them together.
genos_long <- lapply(1:nrow(fdf), function(i) {
  message("Working on ", fdf$file[i])
  call_genos_from_haplotRDS(path = file.path(dir, fdf$file[i])) %>%
    mutate(gtseq_run = fdf$gtseq_run[i]) %>%
    select(gtseq_run, everything())
}) %>%
  bind_rows()
```

    ## Working on POPMix_gtseq.rds

    ## Joining, by = c("id", "locus", "rank")Working on POP_gtseq3.rds
    ## Joining, by = c("id", "locus", "rank")

``` r
# we go ahead and save it in data/processed, with xz compression
saveRDS(genos_long, file = "../data/processed/called_genos.rds", compress = "xz")


#### In the end, let us get a data frame that includes genotypes for all the individuals  ####
# and which explicitly has NAs in places where data are missing, and also 
# has the NMFS_DNA_ID on there
genos_long_explicit_NAs <- genos_long %>%
  select(gtseq_run, id) %>%
  unique() %>%
  unite(col = gid, sep = "_", gtseq_run, id) %>%
  select(gid) %>%
  unlist() %>%
  unname() %>%
  expand.grid(gid = ., locus = unique(genos_long$locus), gene_copy = 1:2, stringsAsFactors = FALSE) %>%
  tbl_df() %>% 
  separate(gid, into = c("gtseq_run", "id"), convert = TRUE) %>%
  left_join(., genos_long) %>%
  arrange(gtseq_run, id, locus, gene_copy)
```

    ## Joining, by = c("gtseq_run", "id", "locus", "gene_copy")

``` r
# and then save that
saveRDS(genos_long_explicit_NAs, file = "../data/processed/called_genos_na_explicit.rds", compress = "xz")
```

``` r
genos_long_explicit_NAs %>%
  group_by(locus) %>%
  filter(!is.na(allele.balance)) %>%
  summarise(loc_depth = sum(depth)) %>%
  arrange(desc(loc_depth))
```

    ## # A tibble: 140 x 2
    ##    locus               loc_depth
    ##    <chr>                   <int>
    ##  1 CM025826.1_6377992     338911
    ##  2 CM025842.1_25291879    250965
    ##  3 CM025830.1_26009542    210885
    ##  4 CM025825.1_43208600    207708
    ##  5 CM025833.1_2596880     196741
    ##  6 CM025834.1_3818587     177467
    ##  7 CM025844.1_11631411    155190
    ##  8 CM025833.1_28093856    144327
    ##  9 CM025825.1_76785270    143547
    ## 10 CM025845.1_12625972    141965
    ## # ... with 130 more rows

## Locus evaluation

How many loci and how many alleles?

``` r
# alleles
genos_long %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique()
```

    ## # A tibble: 437 x 2
    ##    locus               allele                                                   
    ##    <chr>               <chr>                                                    
    ##  1 CM025825.1_24578215 ATTG                                                     
    ##  2 CM025825.1_26537970 CCCA                                                     
    ##  3 CM025825.1_43208600 GAGGCTGTCCATGTTCAATCTCGTAGAG                             
    ##  4 CM025825.1_43208600 GAGGCTGTCCATGTTCAATCACGTAGAG                             
    ##  5 CM025825.1_61422707 CCTCACTGCACGGGGGCGCCCTCTGGTTGCAAAAAGAAGTCTGATTGATAGAAGTC~
    ##  6 CM025825.1_61422707 CCTCACTGTACGGGGGCGCCCTCTGGTTGCAAAAAGAAGTCTGATTGATAGAAGTC~
    ##  7 CM025825.1_76785270 ACAGG                                                    
    ##  8 CM025825.1_7920395  TCCTCCT                                                  
    ##  9 CM025825.1_7920395  TCCTACT                                                  
    ## 10 CM025825.1_79461151 TGG                                                      
    ## # ... with 427 more rows

``` r
# loci
genos_long %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique() %>%
  group_by(locus) %>%
  tally() %>%
  arrange(desc(n))
```

    ## # A tibble: 140 x 2
    ##    locus                   n
    ##    <chr>               <int>
    ##  1 CM025832.1_31942157    10
    ##  2 CM025827.1_14090856     9
    ##  3 CM025846.1_17453041     8
    ##  4 CM025832.1_13857363     7
    ##  5 CM025841.1_24413058     7
    ##  6 CM025844.1_8337516      7
    ##  7 CM025847.1_4559639      7
    ##  8 CM025827.1_13614142     6
    ##  9 CM025844.1_11631411     6
    ## 10 CM025844.1_17570547     6
    ## # ... with 130 more rows

437 alleles across 140 loci with between 1-10 alleles per locus.

Missing data:

how many indiv?

``` r
genos_long %>%
  select(gtseq_run, id) %>%
  unique
```

    ## # A tibble: 286 x 2
    ##    gtseq_run id   
    ##    <chr>     <chr>
    ##  1 POPMix    s1   
    ##  2 POPMix    s10  
    ##  3 POPMix    s100 
    ##  4 POPMix    s101 
    ##  5 POPMix    s102 
    ##  6 POPMix    s103 
    ##  7 POPMix    s104 
    ##  8 POPMix    s105 
    ##  9 POPMix    s106 
    ## 10 POPMix    s107 
    ## # ... with 276 more rows

286

``` r
# missing data across loci
locs_to_toss <- genos_long %>%
  group_by(locus) %>%
  mutate(missingness = ifelse(is.na(allele), 1, 0)) %>%
  summarise(sum(missingness)) %>% # 140 loci x2 x192 indiv = total
  filter(`sum(missingness)`>286) %>% # more than 50% missing data
  select(locus) # drop those loci for now and see how the assignment goes

# just the keepers
genos_locs_filtered <- genos_long %>%
  anti_join(., locs_to_toss)
```

    ## Joining, by = "locus"

``` r
# summary of remaining loci
genos_locs_filtered %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique()
```

    ## # A tibble: 422 x 2
    ##    locus               allele                                                   
    ##    <chr>               <chr>                                                    
    ##  1 CM025825.1_24578215 ATTG                                                     
    ##  2 CM025825.1_26537970 CCCA                                                     
    ##  3 CM025825.1_43208600 GAGGCTGTCCATGTTCAATCTCGTAGAG                             
    ##  4 CM025825.1_43208600 GAGGCTGTCCATGTTCAATCACGTAGAG                             
    ##  5 CM025825.1_61422707 CCTCACTGCACGGGGGCGCCCTCTGGTTGCAAAAAGAAGTCTGATTGATAGAAGTC~
    ##  6 CM025825.1_61422707 CCTCACTGTACGGGGGCGCCCTCTGGTTGCAAAAAGAAGTCTGATTGATAGAAGTC~
    ##  7 CM025825.1_76785270 ACAGG                                                    
    ##  8 CM025825.1_7920395  TCCTCCT                                                  
    ##  9 CM025825.1_7920395  TCCTACT                                                  
    ## 10 CM025825.1_79461151 TGG                                                      
    ## # ... with 412 more rows

``` r
genos_locs_filtered %>%
  filter(!is.na(allele)) %>%
  select(locus, allele) %>%
  unique() %>%
  group_by(locus) %>%
  tally() %>%
  arrange(desc(n))
```

    ## # A tibble: 133 x 2
    ##    locus                   n
    ##    <chr>               <int>
    ##  1 CM025832.1_31942157    10
    ##  2 CM025827.1_14090856     9
    ##  3 CM025846.1_17453041     8
    ##  4 CM025832.1_13857363     7
    ##  5 CM025841.1_24413058     7
    ##  6 CM025844.1_8337516      7
    ##  7 CM025847.1_4559639      7
    ##  8 CM025827.1_13614142     6
    ##  9 CM025844.1_11631411     6
    ## 10 CM025844.1_17570547     6
    ## # ... with 123 more rows

``` r
# QC check
genos_long %>%
  group_by(id, locus) %>% # there should be no more than 2 alleles for a given indiv/locus
  tally() %>%
  filter(n > 2)
```

    ## # A tibble: 13,254 x 3
    ## # Groups:   id [94]
    ##    id    locus                   n
    ##    <chr> <chr>               <int>
    ##  1 s1    CM025825.1_24578215     4
    ##  2 s1    CM025825.1_26537970     4
    ##  3 s1    CM025825.1_43208600     4
    ##  4 s1    CM025825.1_61422707     4
    ##  5 s1    CM025825.1_76785270     4
    ##  6 s1    CM025825.1_7920395      4
    ##  7 s1    CM025825.1_79461151     4
    ##  8 s1    CM025826.1_15801073     4
    ##  9 s1    CM025826.1_22563832     4
    ## 10 s1    CM025826.1_36818699     4
    ## # ... with 13,244 more rows

## Missing data in individuals

Total number of loci = 133 Total number of gene copies = 266

Total number of samples = 286

``` r
# 25% missing data
266*0.25
```

    ## [1] 66.5

``` r
inds_to_toss <- genos_locs_filtered %>%
  group_by(id) %>%
  mutate(missingness = ifelse(is.na(allele), 1, 0)) %>%
  summarise(sum(missingness)) %>%
  arrange(desc(`sum(missingness)`)) %>%
  filter(`sum(missingness)` > 67) # remove samples with >25% missing data

# just the keepers
genos_locs_ind_filtered <- genos_locs_filtered %>%
  anti_join(., inds_to_toss)
```

    ## Joining, by = "id"

How many indivs retained?

``` r
genos_locs_ind_filtered %>%
  select(gtseq_run, id) %>%
  unique()
```

    ## # A tibble: 262 x 2
    ##    gtseq_run id   
    ##    <chr>     <chr>
    ##  1 POPMix    s1   
    ##  2 POPMix    s10  
    ##  3 POPMix    s100 
    ##  4 POPMix    s101 
    ##  5 POPMix    s102 
    ##  6 POPMix    s103 
    ##  7 POPMix    s104 
    ##  8 POPMix    s105 
    ##  9 POPMix    s106 
    ## 10 POPMix    s107 
    ## # ... with 252 more rows

Take a look at that dataset

``` r
genos_locs_ind_filtered  %>%
  unite(gtseq_run, id, col = "sample", remove = F) %>%
  ggplot(aes(x = reorder(locus, depth), y = reorder(sample, depth), fill = log10(depth))) +
  geom_tile()
```

![](02-POPMix-plates_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave("pdf_outputs/POPheatmap_POPmix.pdf")
```

    ## 

    ## Saving 7 x 5 in image

That’s a strange pattern - with NAs for some individuals in some of the
higher read depth loci. Variation in the primer sites across POP groups?

Quick look at a PCA before getting the expected population information
for POP groups/ABLG numbers.

## Taking data through to a PCA

If that’s the case, just take a peak at a PCA for these filtered
genotypes to evaluate whether there’s any obvious structure among these
bycatch samples.

The format I need for the PCA is strata, indiv, markers, genotype
(allele idx)

``` r
# first make integers of the alleles
alle_idxs <- genos_locs_ind_filtered %>% 
  unite(gtseq_run, id, col = "sample", remove = F) %>%
  dplyr::select(sample, locus, gene_copy, allele) %>%
  group_by(locus) %>%
  mutate(alleidx = as.integer(factor(allele, levels = unique(allele)))) %>%
  ungroup() %>%
  arrange(sample, locus, alleidx) # rubias can handle NA's, so no need to change them to 0's
```

Make the df match the requirements for tidy_genomic_data

``` r
long_df <- alle_idxs %>%
  select(-allele, -gene_copy) %>%
  mutate(group = "POP") %>%
  #left_join(., spp_indiv) %>%
  #select(species, everything()) %>%
  rename(INDIVIDUALS = sample, STRATA = group, MARKERS = locus, GT = alleidx)
```

Read in some metadata to take a closer look at those
samples/assignments.

``` r
POPmix_samplesheet <- read_csv("../data/samplesheets/20230308_POPMix_SebastesLarvae.csv", skip = 19)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   I7_Index_ID = col_character(),
    ##   index = col_character(),
    ##   I5_Index_ID = col_character(),
    ##   index2 = col_character(),
    ##   Sample_Project = col_character()
    ## )

``` r
meta_spp <- read_xlsx("../../rockfish-species-id/data_AFSC/POP_northerns_duskies_to_gtseq.xlsx")

# fix the metadata - apparently it was completely wrong.
meta_fixed <- meta_spp %>%
  #filter(SpeciesName %in% c("Sebastes ciliatus", "Sebastes variabilis")) %>%
  mutate(SpeciesName = ifelse(str_detect(AlternateID_s_, "Wes"), "Sebastes variabilis", SpeciesName)) %>%
  mutate(SpeciesName = ifelse(str_detect(AlternateID_s_, "LDRYAK97"), "Sebastes variabilis", SpeciesName)) %>%
    mutate(SpeciesName = ifelse(str_detect(AlternateID_s_, "DLRKOD97"), "Sebastes variabilis", SpeciesName)) 

meta_fixed %>%
  write_csv("../../rockfish-species-id/data_AFSC/POP_northerns_duskies_to_gtseq_FIXEDmetadata20230327.csv")

meta_fixed$ABLG <- as.character(meta_fixed$ABLG)

# join those to get the species
meta_for_gtseqMix <- POPmix_samplesheet %>%
  filter(Sample_Project == "POPMix") %>%
       left_join(., meta_fixed, by = c("Sample_ID" = "ABLG")) %>%
  select(Sample_ID, SpeciesName) %>%
  mutate(id = row_number()) %>%
  mutate(s = "POPMix_s") %>%
  unite(col = "s_num", s, id, sep = "") %>%
  rename(grp = SpeciesName)

meta_for_gtseqMix %>%
  group_by(grp) %>%
  tally()
```

    ## # A tibble: 5 x 2
    ##   grp                     n
    ##   <chr>               <int>
    ## 1 Sebastes alutus        48
    ## 2 Sebastes ciliatus      30
    ## 3 Sebastes polyspinis    48
    ## 4 Sebastes variabilis    64
    ## 5 <NA>                    2

``` r
# read in metadata for POP gtseq3 run with group designations
meta1 <- read_csv("../data/pop_335_population_assignment.csv")
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   ABLG = col_character(),
    ##   pop_assign = col_character(),
    ##   pop = col_character()
    ## )

``` r
samplelist <- read_csv("../data/samplesheets/20230131POPtest4_BFALtest4_KGLibs.csv", skip = 19)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Sample_ID = col_character(),
    ##   Sample_Plate = col_character(),
    ##   Sample_Well = col_character(),
    ##   I7_Index_ID = col_character(),
    ##   index = col_character(),
    ##   I5_Index_ID = col_character(),
    ##   index2 = col_character(),
    ##   Sample_Project = col_character(),
    ##   Description = col_character(),
    ##   s_id = col_character()
    ## )

``` r
metadata <- samplelist %>%
  inner_join(., meta1, by = c("Sample_ID" = "ABLG")) %>%
  mutate(grp = pop_assign) %>%
  select(grp, Sample_ID, s_id) %>%
  mutate(new_id = "POP_") %>%
  unite(col = "s_num", new_id, s_id, sep = "") %>%
  select(Sample_ID, grp, s_num)

all_meta <- metadata %>% 
  bind_rows(., meta_for_gtseqMix)
```

``` r
tmp <- long_df %>%
  left_join(., all_meta, by = c("INDIVIDUALS" = "s_num"))

long_df_w_strata <- tmp %>%
  mutate(STRATA = grp) %>%
  select(STRATA, everything()) %>%
  select(-grp, -Sample_ID)
```

``` r
convert_to_genind_pca(long_df_w_strata, "POPmix_POPgtseq3_comboPCA")
```

    ## `summarise()` has grouped output by 'INDIVIDUALS'. You can override using the
    ## `.groups` argument.Joining, by = c("INDIVIDUALS", "MARKERS")

    ## png 
    ##   2

## Duskies

``` r
# just the duskies
duskies <- long_df_w_strata %>%
  filter(STRATA %in% c("Sebastes variabilis", "Sebastes ciliatus"))

# genind conversion and pca plotting
convert_to_genind_pca(duskies, "duskiesPCA")
```

    ## `summarise()` has grouped output by 'INDIVIDUALS'. You can override using the
    ## `.groups` argument.Joining, by = c("INDIVIDUALS", "MARKERS")

    ## png 
    ##   2

## just the POP groups

``` r
# just the POP
pop_groups <- long_df_w_strata %>%
  filter(STRATA %in% c("A", "B", "C", "D", "B1", "B2", "B3", "B4"))

# genind conversion and pca plotting
convert_to_genind_pca(pop_groups, "POPgroups_PCA")
```

    ## `summarise()` has grouped output by 'INDIVIDUALS'. You can override using the
    ## `.groups` argument.Joining, by = c("INDIVIDUALS", "MARKERS")

    ## png 
    ##   2

### Notes

A few notes about this process:

Remember when merging VCF files that the sample numbers (s1, s100)
cannot be identical across gtseq runs. This requires renaming the sam
files and remaking the bam files for the VCF.

Merge VCF files on Sedna:

/home/dbaetscher/workdir/POP/gtseq/POP_gtseq3_filtered.vcf
/home/dbaetscher/workdir/POP/gtseq/POPMix_gtseq_filtered.vcf

Using some helpful tips from ECA:
<https://eriqande.github.io/eca-bioinf-handbook/handle-vcf.html>

Check out some space/time:

``` sh
srun -c 2 -t 03:00:00 --pty /bin/bash
```

On Sedna, I need to load some modules: bio/bcftools/1.11
lib64/htslib/1.11

before merging, I need to zip the vcf files (bgzip file.vcf) \# merging
requires that the files be indexed bcftools index first3.vcf.gz bcftools
index last3.vcf.gz

error with the POPMix file, need to sort first bcftools sort -O z
–output-file POPMix.vcf.gz POPMix_gtseq_filtered.vcf.gz

Now merge: bcftools merge -Oz POPMix.vcf.gz POP_gtseq3_filtered.vcf.gz
\> POPmerged.vcf.gz

ugh. duplicate sample names. Need to rename the sam files, then remake
the vcf files and then sort, index, and merge.

bcftools sort -O z –output-file POPgtseq3.vcf.gz
POP_gtseq3_renamed_filtered.vcf
