library(plyr)
library(tidyverse)
library(readxl)

hour_01_names <- c("N3", "N4", "N5", "N6", "N7")
hour_12_names <- c("N1", "N2", "N3", "N4", "N5", "N6")
hour_24_names <- c("N2", "N3", "N4", "N5", "N6", "N7")

lipid_types <- c("acy-carnitines", "CE", "cer", "FFA", "PCandSM",
                 "PE", "PG", "PI", "PS", "QC", "TAG1", "TAG2")

metob_types <- c("neg1", "neg2", "pos1", "pos2")

read_lipids <- function(prefix, suffix, tp, tp_names, types, incl_suffix = F) {
  lapply(types, function(type){
    
    all_names <- lapply(tp_names, function(n) {
      veh_name = paste0("veh_", n)
      
      veh = read_xlsx(paste0(prefix, "/", type, "/", "Veh-", suffix,
                             tp, n, "_", type, "_results.xlsx")) %>%
        rename(lipid = `Lipid name`) %>%
        rename(!!veh_name := `Total Intensity`) %>%
        mutate(type = type) %>%
        na.omit
      
      abe_name = paste0("Abe_", n)
      
      abe = read_xlsx(paste0(prefix, "/", type, "/", "Ab-", suffix,
                             tp, n, "_", type, "_results.xlsx")) %>%
        rename(lipid = `Lipid name`) %>%
        rename(!!abe_name := `Total Intensity`) %>%
        mutate(type = type) %>%
        na.omit
      
      blank = read_xlsx(paste0(prefix, "/", type, "/", "Blank_",
                               if_else(incl_suffix, paste0(suffix, "_"),
                                       ""),
                               type, "_results.xlsx")) %>%
        rename(lipid = `Lipid name`) %>%
        rename(Blank = `Total Intensity`)
      
      merge(veh, abe) %>% merge(blank) %>% as_tibble
      
    })
    
    Reduce(function(x, y) merge(x, y), all_names)
  }) %>% bind_rows() %>% as_tibble()
}

cl_01h <- read_lipids("raw_data/cells-lipids", "cl", "_1h_",
                      hour_01_names, lipid_types)
cl_12h <- read_lipids("raw_data/cells-lipids", "cl", "_12h_",
                      hour_12_names, lipid_types)
cl_24h <- read_lipids("raw_data/cells-lipids", "cl", "_24h_",
                      hour_24_names, lipid_types)

ml_01h <- read_lipids("raw_data/media-lipids", "medium", "_1h_",
                      hour_01_names, lipid_types)
ml_12h <- read_lipids("raw_data/media-lipids", "medium", "_12h_",
                      hour_12_names, lipid_types)
ml_24h <- read_lipids("raw_data/media-lipids", "medium", "_24h_",
                      hour_24_names, lipid_types)

cm_01h <- read_lipids("raw_data/polar", "cm", "_1h_",
                      hour_01_names, metob_types, T)
cm_12h <- read_lipids("raw_data/polar", "cm", "_12h_",
                      hour_12_names, metob_types, T)
cm_24h <- read_lipids("raw_data/polar", "cm", "_24h_",
                      hour_24_names, metob_types, T)

mm_01h <- read_lipids("raw_data/polar", "mm", "_1h_",
                      hour_01_names, metob_types, T)
mm_12h <- read_lipids("raw_data/polar", "mm", "_12h_",
                      hour_12_names, metob_types, T)
mm_24h <- read_lipids("raw_data/polar", "mm", "_24h_",
                      hour_24_names, metob_types, T)


colnames(cl_24h) <- c("lipid","Transition", "type", "Blank",  "24veh_N2","24Abe_N2", "24veh_N3", "24Abe_N3", "24veh_N4", "24Abe_N4", "24veh_N5", "24Abe_N5", "24veh_N6", "24Abe_N6", "24veh_N7", "24Abe_N7")

### EdgeR analysis

## Define edgeR variables

library(edgeR)

# Groups for 01h
gr_01h = c("Blank",
           "VH", "AB", "VH", "AB", "VH", "AB",
           "VH", "AB", "VH", "AB") %>%
  factor(levels = c("Blank", "AB", "VH"))

design_01h = model.matrix(~gr_01h)

contrasts_01h = makeContrasts(
  H = gr_01hAB - gr_01hVH,
  levels = design_01h
)

# Groups for 12h
gr_12h = c("Blank",
           "VH", "AB", "VH", "AB", "VH", "AB",
           "VH", "AB", "VH", "AB", "VH", "AB") %>%
  factor(levels = c("Blank", "AB", "VH"))

design_12h = model.matrix(~gr_12h)

contrasts_12h = makeContrasts(
  H = gr_12hAB - gr_12hVH,
  levels = design_12h
)

# Groups for 24h
gr_24h = c("Blank",
           "VH", "AB", "VH", "AB", "VH", "AB",
           "VH", "AB", "VH", "AB", "VH", "AB") %>%
  factor(levels = c("Blank", "AB", "VH"))

design_24h = model.matrix(~gr_24h)

contrasts_24h = makeContrasts(
  H = gr_24hAB - gr_24hVH,
  levels = design_24h
)


## Define function for statistical significance and perform initial analysis


perform_analysis_raw <- function(counts, design_mat, gr) {
  
  data.edgeR <- DGEList(counts = counts %>%
                          na.omit %>%
                          mutate(lipid = make.unique(lipid)) %>%
                          select(-Transition, -type) %>%
                          column_to_rownames("lipid"),
                        group = gr
  )
  
  data.edgeR <- calcNormFactors(data.edgeR, method="TMM")
  data.edgeR <- estimateCommonDisp(data.edgeR, design=design_mat)
  data.edgeR
}

calculate_significance <- function(dge, contrast) {
  dge %>%
    glmFit() %>%
    glmLRT(contrast = contrast)
}



cl_01h %>%
  perform_analysis_raw(design_01h, gr_01h) %>%
  calculate_significance(contrasts_01h) %>%
  topTags(30) %>% as.data.frame() %>% as_tibble()

cl_24h %>%
  perform_analysis_raw(design_24h, gr_24h) %>%
  calculate_significance(contrasts_24h) %>%
  topTags(30) %>% as.data.frame() %>% as_tibble()


source("ggbiplot.R")

get_DE_lipids <- function(counts, design_mat, gr, contrasts, p.value = 0.1) {
  dls <-
    counts %>%
    perform_analysis_raw(design_mat, gr) %>%
    calculate_significance(contrasts) %>%
    decideTestsDGE(p.value = p.value)
  
  rownames(dls)[dls %>% as.logical()]
}


### 01 hr case #####

cl_01h_edger <-
  cl_01h %>%
  get_DE_lipids(design_01h, gr_01h, contrasts_01h)

if(length(cl_01h_edger) == 0) {
  cat("No significant lipids for ", title)
  return()
}

if(length(cl_01h_edger) == 1) {
  cat("Single significant lipids for ", title, " is ", cl_01h_edger[1])
  return()
}

cl_01h %>%
  na.omit %>%
  mutate(lipid = make.unique(lipid)) %>%
  filter(lipid %in% cl_01h_edger) %>%
  select(-Transition, -type, -Blank) %>%
  column_to_rownames("lipid") %>%
  as.matrix() %>%
  t -> df1h

  df1h%>%prcomp(center = T, scale = T) ->
  prcomp_data1

groups1 = NULL

cl_01h %>%
  select(-Transition, -type, -Blank, -lipid) %>%
  colnames() ->
  labels.tmp

groups1 = substr(labels.tmp, 1, 3)

if (!is.null(labels)) {
  labels1 = labels.tmp
}


###24 hr case #####

cl_1h_24h <- merge(cl_01h, cl_24h, all=TRUE, sort=FALSE)

# Groups for 1h
gr_01_24h = c("Blank",
           "VH1", "AB1", "VH1", "AB1", "VH1", "AB1",
           "VH1", "AB1", "VH1", "AB1",
           "VH24", "AB24", "VH24", "AB24", "VH24", "AB24",
           "VH24", "AB24", "VH24", "AB24", "VH24", "AB24") %>%
  factor(levels = c("Blank", "AB1", "VH1", "AB24", "VH24"))

design_01_24h = model.matrix(~gr_01_24h)

contrasts_01_24h = makeContrasts(
  H = gr_01_24hAB1 - gr_01_24hVH1, 
  H = gr_01_24hAB24 - gr_01_24hVH24,
  levels = design_01_24h
)

cl_1h_24h %>%
  perform_analysis_raw(design_01_24h, gr_01_24h) %>%
  calculate_significance(contrasts_01_24h) %>%
  topTags(30) %>% as.data.frame() %>% as_tibble()


perform_analysis_raw <- function(counts, design_mat, gr) {
  
  data.edgeR <- DGEList(counts = counts %>%
                          na.omit %>%
                          mutate(lipid = make.unique(lipid)) %>%
                          select(-Transition, -type) %>%
                          column_to_rownames("lipid"),
                        group = gr
  )
  
  data.edgeR <- calcNormFactors(data.edgeR, method="TMM")
  data.edgeR <- estimateCommonDisp(data.edgeR, design=design_mat)
  data.edgeR
}

calculate_significance <- function(dge, contrast) {
  dge %>%
    glmFit() %>%
    glmLRT(contrast = contrast)
}


get_DE_lipids <- function(counts, design_mat, gr, contrasts, p.value = 0.1) {
  dls <-
    counts %>%
    perform_analysis_raw(design_mat, gr) %>%
    calculate_significance(contrasts) %>%
    decideTestsDGE(p.value = p.value)
  
  rownames(dls)[dls %>% as.logical()]
}

cl_1h_24h_edger <-
  cl_1h_24h %>% get_DE_lipids(design_01_24h, gr_01_24h, contrasts_01_24h)

if(length(cl_1h_24h_edger) == 0) {
  cat("No significant lipids for ", title)
  return()
}

if(length(cl_1h_24h_edger) == 1) {
  cat("Single significant lipids for ", title, " is ", cl_1h_24h_edger[1])
  return()
}

cl_1h_24h %>%
  na.omit %>%
  mutate(lipid = make.unique(lipid)) %>%
  filter(lipid %in% cl_1h_24h_edger) %>%
  select(-Transition, -type, -Blank) %>%
  column_to_rownames("lipid") %>%
  as.matrix() %>%
  t -> df_1h_24h

df_1h_24h%>% prcomp(center = T, scale = T) -> prcomp_data_1h_24h

df4.groups <- c("1h_veh", "1h_Abe", "1h_veh", "1h_Abe", "1h_veh", "1h_Abe", "1h_veh", "1h_Abe", "1h_veh", "1h_Abe",
                "24h_veh", "24h_Abe","24h_veh", "24h_Abe", "24h_veh", "24h_Abe", "24h_veh", "24h_Abe", "24h_veh", "24h_Abe", "24h_veh", "24h_Abe")

df4.labels <- c("veh_N3", "Abe_N3", "veh_N4", "Abe_N4", "veh_N5", "Abe_N5", "veh_N6", "Abe_N6", "veh_N7", "Abe_N7",
                "24veh_N2","24Abe_N2", "24veh_N3", "24Abe_N3", "24veh_N4", "24Abe_N4", "24veh_N5", "24Abe_N5", "24veh_N6", "24Abe_N6", "24veh_N7", "24Abe_N7")


prcomp_data_1h_24h %>%
  ggbiplot(ellipse = TRUE,
           labels = df4.labels,
           groups = df4.groups,
           var.axes = FALSE
  )


# combined <- merge(cl_01h, cl_24h, all=TRUE, sort=FALSE)

###########################

###combined case #####

cl_24h_edger <-
  cl_24h %>%
  get_DE_lipids(design_24h, gr_24h, contrasts_24h)

if(length(cl_24h_edger) == 0) {
  cat("No significant lipids for ", title)
  return()
}

if(length(cl_24h_edger) == 1) {
  cat("Single significant lipids for ", title, " is ", cl_24h_edger[1])
  return()
}

cl_24h %>%
  na.omit %>%
  mutate(lipid = make.unique(lipid)) %>%
  filter(lipid %in% cl_24h_edger) %>%
  select(-Transition, -type, -Blank) %>%
  column_to_rownames("lipid") %>%
  as.matrix() %>%
  t -> df24h


##########################

combined %>% prcomp(center = T, scale = T) -> prcomp_combined

# write.csv(df1h, "df1h.csv")
# write.csv(df24h, "df24h.csv")

library("dplyr")
df1h <- as.data.frame(df1h)
df24h <- as.data.frame(df24h)

df4 <- bind_rows(df1h, df24h) 
df4 <- data.matrix(df4)



df3 <- read.csv("df1h_24h_common_rows.csv", row.names = 1)
df3 <- data.matrix(df3)

df4 %>% prcomp(center = T, scale = T)
%>% select_if(~ !any(is.na(.))) 
%>% prcomp(center = T, scale = T) -> prcomp_data_combined4
df4.groups <- c("1h_veh", "1h_Abe", "1h_veh", "1h_Abe", "1h_veh", "1h_Abe", "1h_veh", "1h_Abe", "1h_veh", "1h_Abe",
                "24h_veh", "24h_Abe","24h_veh", "24h_Abe", "24h_veh", "24h_Abe", "24h_veh", "24h_Abe", "24h_veh", "24h_Abe", "24h_veh", "24h_Abe")

df4.labels <- c("veh_N3", "Abe_N3", "veh_N4", "Abe_N4", "veh_N5", "Abe_N5", "veh_N6", "Abe_N6", "veh_N7", "Abe_N7",
                "24veh_N2","24Abe_N2", "24veh_N3", "24Abe_N3", "24veh_N4", "24Abe_N4", "24veh_N5", "24Abe_N5", "24veh_N6", "24Abe_N6", "24veh_N7", "24Abe_N7")

ggbiplot(prcomp_combined, ellipse = TRUE,
                       labels = df4.labels,
                       groups = df4.groups,
                       var.axes = FALSE)



groups2 = NULL

cl_24h %>%
  select(-Transition, -type, -Blank, -lipid) %>%
  colnames() ->
  labels.tmp

groups2 = substr(labels.tmp, 1, 3)

if (!is.null(labels)) {
  labels2 = labels.tmp
}








# df1 <- read.csv("ml_01h.csv")
# df2 <- read.csv("ml_24h.csv")
# df3 <- cbind(df1, df2)
# 
# df3.pca <- prcomp(df3[,c(10:19,30:38)], center = TRUE,scale. = TRUE)
# summary(df3.pca)
# 
# ggbiplot(df3.pca,ellipse = TRUE, labels=rownames(df3))

# a1<- ggbiplot(prcomp_data1,choices = 1:2, ellipse = TRUE,
#          labels = labels1,
#          groups = groups1,
#          var.axes = FALSE)
# 
# 
# a2<- ggbiplot(prcomp_data2, ellipse = TRUE,
#          labels = labels2,
#          groups = groups2,
#          var.axes = FALSE)

