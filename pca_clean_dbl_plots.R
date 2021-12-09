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

# Giving_separate column names for two combining dataframes

colnames(cl_01h) <- c("lipid","Transition", "type", "Blank", "01veh_N3", "01Abe_N3", "01veh_N4", "01Abe_N4", "01veh_N5", "01Abe_N5", "01veh_N6", "01Abe_N6", "01veh_N7", "01Abe_N7")
colnames(cl_24h) <- c("lipid","Transition", "type", "Blank",  "24veh_N2","24Abe_N2", "24veh_N3", "24Abe_N3", "24veh_N4", "24Abe_N4", "24veh_N5", "24Abe_N5", "24veh_N6", "24Abe_N6", "24veh_N7", "24Abe_N7")


#### Merging

cl_01_24h <- merge(cl_01h, cl_24h, all=TRUE, sort=FALSE)


library(edgeR)

# Groups for combined 1h and 24 hrs

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

write_summary_and_results <- function(tp, design, gr, contrasts, name) {
  tp %>%
    perform_analysis_raw(design, gr) %>%
    calculate_significance(contrasts) %>%
    topTags(30000) %>%
    as.data.frame() %>%
    rownames_to_column("lipid") %>%
    merge(tp) %>%
    as_tibble() %>%
    arrange(FDR) -> results
  
  write_csv(results, paste0("results/", name, "_full.csv"))
  
  results %>%
    group_by(type) %>%
    summarise(ab = sum(logFC < 0 & FDR < 0.1),
              ve = sum(logFC > 0 & FDR < 0.1)) %>%
    write_csv(paste0("results/", name, "_summary.csv"))
}

#### PCA for the combined 01 and 24 hrs


source("ggbiplot.R")

get_DE_lipids <- function(counts, design_mat, gr, contrasts, p.value = 0.1) {
  dls <-
    counts %>%
    perform_analysis_raw(design_mat, gr) %>%
    calculate_significance(contrasts) %>%
    decideTestsDGE(p.value = p.value)
  
  rownames(dls)[dls %>% as.logical()]
}

make_pca_plot <- function(tp, design_mat, gr, contrasts,
                          title = "PCA plot",
                          ellipse = T, var.axes = F,
                          labels = NULL) {
  
  tp_edger <-
    tp %>%
    get_DE_lipids(design_mat, gr, contrasts)
  
  if(length(tp_edger) == 0) {
    cat("No significant lipids for ", title)
    return()
  }
  
  if(length(tp_edger) == 1) {
    cat("Single significant lipids for ", title, " is ", tp_edger[1])
    return()
  }
  
  tp %>%
    na.omit %>%
    mutate(lipid = make.unique(lipid)) %>%
    filter(lipid %in% tp_edger) %>%
    select(-Transition, -type, -Blank) %>%
    column_to_rownames("lipid") %>%
    as.matrix() %>%
    t %>%
    prcomp(center = T, scale = T) ->
    prcomp_data
  
  groups = NULL
  
  tp %>%
    select(-Transition, -type, -Blank, -lipid) %>%
    colnames() ->
    labels.tmp
  
  groups = substr(labels.tmp, 1, 5)
  
  if (!is.null(labels)) {
    labels = labels.tmp
  }
  
  prcomp_data %>%
    ggbiplot(ellipse = ellipse,
             labels = labels,
             groups = groups,
             var.axes = var.axes
    ) +
    ggtitle(title) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect("transparent", "black"),
          legend.key = element_blank())
}


svglite::svglite("figures/primary_pca_cl_combined_01_and_24h.svg")
cl_01_24h %>%
  make_pca_plot(design_01_24h, gr_01_24h, contrasts_01_24h, "Cell Lipids at 01 Hr and 24 Hrs Together")
dev.off()


cl_01_24h %>%
  make_pca_plot(design_01_24h, gr_01_24h, contrasts_01_24h, "Cell Lipids at 01 Hr and 24 Hrs Together")
