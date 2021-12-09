library(plyr)
library(tidyverse)
library(RColorBrewer)

cl_01 <- read.csv(file = 'cl_01h_full.csv')
cl_12 <- read.csv(file = 'cl_12h_full.csv')
cl_24 <- read.csv(file = 'cl_24h_full.csv')

cm_01 <- read.csv(file = 'cm_01h_full.csv')
cm_12 <- read.csv(file = 'cm_12h_full.csv')
cm_24 <- read.csv(file = 'cm_24h_full.csv')

ml_01 <- read.csv(file = 'ml_01h_full.csv')
ml_12 <- read.csv(file = 'ml_12h_full.csv')
ml_24 <- read.csv(file = 'ml_24h_full.csv')

mm_01 <- read.csv(file = 'mm_01h_full.csv')
mm_12 <- read.csv(file = 'mm_12h_full.csv')
mm_24 <- read.csv(file = 'mm_24h_full.csv')


################## CL 01 hr ###################################
df <- cl_01

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6,
         -veh_N7,	-Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Cell Lipids - 01 hr", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     breaks = seq(-2.5, 1.5, length.out = 100),
                     legend_breaks = c(-2, -1, 0, 1), 
                     legend_labels = c(-2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 10)

################## CL 12 hr ###################################
df <- cl_12

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6, -veh_N1, -veh_N2,
         -Abe_N1, -Abe_N2) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Cell Lipids - 12 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     breaks = seq(-2.2, 1.3, length.out = 100),
                     legend_breaks = c(-2, -1, 0, 1), 
                     legend_labels = c(-2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 10)


################## CL 24 hrs ###################################
df <- cl_24

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6, 
         -veh_N2, -veh_N7,
         -Abe_N2, -Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Cell Lipids - 24 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     breaks = seq(-2, 1, length.out = 100),
                     legend_breaks = c(-2, -1, 0, 1), 
                     legend_labels = c(-2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 10)


################## CM 01 hr ###################################
df <- cm_01

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6,
         -veh_N7,	-Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Cell Metabolites - 01 hr", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     # breaks = seq(-2.5, 1.5, length.out = 100),
                     legend_breaks = c(-2, -1, 0, 1), 
                     legend_labels = c(-2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 10)

################## CM 12 hr ###################################
df <- cm_12

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6,
         -veh_N1, -veh_N2, -Abe_N1, -Abe_N2) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Cell Metabolites - 12 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     # breaks = seq(-2.5, 1.5, length.out = 100),
                     legend_breaks = c(-1, 0, 1), 
                     legend_labels = c(-1, 0, 1),
                     fontsize = 10, fontsize_row = 10)

################## CM 24 hrs ###################################
df <- cm_24

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6, 
         -veh_N2, -veh_N7,
         -Abe_N2, -Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Cell Metabolites - 24 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     # breaks = seq(-2, 1, length.out = 100),
                     legend_breaks = c( -1, 0, 1,2), 
                     legend_labels = c( -1, 0, 1,2),
                     fontsize = 10, fontsize_row = 10)


################## ML 01 hr ###################################
df <- ml_01

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6,
         -veh_N7,	-Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

##### No siginificant lipids

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Media Lipids - 01 hr", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     # breaks = seq(-2.5, 1.5, length.out = 100),
                     legend_breaks = c(-2, -1, 0, 1), 
                     legend_labels = c(-2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 10)

################## ML 12 hr ###################################
df <- ml_12

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6,
         -veh_N1, -veh_N2, -Abe_N1, -Abe_N2) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Media Lipids - 12 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE)#, 
                     # breaks = seq(-2.5, 1.5, length.out = 100),
                     # legend_breaks = c(-1, 0, 1), 
                     # legend_labels = c(-1, 0, 1),
                     # fontsize = 10, fontsize_row = 10)

################## ML 24 hrs ###################################
df <- ml_24

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6, 
         -veh_N2, -veh_N7,
         -Abe_N2, -Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Media Lipids - 24 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     # breaks = seq(-2, 1, length.out = 100),
                     legend_breaks = c( -1, 0, 1,2), 
                     legend_labels = c( -1, 0, 1,2),
                     fontsize = 10, fontsize_row = 10)

################## MM 01 hr ###################################
df <- mm_01

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6,
         -veh_N7,	-Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Media Metabolites - 01 hr", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE)#, 
                     # # breaks = seq(-2.5, 1.5, length.out = 100),
                     # legend_breaks = c(-2, -1, 0, 1), 
                     # legend_labels = c(-2, -1, 0, 1),
                     # fontsize = 10, fontsize_row = 10)

################## MM 12 hr ###################################
df <- mm_12

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6,
         -veh_N1, -veh_N2, -Abe_N1, -Abe_N2) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Media Metabolites - 12 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
# breaks = seq(-2.5, 1.5, length.out = 100),
                     legend_breaks = c(-1, 0, 1),
                     legend_labels = c(-1, 0, 1),
                     fontsize = 10, fontsize_row = 10)


################## MM 24 hrs ###################################
df <- mm_24

df %>% mutate(lipids = make.unique(lipid)) %>%
  filter(FDR<0.1) %>%
  select(-lipid,	-logFC,	-logCPM,	-LR,	-PValue,	-FDR,	-Transition,	-type,	-Blank,
         -veh_N3,	-Abe_N3,	-veh_N4,	-Abe_N4,	-veh_N5,	-Abe_N5,	-veh_N6,	-Abe_N6, 
         -veh_N2, -veh_N7,
         -Abe_N2, -Abe_N7) %>%
  rowwise() %>%
  mutate(mean = mean(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, log10) %>%
  mutate_if(is.numeric, list(~ . - mean)) %>%
  mutate(sd = sd(c(abe_average, veh_average))) %>%
  mutate_if(is.numeric, list(~ . / sd)) %>%
  select( -mean, -sd) ->
  vals

vals %>%
  column_to_rownames("lipids") %>%
  as.matrix() %>%
  pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                               "RdBu")))(100), 
                     main = "Media Lipids - 24 hrs", 
                     cluster_rows = TRUE,
                     cluster_cols = FALSE, 
                     # breaks = seq(-2, 1, length.out = 100),
                     legend_breaks = c( -2, -1, 0, 1), 
                     legend_labels = c( -2, -1, 0, 1),
                     fontsize = 10, fontsize_row = 10)