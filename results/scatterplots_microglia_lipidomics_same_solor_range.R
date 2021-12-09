library(ggplot2)
library(ggrepel)
library(tidyverse)

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

###
### Which range of log FC should be used for the gradient?
###

###cl_01

cl_01_sig_p <- filter(cl_01, cl_01$FDR<0.1)

##Log10 scale
mid<-0
high <- filter(cl_01_sig_p, cl_01_sig_p$logFC>1)
low <- filter(cl_01_sig_p, -1 > cl_01_sig_p$logFC)
mid1 <- filter(cl_01_sig_p, 1 > cl_01_sig_p$logFC & cl_01_sig_p$logFC >= -1 )
cl_01 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=mid1, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-1, 1), 
                                             breaks=c( -1, -0.5, 0, 0.5, 1),
                                             labels = c(-1, -0.5, 0, 0.5, 1))+
  geom_point(data=high, alpha=1, color="red", size = 3) +
  geom_point(data=low, alpha=1, color="blue", size = 3) +
  theme_bw()+
  ggtitle("Cell Lipids - 01 hr") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = cl_01_sig_p, aes(label = lipid))

###cl_12

cl_12_sig_p <- filter(cl_12, cl_12$FDR<0.1)

##Log10 scale
mid<-0
high <- filter(cl_12_sig_p, cl_12_sig_p$logFC>1)
low <- filter(cl_12_sig_p, -1 > cl_12_sig_p$logFC)
mid1 <- filter(cl_12_sig_p, 1 > cl_12_sig_p$logFC & cl_12_sig_p$logFC >= -1 )
cl_12 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=mid1, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-1, 1), 
                                             breaks=c( -1, -0.5, 0, 0.5, 1),
                                             labels = c(-1, -0.5, 0, 0.5, 1))+
  geom_point(data=high, alpha=1, color="red", size = 3) +
  geom_point(data=low, alpha=1, color="blue", size = 3) +
  theme_bw()+
  ggtitle("Cell Lipids - 12 hr") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = cl_12_sig_p, aes(label = lipid))


###cl_24

cl_24_sig_p <- filter(cl_24, cl_24$FDR<0.1)

##Log10 scale
mid<-0
high <- filter(cl_24_sig_p, cl_24_sig_p$logFC>1)
low <- filter(cl_24_sig_p, -1 > cl_24_sig_p$logFC)
mid1 <- filter(cl_24_sig_p, 1 > cl_24_sig_p$logFC & cl_24_sig_p$logFC >= -1 )
cl_24 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=mid1, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-1, 1), 
                                             breaks=c( -1, -0.5, 0, 0.5, 1),
                                             labels = c(-1, -0.5, 0, 0.5, 1))+
  geom_point(data=high, alpha=1, color="red", size = 3) +
  geom_point(data=low, alpha=1, color="blue", size = 3) +
  theme_bw()+
  ggtitle("Cell Lipids - 24 hr") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = cl_24_sig_p, aes(label = lipid))

###cm_01

cm_01_sig_p <- filter(cm_01, cm_01$FDR<0.1)

##Log10 scale
mid<-0
cm_01 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=cm_01_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-2.2, 1.2), 
                                             breaks=c(-2, -1.5, -1, -0.5, 0, 0.5, 1),
                                             labels = c(-2, -1.5, -1, -0.5, 0, 0.5, 1))+
  theme_bw()+
  ggtitle("Cell Metabolites - 01 hr") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = cm_01_sig_p, aes(label = lipid)) 

###cm_12

cm_12_sig_p <- filter(cm_12, cm_12$FDR<0.1)

##Log10 scale
mid<-0.3
cm_12 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=cm_12_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-0.6, 1.2), 
                                             breaks=c(-0.5, 0, 0.5, 1),
                                             labels = c(-0.5, 0, 0.5, 1))+
  theme_bw()+
  ggtitle("Cell Metabolites - 12 hrs") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = cm_12_sig_p, aes(label = lipid)) 


###cm_24

cm_24_sig_p <- filter(cm_24, cm_24$FDR<0.1)

##Log10 scale
mid<-0.5
cm_24 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=cm_24_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-0.5, 1.5), 
                                             breaks=c(-0.5, 0, 0.5, 1, 1.5),
                                             labels = c(-0.5, 0, 0.5, 1, 1.5))+
  theme_bw()+
  ggtitle("Cell Metabolites - 24 hrs") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = cm_24_sig_p, aes(label = lipid)) 


###ml_01

ml_01_sig_p <- filter(ml_01, ml_01$FDR<0.1)

##Log10 scale
mid<-0
ml_01 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=ml_01_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-1.7, 0.6), 
                                             breaks=c(-1.5, -1, -0.5, 0, 0.5),
                                             labels = c(-1.5, -1, -0.5, 0, 0.5))+
  theme_bw()+
  ggtitle("Media Lipids - 01 hr") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = ml_01_sig_p, aes(label = lipid)) 


###ml_12

ml_12_sig_p <- filter(ml_12, ml_12$FDR<0.1)

##Log10 scale
mid<-0
ml_12 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=ml_12_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-0.7, 0.7), 
                                             breaks=c(-0.5, 0, 0.5),
                                             labels = c(-0.5, 0, 0.5))+
  theme_bw()+
  ggtitle("Media Lipids - 12 hrs") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = ml_12_sig_p, aes(label = lipid)) 


###ml_24

ml_24_sig_p <- filter(ml_24, ml_24$FDR<0.1)

##Log10 scale
mid<-0
ml_24 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=ml_24_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-1, 1), 
                                             breaks=c(-1, 0, 1),
                                             labels = c(-1, 0, 1))+
  theme_bw()+
  ggtitle("Media Lipids - 24 hr") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = ml_24_sig_p, aes(label = lipid))

###mm_01

mm_01_sig_p <- filter(mm_01, mm_01$FDR<0.1)

##Log10 scale
mid<- -2.5
mm_01 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=mm_01_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-5, 0), 
                                             breaks=c(-5, -2.5, 0),
                                             labels = c(-5, -2.5, 0))+
  theme_bw()+
  ggtitle("Media Metabolites - 01 hr") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = mm_01_sig_p, aes(label = lipid))


###mm_12

mm_12_sig_p <- filter(mm_12, mm_12$FDR<0.1)

##Log10 scale
mid<-0.5
mm_12 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=mm_12_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-0.5, 1.5), 
                                             breaks=c(-0.5, 0, 0.5, 1, 1.5),
                                             labels = c(-0.5, 0, 0.5, 1, 1.5))+
  theme_bw()+
  ggtitle("Media Metabolites - 12 hrs") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = mm_12_sig_p, aes(label = lipid))


###mm_24

mm_24_sig_p <- filter(mm_24, mm_24$FDR<0.1)

##Log10 scale
mid<-0.5
mm_24 %>%
  ggplot(aes(x=log10(abe_average),y=log10(veh_average), color=logFC)) + 
  geom_point(alpha=1, color="darkgrey", size = 3) +
  
  geom_point(data=mm_24_sig_p, 
             aes(x=log10(abe_average),y=log10(veh_average)),
             size=3) + scale_color_gradient2(midpoint=mid, low="blue",mid="#DC9313",
                                             high="red", space ="Lab", limits = c(-0.5, 1.5), 
                                             breaks=c(-0.5, 0, 0.5, 1, 1.5),
                                             labels = c(-0.5, 0, 0.5, 1, 1.5))+
  theme_bw()+
  ggtitle("Media Metabolites - 24 hrs") +
  theme(text = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggrepel::geom_text_repel(data = mm_24_sig_p, aes(label = lipid))