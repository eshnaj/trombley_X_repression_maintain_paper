library(tidyr)
library(ggplot2)
library(data.table)

#WT, dpy-21, cec-4, dpy-21; cec-4
all_fpkm <- fpkm(DErun)
d21_fpkm <- all_fpkm[,c("9501-EJ-1",
                        "9501-EJ-2",
                        "9501-EJ-3",
                        "9501-EJ-4",
                        "9501-EJ-5",
                        "9501-EJ-6",
                        "9501-EJ-10",
                        "9501-EJ-11",
                        "9501-EJ-12",
                        "9501-EJ-13",
                        "9501-EJ-14",
                        "9501-EJ-15")]
summary(rowSums(d21_fpkm[,c(1:12)]) == 0)
d21_fpkm <- d21_fpkm[rowSums(d21_fpkm[,c(1:12)]) > 0,]
summary(rowSums(d21_fpkm[,c(1:12)]) == 0)
d21_fpkm <- t(d21_fpkm)

PCA_data <- prcomp(d21_fpkm, center = TRUE)
plot_PCA_data <- PCA_data[["x"]]
plot_PCA_data <- as.data.frame(plot_PCA_data)
plot_PCA_data$label <- rownames(plot_PCA_data)
plot_PCA_data$genotype <- c(rep("N2", 3), 
                            rep("cec4", 3), 
                            rep("dpy21", 3), 
                            rep("dpy21_cec4", 3) )

summary(PCA_data)

dir.create("PCA_plots")

ggplot(plot_PCA_data, 
       aes(x = PC1, y = PC2, color = genotype)) + 
  geom_point(size=3) + 
  labs(x = paste0("PC1 (", summary(PCA_data)$importance[2,][1]*100, "% variance)"),
       y = paste0("PC2 (", summary(PCA_data)$importance[2,][2]*100, "% variance)")) +
  scale_color_manual(values = c("#007991","#439A86","#BCD8C1","#222E50" )) +
  theme(title = element_text(size = 20), 
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18), 
        axis.ticks = element_blank(), 
        axis.text = element_blank())

ggsave("PCA_plots/PC1_PC2_d21.png", 
       height = 4, 
       width = 5)

#auxin-treated samples
d27aux_fpkm <- all_fpkm[,c("9501-EJ-25",
                                              "9501-EJ-26",
                                              "9501-EJ-36",
                                              "9501-EJ-23",
                                              "9501-EJ-24",
                                              "9501-EJ-37",
                                              "9501-EJ-31",
                                              "9501-EJ-32",
                                              "9501-EJ-33",
                                              "9501-EJ-28",
                                              "9501-EJ-29",
                                              "9501-EJ-30")]
summary(rowSums(d27aux_fpkm) == 0)
d27aux_fpkm <- d27aux_fpkm[rowSums(d27aux_fpkm) > 0,]
summary(rowSums(d27aux_fpkm) == 0)

d27aux_fpkm <- t(d27aux_fpkm)

PCA10_data <- prcomp(d27aux_fpkm, center = TRUE, scale. = TRUE)
plot_PCA10_data <- PCA10_data[["x"]]
plot_PCA10_data <- as.data.frame(plot_PCA10_data)
plot_PCA10_data$label <- rownames(plot_PCA10_data)
plot_PCA10_data$genotype <- c(rep("tir1_dpy27AID_aux", 3), 
                              rep("tir1_dpy27AID_noaux", 3),
                              rep("tir1_dpy27AID_cec4_aux", 3),
                              rep("tir_dpy27AID_cec4_noaux", 3))

summary(PCA10_data)

#version for publication
ggplot(plot_PCA10_data, aes(x = PC1, y = PC2, color = genotype)) + 
  geom_point(size=3) + 
  labs(x = paste0("PC1 (", summary(PCA10_data)$importance[2,][1]*100, "% variance)"), 
       y = paste0("PC2 (", summary(PCA10_data)$importance[2,][2]*100, "% variance)")) + 
  scale_color_manual(values = c("#820A5C","#EC5027","#D24EA8","#B02E0C" )) +
  theme(title = element_text(size = 20), 
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18),
        axis.ticks = element_blank(),
        axis.text = element_blank())

ggsave("d27aux_PC1_PC2.png", 
       height = 4, 
       width =6.5)