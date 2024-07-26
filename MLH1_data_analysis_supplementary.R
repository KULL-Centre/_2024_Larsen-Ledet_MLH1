#Load in required packages
library(tidyverse)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)
library(ggExtra)
library(ggpubr)
library(ggpmisc)
library(hexbin)

#Change the path to your directory
path <- ""

#Read in dataframes
MLH1_data <- read.csv(file.path(path, "MLH1_data.csv"))
Tiles <- read.csv(file.path(path, "Tile_1-7_data.csv"))
Replicate_correlation <- read.csv(file.path(path, "Replicate_correlation.csv"))
rSASA_Secondary <- read.csv(file.path(path, "rSASA_Secondary.csv"))


#Fig. S4B
ggsave(file.path(path, "Fig.S4B.pdf"), ggplot(Tiles, aes(x = score_abundance, fill = Tile)) +
  geom_density(alpha = 0.3) + 
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1.7), color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("label", x = 0, y = 1.7, label = "WT", color = "black", size = 3) +
  scale_x_continuous(limits = c(-3, 2), breaks = c(-3, -2, -1, 0, 1, 2)) +
  scale_y_continuous(limits = c(0, 1.7), breaks = c(0, 0.5, 1, 1.5)) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.3))) +
  xlab("Abundance score") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 6, height = 6)


#Fig. S4C
ggsave(file.path(path, "Fig.S4C.pdf"), ggplot(Tiles, aes(x = score_interaction, fill = Tile)) +
  geom_density(alpha = 0.3) + 
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 7), color = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("label", x = 0, y = 7, label = "WT", color = "black", size = 3) +
  scale_x_continuous(limits = c(-2, 3), breaks = c(-2, -1, 0, 1, 2, 3)) +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6)) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.3))) +
  xlab("Interaction score") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 6, height = 6)


#Fig. S5
ggsave(file.path(path, "Fig.S51.pdf"), ggplot(Replicate_correlation, aes(score_abundance_rep1, score_abundance_rep2)) +
  geom_point() +
  stat_cor(method = "pearson", aes(label = ..r.label..)) +
  stat_poly_line() + 
  xlab("Abundance_rep1") + 
  ylab("Abundance_rep2") + 
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw(), width = 6, height = 6)
ggsave(file.path(path, "Fig.S52.pdf"), ggplot(Replicate_correlation, aes(score_abundance_rep1, score_abundance_rep3)) +
  geom_point() +
  stat_cor(method = "pearson", aes(label = ..r.label..)) +
  stat_poly_line() + 
  xlab("Abundance_rep1") + 
  ylab("Abundance_rep3") + 
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw(), width = 6, height = 6)
ggsave(file.path(path, "Fig.S53.pdf"), ggplot(Replicate_correlation, aes(score_abundance_rep2, score_abundance_rep3)) +
  geom_point() +
  stat_cor(method = "pearson", aes(label = ..r.label..)) +
  stat_poly_line() + 
  xlab("Abundance_rep2") + 
  ylab("Abundance_rep3") + 
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw(), width = 6, height = 6)

ggsave(file.path(path, "Fig.S54.pdf"), ggplot(Replicate_correlation, aes(score_interaction_rep1, score_interaction_rep2)) +
  geom_point() +
  stat_cor(method = "pearson", aes(label = ..r.label..)) +
  stat_poly_line() + 
  xlab("Interaction_rep1") + 
  ylab("Interaction_rep2") + 
  theme_bw(), width = 6, height = 6)
ggsave(file.path(path, "Fig.S55.pdf"), ggplot(Replicate_correlation, aes(score_interaction_rep1, score_interaction_rep3)) +
  geom_point() +
  stat_cor(method = "pearson", aes(label = ..r.label..)) +
  stat_poly_line() + 
  xlab("Interaction_rep1") + 
  ylab("Interaction_rep3") + 
  theme_bw(), width = 6, height = 6)
ggsave(file.path(path, "Fig.S56.pdf"), ggplot(Replicate_correlation, aes(score_interaction_rep2, score_interaction_rep3)) +
  geom_point() +
  stat_cor(method = "pearson", aes(label = ..r.label..)) +
  stat_poly_line() + 
  xlab("Interaction_rep2") + 
  ylab("Interaction_rep3") + 
  theme_bw(), width = 6, height = 6)


#Fig. S6A
order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Pro', 'Ala', 'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Ter')
wt <- MLH1_data %>% 
  select(position, wt_aa) %>%
  distinct(position, wt_aa) %>% 
  arrange(position)
heatmap_SE_abundance <- MLH1_data %>%
  select(position, sub_aa, SE_abundance) %>%
  spread(position, SE_abundance) %>%
  arrange(factor(sub_aa, order_aa)) %>%
  column_to_rownames("sub_aa") %>%
  as.matrix()
heatmap_SE_abundance_wt <- heatmap_SE_abundance
for (i in 1:270){
  wt_res <- wt[i,'wt_aa']
  heatmap_SE_abundance_wt[wt_res, i] <- 10
}
Fig.S6A <- draw(Heatmap(heatmap_SE_abundance_wt,                                         
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(seq(0.007654119, 0.4539574, length.out = 199), 10), c(colorRampPalette(c("white", "darkgreen"))(199), "yellow")),
             name = ("Abundance SE"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         legend_width = unit(20, "cm"),
                                         labels_gp = gpar(fontsize = 20),
                                         title_gp = gpar(fontsize = 20, fontface = "bold"),
                                         border = "black", 
                                         at = c(0, 0.15, 0.3, 0.45),
                                         labels = c("0", "0.15", "0.3", "0.45")),
             row_names_gp = gpar(fontsize = 20),
             column_names_gp = gpar(fontsize = 16),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 0.5),
             rect_gp = gpar(col = "black", lwd = 0.5),
             na_col = "dimgrey"), 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Wild-type amino acid", "Missing variant"),
                                          legend_gp = gpar(fill = c("yellow", "dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 20))))
pdf(file.path(path, "Fig.S6A.pdf"), width = 50, height = 10)
Fig.S6A <- Fig.S6A
draw(Fig.S6A)
dev.off()


#Fig. S6B
heatmap_SE_interaction <- MLH1_data %>%
  select(position, sub_aa, SE_interaction) %>%
  spread(position, SE_interaction) %>%
  arrange(factor(sub_aa, order_aa)) %>%
  column_to_rownames("sub_aa") %>%
  as.matrix()
heatmap_SE_interaction_wt <- heatmap_SE_interaction
for (i in 1:270){
  wt_res <- wt[i,'wt_aa']
  heatmap_SE_interaction_wt[wt_res, i] <- 10
}
Fig.S6B <- draw(Heatmap(heatmap_SE_interaction_wt,                                         
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(seq(0.0129505, 0.298587, length.out = 199), 10), c(colorRampPalette(c("white", "darkgreen"))(199), "yellow")),
             name = ("Interaction SE"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         legend_width = unit(20, "cm"),
                                         labels_gp = gpar(fontsize = 20),
                                         title_gp = gpar(fontsize = 20, fontface = "bold"),
                                         border = "black", 
                                         at = c(0, 0.1, 0.2),
                                         labels = c("0", "0.1", "0.2")),
             row_names_gp = gpar(fontsize = 20),
             column_names_gp = gpar(fontsize = 16),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 0.5),
             rect_gp = gpar(col = "black", lwd = 0.5),
             na_col = "dimgrey"), 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Wild-type amino acid", "Missing variant"),
                                          legend_gp = gpar(fill = c("yellow", "dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 20))))
pdf(file.path(path, "Fig.S6B.pdf"), width = 50, height = 10)
Fig.S6B <- Fig.S6B
draw(Fig.S6B)
dev.off()


#Fig. S7A
rSASA_Secondary$Secondary_structure <- factor(rSASA_Secondary$Secondary_structure , levels = c("Loops", "Helices", "Sheets"))
ggsave(file.path(path, "Fig.S7A.pdf"), ggplot() +
  geom_point(rSASA_Secondary, mapping = aes(score_abundance_median, rSASA), size = 2.5) +
  geom_point(rSASA_Secondary, mapping = aes(score_abundance_median, rSASA, color = Secondary_structure), size = 1.5) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_color_manual(values = c("#73A222", "#FFF670", "#CF47FF"), labels=c("Loops", "Helices", "Sheets")) +
  xlab("Median abundance score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6) 


#Fig. S7B
ggsave(file.path(path, "Fig.S7B.pdf"), ggplot() +
  geom_point(rSASA_Secondary, mapping = aes(score_interaction_median, rSASA), size = 2.5) +
  geom_point(rSASA_Secondary, mapping = aes(score_interaction_median, rSASA, color = Secondary_structure), size = 1.5) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_color_manual(values = c("#73A222", "#FFF670", "#CF47FF"), labels=c("Loops", "Helices", "Sheets")) +
  xlab("Median interaction score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S8A
LOI_GOI_variants <- MLH1_data %>%
  filter(variant == "p.Leu547Asp" |
         variant == "p.Leu547Arg" |
         variant == "p.Leu540Glu" |
         variant == "p.Leu540Arg" |
         variant == "p.Gln544Trp" |
         variant == "p.Gly602Leu" |
         variant == "p.Lys604Cys")
ggsave(file.path(path, "Fig.S8A.pdf"), ggplot() + 
  geom_hex(MLH1_data, mapping = aes(score_abundance, score_interaction), bins = 70) +
  geom_hex(LOI_GOI_variants, mapping = aes(score_abundance, score_interaction), bins = 70, fill = "red") +
  scale_fill_gradient(low = "black", high = "magenta", breaks = c(20, 40, 60)) +
  scale_x_continuous(limits = c(-1.27,1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3)) +
  guides(fill = guide_colourbar(title = "Count")) +
  xlab("Abundance score") +
  ylab("Interaction score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S12A
LOS_WTI_variants <- MLH1_data %>%
  filter(variant == "p.Val531Asp" |
         variant == "p.Trp538Ala" |
         variant == "p.Leu540Ser" |
         variant == "p.Thr553Phe" |
         variant == "p.Gln562Tyr")
ggsave(file.path(path, "Fig.S12A.pdf"), ggplot() + 
  geom_hex(MLH1_data, mapping = aes(score_abundance, score_interaction), bins = 70) +
  geom_hex(LOS_WTI_variants, mapping = aes(score_abundance, score_interaction), bins = 70, fill = "red") +
  scale_fill_gradient(low = "black", high = "magenta", breaks = c(20, 40, 60)) +
  scale_x_continuous(limits = c(-1.27,1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3)) +
  guides(fill = guide_colourbar(title = "Count")) +
  xlab("Abundance score") +
  ylab("Interaction score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)

