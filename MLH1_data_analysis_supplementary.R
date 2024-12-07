#Load in required packages
library(tidyverse)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)
library(ggExtra)
library(ggpubr)
library(pROC)
library(ggpmisc)
library(hexbin)

#Change the path to your directory
path <- "your_path"

#Read in dataframes
MLH1_data <- read.csv(file.path(path, "MLH1_data.csv"))
Tiles <- read.csv(file.path(path, "Tile_1-7_data.csv"))
Replicate_correlation <- read.csv(file.path(path, "Replicate_correlation.csv"))
rSASA_Secondary <- read.csv(file.path(path, "rSASA_Secondary.csv"))
Rosetta_Gemme_data <- read.csv(file.path(path, "Rosetta_Gemme_data.csv"))
ClinVar_gnomAD_data <- read.csv(file.path(path, "ClinVar_gnomAD_data.csv"))
Benchmark_Hinrichsen <- read_csv(file.path(path, "Benchmark_Hinrichsen.csv"))
Benchmark_Abildgaard <- read_csv(file.path(path, "Benchmark_Abildgaard.csv"))
Benchmark_Kosinski <- read_csv(file.path(path, "Benchmark_Kosinski.csv"))


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
ggsave(file.path(path, "Fig.S7A.pdf"), ggplot(rSASA_Secondary) +
  geom_point(mapping = aes(score_abundance_median, rSASA), size = 2.5) +
  geom_point(mapping = aes(score_abundance_median, rSASA, color = Secondary_structure), size = 1.5) +
  stat_cor(method = "spearman", mapping = aes(score_abundance_median, rSASA, label = after_stat(r.label))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_color_manual(values = c("#73A222", "#FFF670", "#CF47FF"), labels=c("Loops", "Helices", "Sheets")) +
  xlab("Median abundance score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6) 


#Fig. S7B
ggsave(file.path(path, "Fig.S7B.pdf"), ggplot(rSASA_Secondary) +
  geom_point(mapping = aes(score_interaction_median, rSASA), size = 2.5) +
  geom_point(mapping = aes(score_interaction_median, rSASA, color = Secondary_structure), size = 1.5) +
  stat_cor(method = "spearman", mapping = aes(score_interaction_median, rSASA, label = after_stat(r.label))) +
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
  scale_x_continuous(limits = c(-1.27, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3)) +
  guides(fill = guide_colourbar(title = "Count")) +
  xlab("Abundance score") +
  ylab("Interaction score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S11A
ggsave(file.path(path, "Fig.S11A.pdf"), ggplot(Rosetta_Gemme_data) + 
  geom_hex(mapping = aes(score_abundance, Rosetta_ddG), bins = 70) +
  stat_cor(method = "spearman", mapping = aes(score_abundance, Rosetta_ddG, label = after_stat(r.label))) +
  scale_fill_gradient(low = "black", high = "magenta", breaks = c(5, 15, 25)) +
  scale_x_continuous(limits = c(-1.27, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  guides(fill = guide_colourbar(title = "Count")) +
  xlab("Abundance score") +
  ylab(expression(Rosetta~Delta*Delta*G)) +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S11B
ggsave(file.path(path, "Fig.S11B.pdf"), ggplot(Rosetta_Gemme_data) + 
  geom_hex(mapping = aes(score_abundance, Gemme_ddE), bins = 70) +
  stat_cor(method = "spearman", mapping = aes(score_abundance, Gemme_ddE, label = after_stat(r.label))) +
  scale_fill_gradient(low = "black", high = "magenta", breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(-1.27, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  guides(fill = guide_colourbar(title = "Count")) +
  xlab("Abundance score") +
  ylab(expression(GEMME~Delta*Delta*E)) +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S11C
ggsave(file.path(path, "Fig.S11C.pdf"), ggplot(Rosetta_Gemme_data) + 
  geom_hex(mapping = aes(score_interaction, Rosetta_ddG), bins = 70) +
  stat_cor(method = "spearman", mapping = aes(score_interaction, Rosetta_ddG, label = after_stat(r.label))) +
  scale_fill_gradient(low = "black", high = "magenta", breaks = c(20, 40, 60)) +
  guides(fill = guide_colourbar(title = "Count")) +
  xlab("Interaction score") +
  ylab(expression(Rosetta~Delta*Delta*G)) +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S11D
ggsave(file.path(path, "Fig.11D.pdf"), ggplot(Rosetta_Gemme_data) + 
  geom_hex(mapping = aes(score_interaction, Gemme_ddE), bins = 70) +
  stat_cor(method = "spearman", mapping = aes(score_interaction, Gemme_ddE, label = after_stat(r.label))) +
  scale_fill_gradient(low = "black", high = "magenta", breaks = c(10, 20, 30)) +
  guides(fill = guide_colourbar(title = "Count")) +
  xlab("Interaction score") +
  ylab(expression(GEMME~Delta*Delta*E)) +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S12
ROC <- Rosetta_Gemme_data %>%
  mutate(cutoff_abundance = ifelse(score_abundance > -0.5, 1, 0),
         cutoff_interaction = ifelse(score_interaction > -0.5, 1, 0))
ggsave(file.path(path, "Fig.S12.pdf"), ggroc(list(rosetta_abundance = roc(ROC$cutoff_abundance, ROC$Rosetta_ddG),
  gemme_abundance = roc(ROC$cutoff_abundance, ROC$Gemme_ddE),
  rosetta_interaction = roc(ROC$cutoff_interaction, ROC$Rosetta_ddG),
  gemme_interaction = roc(ROC$cutoff_interaction, ROC$Gemme_ddE)),
  legacy.axes = TRUE, linewidth = 2, aes = c("color", "linetype")) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", labels = c(bquote(Rosetta~Delta*Delta*G: "Abundance score " (AUC == .(round(auc(ROC$cutoff_abundance, ROC$Rosetta_ddG), 2)))),
                                       bquote(GEMME~Delta*Delta*E: "Abundance score " (AUC == .(round(auc(ROC$cutoff_abundance, ROC$Gemme_ddE), 2)))),
                                       bquote(Rosetta~Delta*Delta*G: "Interaction score " (AUC == .(round(auc(ROC$cutoff_interaction, ROC$Rosetta_ddG), 2)))),
                                       bquote(GEMME~Delta*Delta*E: "Interaction score " (AUC == .(round(auc(ROC$cutoff_interaction, ROC$Gemme_ddE), 2))))),
                            values = c(1, 1, 1, 1)) + 
  scale_color_manual("", labels = c(bquote(Rosetta~Delta*Delta*G: "Abundance score " (AUC == .(round(auc(ROC$cutoff_abundance, ROC$Rosetta_ddG), 2)))),
                                    bquote(GEMME~Delta*Delta*E: "Abundance score " (AUC == .(round(auc(ROC$cutoff_abundance, ROC$Gemme_ddE), 2)))),
                                    bquote(Rosetta~Delta*Delta*G: "Interaction score " (AUC == .(round(auc(ROC$cutoff_interaction, ROC$Rosetta_ddG), 2)))),
                                    bquote(GEMME~Delta*Delta*E: "Interaction score " (AUC == .(round(auc(ROC$cutoff_interaction, ROC$Gemme_ddE), 2))))),
                            values = c("orange", "magenta3", "limegreen", "dodgerblue2")) +
  geom_abline(linetype = "dashed", linewidth = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0, aspect.ratio = 1), width = 9, height = 8)


#Fig. S13
ggsave(file.path(path, "Fig.S13.pdf"), ggMarginal(ggplot(Rosetta_Gemme_data, aes(score_abundance, Gemme_ddE)) +
  geom_point(aes(color = score_interaction), size = 3, alpha = ifelse(Rosetta_Gemme_data$score_interaction > -0.2 & Rosetta_Gemme_data$score_interaction < 0.2, 0.2, 1)) + 
  scale_colour_gradient2(limits = c(-1.902492, 1.902492), oob = scales::squish, breaks = c(-1.6, 0, 1.6), labels = c("Weak interaction", "WT interaction", "Strong interaction")) +
  labs(color = "Interaction score") +
  scale_x_continuous(limits = c(-1.27, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  xlab("Abundance score") +  
  ylab(expression(GEMME~Delta*Delta*E)) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "grey80", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank()), type = "density"), width = 9, height = 8)


#Fig. S14
Clin_ROC <- ClinVar_gnomAD_data %>%
  filter(ClinVar != "Uncertain significance") %>%
  mutate(cutoff_clin = ifelse(ClinVar == "Benign" | ClinVar == "Likely benign", 1, 0)) %>%
  left_join(Rosetta_Gemme_data %>% select(variant, Rosetta_ddG, Gemme_ddE), by = "variant")
ggsave(file.path(path, "Fig.S14.pdf"), ggroc(list(abundance_clin = roc(Clin_ROC$cutoff_clin, Clin_ROC$score_abundance),
  interaction_clin = roc(Clin_ROC$cutoff_clin, Clin_ROC$score_interaction),
  rosetta_clin = roc(Clin_ROC$cutoff_clin, Clin_ROC$Rosetta_ddG),
  gemme_clin = roc(Clin_ROC$cutoff_clin, Clin_ROC$Gemme_ddE)),
  legacy.axes = TRUE, linewidth = 2, aes = c("color", "linetype")) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", labels = c(bquote(Abundance : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$score_abundance), 2)))),
                                       bquote(Interaction : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$score_interaction), 2)))),
                                       bquote(Rosetta~Delta*Delta*G : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$Rosetta_ddG), 2)))),
                                       bquote(GEMME~Delta*Delta*E : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$Gemme_ddE), 2))))),
                            values = c(1, 1, 1, 1)) + 
         scale_color_manual("", labels = c(bquote(Abundance : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$score_abundance), 2)))),
                                           bquote(Interaction : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$score_interaction), 2)))),
                                           bquote(Rosetta~Delta*Delta*G : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$Rosetta_ddG), 2)))),
                                           bquote(GEMME~Delta*Delta*E : "ClinVar " (AUC == .(round(auc(Clin_ROC$cutoff_clin, Clin_ROC$Gemme_ddE), 2))))),
                                values = c("orange", "magenta3", "limegreen", "dodgerblue2")) +
  geom_abline(linetype = "dashed", linewidth = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0, aspect.ratio = 1), width = 9, height = 8)


#Fig. S16A
LOS_WTI_variants <- MLH1_data %>%
  filter(variant == "p.Val531Asp" |
         variant == "p.Trp538Ala" |
         variant == "p.Leu540Ser" |
         variant == "p.Thr553Phe" |
         variant == "p.Gln562Tyr")
ggsave(file.path(path, "Fig.S16A.pdf"), ggplot() + 
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


#Fig. S17A
Benchmark_Hinrichsen_sign <- Benchmark_Hinrichsen %>%
  filter(pval < 0.05)
ggsave(file.path(path, "Fig.S17A.pdf"), ggplot(Benchmark_Hinrichsen, mapping = aes(score_abundance, pval)) +
  geom_point() +
  geom_point(Benchmark_Hinrichsen_sign, mapping = aes(score_abundance, pval), color = "red") +
  stat_cor(method = "spearman", mapping = aes(score_abundance, pval, label = after_stat(r.label))) +
  stat_poly_line() + 
  scale_x_continuous(limits = c(-1.27, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  xlab("Abundance score") + 
  ylab("P-values for expression relative to WT") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S17B
ggsave(file.path(path, "Fig.S17B.pdf"), ggplot(Benchmark_Abildgaard, mapping = aes(score_abundance, percentage)) +
  geom_point() +
  geom_errorbarh(aes(xmin = score_abundance-SD_score_abundance, xmax = score_abundance+SD_score_abundance), alpha = 0.4) +
  geom_errorbar(aes(ymin = percentage-SD_percentage, ymax = percentage+SD_percentage), alpha = 0.4) +
  scale_x_continuous(limits = c(-1.27, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 145), breaks = c(0, 25, 50, 75, 100)) +
  stat_cor(method = "spearman", mapping = aes(score_abundance, percentage, label = after_stat(r.label))) +
  stat_poly_line() + 
  xlab("Abundance score") + 
  ylab("Steady-state levels (% of WT)") +  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), width = 7, height = 6)


#Fig. S17C
ggsave(file.path(path, "Fig.S17C.pdf"), ggplot(Benchmark_Kosinski, mapping = aes(x = reorder(variant, -score_interaction), y = score_interaction, fill = int)) + 
  geom_bar(stat = "identity") + 
  guides(fill = guide_legend("Interaction")) +
  scale_fill_manual(breaks = c("Yes", "No", NA), values = c("darkgreen", "darkred", "darkgrey")) +
  xlab("Variant") + 
  ylab("Interaction score") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90)), width = 7, height = 6)
