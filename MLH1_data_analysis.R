#Load in required packages
library(tidyverse)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)
library(ggExtra)

#Change the path to your directory
path <- ""

#Read in dataframes
MLH1_data <- read.csv(file.path(path, "MLH1_data.csv"))
Rosetta_Gemme_data <- read.csv(file.path(path, "Rosetta_Gemme_data.csv"))
ClinVar_gnomAD_data <- read.csv(file.path(path, "ClinVar_gnomAD_data.csv"))


#Fig. 2C (left panel)
L550P_abundance <- MLH1_data$score_abundance[MLH1_data$variant == "p.Leu550Pro"]
ggplot(MLH1_data, aes(score_abundance)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "lightgreen") +
  geom_segment(mapping = aes(x = 0, y = 0, xend = 0, yend = 500), color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_segment(mapping = aes(x = L550P_abundance, y = 0, xend = L550P_abundance, yend = 500), color = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("label", x = 0, y = 500, label = "WT" , color = "black", size = 3) +
  annotate("label", x = L550P_abundance, y = 500, label = "L550P" , color = "red", size = 3) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 500), breaks = c(0, 100, 200, 300, 400, 500)) +
  xlab("Abundance score") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#Fig. 2C (right panel)
L550P_interaction <- MLH1_data$score_interaction[MLH1_data$variant == "p.Leu550Pro"]
ggplot(MLH1_data, aes(score_interaction)) + 
  geom_histogram(binwidth = 0.115, color = "black", fill = "lightgreen") +
  geom_segment(mapping = aes(x = 0, y = 0, xend = 0, yend = 1700), color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_segment(mapping = aes(x = L550P_interaction, y = 0, xend = L550P_interaction, yend = 1700), color = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("label", x = 0, y = 1700, label = "WT" , color = "black", size = 3) +
  annotate("label", x = L550P_interaction, y = 1700, label = "L550P" , color = "red", size = 3) +
  scale_x_continuous(limits = c(-2, 3), breaks = c(-2, -1, 0, 1, 2, 3)) +
  scale_y_continuous(limits = c(0, 1700), breaks = c(0, 500, 1000, 1500)) +
  xlab("3-AT interaction score") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#Fig. 3A
order_aa <- c('His', 'Lys', 'Arg', 'Asp', 'Glu', 'Cys', 'Met', 'Asn', 'Gln', 'Ser', 'Thr', 'Pro', 'Ala', 'Ile', 'Leu', 'Val', 'Phe', 'Trp', 'Tyr', 'Gly', 'Ter')
wt <- MLH1_data %>% 
  select(position, wt_aa) %>%
  distinct(position, wt_aa) %>% 
  arrange(position)
heatmap_abundance <- MLH1_data %>%
  select(position, sub_aa, score_abundance) %>%
  spread(position, score_abundance) %>%
  arrange(factor(sub_aa, order_aa)) %>%
  column_to_rownames("sub_aa") %>%
  as.matrix()
heatmap_abundance_wt <- heatmap_abundance
for (i in 1:270){
  wt_res <- wt[i,'wt_aa']
  heatmap_abundance_wt[wt_res, i] <- 10
}
draw(Heatmap(heatmap_abundance_wt,                                         
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(seq(-1.240268, 1.240268, length.out = 199), 10), c(colorRampPalette(c("darkred", "white", "darkblue"))(199), "yellow")),
             name = ("Abundance score"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         legend_width = unit(10, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-1, 0, 1),
                                         labels = c("Low\nabundance", "WT\nabundance", "High\nabundance")),
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 7),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 0.5),
             rect_gp = gpar(col = "black", lwd = 0.5),
             na_col = "dimgrey"), 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Wild-type amino acid", "Missing variant"),
                                          legend_gp = gpar(fill = c("yellow", "dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 12))))


#Fig. 3C
heatmap_interaction <- MLH1_data %>%
  select(position, sub_aa, score_interaction) %>%
  spread(position, score_interaction) %>%
  arrange(factor(sub_aa, order_aa)) %>%
  column_to_rownames("sub_aa") %>%
  as.matrix()
heatmap_interaction_wt <- heatmap_interaction
for (i in 1:270){
  wt_res <- wt[i,'wt_aa']
  heatmap_interaction_wt[wt_res, i] <- 10
}
draw(Heatmap(heatmap_interaction_wt,                                         
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(seq(-1.902492, 1.902492, length.out = 199), 10), c(colorRampPalette(c("darkred", "white", "darkblue"))(199), "yellow")),
             name = ("Interaction score"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         legend_width = unit(10, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-2, 0, 2),
                                         labels = c("Weak\ninteraction", "WT\ninteraction", "Strong\ninteraction")),
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 7),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 0.5),
             rect_gp = gpar(col = "black", lwd = 0.5),
             na_col = "dimgrey"), 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Wild-type amino acid", "Missing variant"),
                                          legend_gp = gpar(fill = c("yellow", "dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 12))))


#Fig. 4A
ggMarginal(ggplot() + 
             geom_point(Rosetta_Gemme_data, mapping = aes(Rosetta_ddG, Gemme_ddE), color = "black", size = 1) + 
             geom_point(Rosetta_Gemme_data, mapping = aes(Rosetta_ddG, Gemme_ddE, color = score_abundance), size = 0.5) + 
             labs(color = "Abundance score") + 
             scale_colour_gradient2(limits=c(-1.240268, 1.240268), oob = scales::squish, breaks = c(-1, 0, 1), labels = c("Low abundance", "WT abundance", "High abundance")) + 
             xlab(expression(Rosetta~Delta*Delta*G)) +  
             ylab(expression(GEMME~Delta*Delta*E)) + 
             theme_bw() + 
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), type = "density")


#Fig. 4B
ggMarginal(ggplot() + 
             geom_point(Rosetta_Gemme_data, mapping = aes(Rosetta_ddG, Gemme_ddE), color = "black", size = 1) + 
             geom_point(Rosetta_Gemme_data, mapping = aes(Rosetta_ddG, Gemme_ddE, color = score_interaction), size = 0.5) + 
             labs(color = "Interaction score") + 
             scale_colour_gradient2(limits=c(-1.902492, 1.902492), oob = scales::squish, breaks = c(-1.6, 0, 1.6), labels = c("Weak interaction", "WT interaction", "Strong interaction")) +
             xlab(expression(Rosetta~Delta*Delta*G)) +  
             ylab(expression(GEMME~Delta*Delta*E)) + 
             theme_bw() + 
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()), type = "density")


#Fig. 4C
VUS_variants <- ClinVar_gnomAD_data %>%
  filter(ClinVar == "Uncertain significance")
ClinVar_variants <- ClinVar_gnomAD_data %>%
  filter(ClinVar == "Benign" |
         ClinVar == "Pathogenic" |
         ClinVar == "Likely benign" |
         ClinVar == "Likely pathogenic")
ggplot() +
  geom_point(VUS_variants, mapping = aes(x = log10(gnomAD), y = score_abundance), color = "gray90", size = 2) +
  geom_point(ClinVar_variants, mapping = aes(x = log10(gnomAD), y = score_abundance), color = "black", size = 2) +
  geom_point(ClinVar_variants, mapping = aes(x = log10(gnomAD), y = score_abundance, color = ClinVar), size = 1.5) +
  scale_color_manual("ClinVar", values = c("darkgreen", "green", "yellow", "red")) +
  scale_x_continuous(labels = function(x) parse(text = sprintf("10^%d", x)), breaks = c(-6, -5, -4, -3)) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 1)) +
  xlab("gnomAD allele frequency (log)") +
  ylab("Abundance score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#Fig. 4D
ggplot() +
  geom_point(VUS_variants, mapping = aes(x = log10(gnomAD), y = score_interaction), color = "gray90", size = 2) +
  geom_point(ClinVar_variants, mapping = aes(x = log10(gnomAD), y = score_interaction), color = "black", size = 2) +
  geom_point(ClinVar_variants, mapping = aes(x = log10(gnomAD), y = score_interaction, color = ClinVar), size = 1.5) +
  scale_color_manual("ClinVar", values = c("darkgreen", "green", "yellow", "red")) +
  scale_x_continuous(labels = function(x) parse(text = sprintf("10^%d", x)), breaks = c(-6, -5, -4, -3)) +
  xlab("gnomAD allele frequency (log)") +
  ylab("Interaction score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
