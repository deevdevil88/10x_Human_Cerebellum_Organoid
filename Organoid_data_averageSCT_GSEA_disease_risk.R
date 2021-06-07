library(Seurat)
library(tidyverse)
library(pheatmap)
library(cowplot)
oragnoid <- readRDS(file="./20200918vs2_v3update_newganoidscombinednames.rds")
oragnoid[["cell_type"]] <- Idents(oragnoid)
Idents(oragnoid) <- gsub(pattern = " ", replacement = "_", x = Idents(oragnoid))
Avg_organoid <- AverageExpression(oragnoid, return.seurat = T)

Avg_sct_log <- GetAssayData(Avg_organoid,slot = "data", assay = "SCT")
Avg_sct_log <- as.data.frame(Avg_sct_log)
Avg_sct_log <- Avg_sct_log %>% rownames_to_column(var = "NAME")

Avg_sct_log$DESCRIPTION <- "NA"
Avg_sct_log_2 <- Avg_sct_log[, c(1,14,2:13)]
write.table(Avg_sct_log_2, file = "Organoid_Avg_log_SCT_12PCW.txt", sep="\t",quote = F,row.names = F)





#disease_risk_genes <- read_csv(file="./AldingerRiskGenes.csv")

disease_risk_genes <- read_csv(file="./SupplementaryTable33_RiskGenes.csv")

write.table(disease_risk_genes, file="disease_genes_v2.gmx", sep="\t", quote = F,row.names = F)


disease_risk_genes_v2 <- read_csv(file="./ALTERNATE_SupplementaryTable33_RiskGenes.csv")

write.table(disease_risk_genes_v2, file="disease_genes_v2_alternate.gmx", sep="\t", quote = F,row.names = F)


# Main Gene sets 

GSEA_results <- read.table(file="./GSEA_27oct_main/GSEA_disease_risk_organoid_celltypes_v2.txt", sep="\t",header=T)
GSEA_results <- GSEA_results[,1:9]
GSEA_results$p_adj <- p.adjust(GSEA_results$NOM.p.val, method = "BH")
GSEA_results$bon <- p.adjust(GSEA_results$NOM.p.val, method="bonferroni")
GSEA_results$log10_padj <- -log10(GSEA_results$p_adj)
GSEA_results$bon_log10_padj <- -log10(GSEA_results$bon)
GSEA_results

GSEA_results$direction <- ifelse(GSEA_results$NES > 0,"Up","Down")
#GSEA_results$Celltype <- gsub("_"," ",GSEA_results$Celltype)

write.table(GSEA_results, file="GSEA_disease_risk_organoid_celltypes_final_v2.txt", sep="\t",quote=F)




 DefaultAssay(Avg_organoid) <- "SCT"
Autism_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$AUTISM, slot = "scale.data")
Autism_heatmap_data_t <- t(Autism_heatmap_data)

Autism_heatmap <- pheatmap(Autism_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 5, fontsize = 3, clustering_method = "ward.D2",border_color = NA)

Autism_results <- GSEA_results %>% filter(Disease == "AUTISM") 

level_order <- c('Glutamatergic DCN', 'Granule cells-2', 'Endothelial','Granule cells-1','GC precusors','Rhombic lip','Ciliated cells','Unknown','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
Autism_bar_plot <-  ggplot(Autism_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
Autism_bar_plot <- Autism_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
Autism_bar_plot <- Autism_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
Autism_bar_plot <- Autism_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
Autism_bar_plot <- Autism_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -1.0, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
Autism_bar_plot

bottom_row <- plot_grid(Autism_heatmap[[4]], nrow = 1, ncol = 1, scale = 0.95)
top_row <- plot_grid(Autism_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
Autism_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="tl", rel_heights = c(0.05,0.9))
Autism_final

# now add the title
#ast <- ggdraw() +draw_text(text = " *",x =0.1,y=-0.08) + theme(plot.margin = margin(0,0,0.01,0.01))
title <- ggdraw() + 
  draw_label(
    "Autism spectrum disorder",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
Autism_final_v2 <- plot_grid(
  title,Autism_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
Autism_final_v2
ggsave(filename = "Autism_celltype_trend_bonf_v2.eps",plot = Autism_final_v2,width = 10.5, height = 18.5, units = "cm")

ggsave(filename = "Autism_celltype_trend_bonf_v2.pdf",plot = Autism_final_v2,width = 10.5, height = 18.5, units = "cm")


# CBLM 
CBLM_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$CBLM, slot = "scale.data")
CBLM_heatmap_data_t <- t(CBLM_heatmap_data)

CBLM_heatmap <- pheatmap(CBLM_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 5, fontsize = 4, clustering_method = "ward.D2",legend = T,border_color = NA)

CBLM_results <- GSEA_results %>% filter(Disease == "CBLM") 

level_order <- c('Endothelial','Rhombic lip','Granule cells-1','Granule cells-2','Glutamatergic DCN','GC precusors','Ciliated cells','Unknown','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
CBLM_bar_plot <-  ggplot(CBLM_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
CBLM_bar_plot <- CBLM_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
CBLM_bar_plot <- CBLM_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
CBLM_bar_plot <- CBLM_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
CBLM_bar_plot <- CBLM_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
CBLM_bar_plot

bottom_row <- plot_grid(CBLM_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.945)
top_row <- plot_grid(CBLM_bar_plot, nrow = 1, ncol = 1, scale = 0.985)
CBLM_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
CBLM_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.34,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Cerebellar malformations",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
CBLM_final_v2 <- plot_grid(
  title,CBLM_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
CBLM_final_v2
ggsave(filename = "CBLM_celltype_trend_bonf_v2.eps",plot = CBLM_final_v2,width = 8.0, height = 11, units = "cm")

ggsave(filename = "CBLM_celltype_trend_bonf_v2.pdf",plot = CBLM_final_v2,width = 8.0, height = 11, units = "cm")





# SCA/EA
SCA_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$`SCA/EA`, slot = "scale.data")
SCA_heatmap_data_t <- t(SCA_heatmap_data)

SCA_heatmap <- pheatmap(SCA_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 4, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

SCA_results <- GSEA_results %>% filter(Disease == "SCA/EA") 

level_order <- c('Ciliated cells','Choroid plexus','Roof plate','Granule cells-1','Granule cells-2','Bergmann glia','GC precusors','Rhombic lip','Glutamatergic DCN','Endothelial','Purkinje neurons','Unknown')
SCA_bar_plot <-  ggplot(SCA_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
SCA_bar_plot <- SCA_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
SCA_bar_plot <- SCA_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
SCA_bar_plot <- SCA_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
SCA_bar_plot <- SCA_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
SCA_bar_plot

bottom_row <- plot_grid(SCA_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.92)
top_row <- plot_grid(SCA_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
SCA_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
SCA_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "SCA/EA",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
SCA_final_v2 <- plot_grid(
  title,SCA_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
SCA_final_v2
ggsave(filename = "SCA_celltype_trend_bonf_v2.eps",plot = SCA_final_v2,width = 6.5, height = 10, units = "cm")

ggsave(filename = "SCA_celltype_trend_bonf_v2.pdf",plot = SCA_final_v2,width = 6.5, height = 10, units = "cm")


muultipanel_A <- plot_grid(CBLM_final_v2,SCA_final_v2, align = "v", nrow = 2, rel_heights = c(0.7,0.45) )
muultipanel_A

panel_B <- plot_grid(muultipanel_A,Autism_final_v2, align="h", ncol = 2, rel_widths = c(0.7,0.7))
panel_B

ggsave(filename = "disease_risk_main.eps",plot = panel_B, width = 18.5, height = 20, units = "cm")

ggsave(filename = "disease_risk_main.pdf",plot = panel_B, width = 18.5, height = 20, units = "cm")


# ID
ID_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$ID, slot = "scale.data")
ID_heatmap_data_t <- t(ID_heatmap_data)

ID_heatmap <- pheatmap(ID_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 2.0, cutree_rows = 5, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

ID_results <- GSEA_results %>% filter(Disease == "ID") 

level_order <- c('Granule cells-1','Granule cells-2','GC precusors','Glutamatergic DCN','Rhombic lip','Unknown','Endothelial','Ciliated cells','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
ID_bar_plot <-  ggplot(ID_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
ID_bar_plot <- ID_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
ID_bar_plot <- ID_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
ID_bar_plot <- ID_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
ID_bar_plot <- ID_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0, 0, -1.0, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
ID_bar_plot

bottom_row <- plot_grid(ID_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.945)
top_row <- plot_grid(ID_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
ID_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="tl", rel_heights = c(0.05,0.9))
ID_final

# now add the title
#ast <- ggdraw() +draw_text(text = "* *         *",x =0.35,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Intellectual disability",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
ID_final_v2 <- plot_grid(
  title,ID_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
ID_final_v2

ggsave(filename = "ID_celltype_trend_bonf_v2.eps",plot = ID_final_v2,width = 10.5, height = 19.5, units = "cm")

ggsave(filename = "ID_celltype_trend_bonf_v2.pdf",plot = ID_final_v2,width = 10.5, height = 19.5, units = "cm")





muultipanel_A <- plot_grid(CBLM_final_v2,SCA_final_v2, align = "v", nrow = 2, rel_heights = c(0.6,0.45) )
muultipanel_A

panel_B <- plot_grid(Autism_final_v2, ID_final_v2, align="h", ncol = 2, rel_widths = c(0.6,0.6), axis = "tlrb", scale = c(0.99,0.99))
panel_B


panel_C <- plot_grid(muultipanel_A,panel_B,align = "h", ncol = 2,rel_widths = c(0.35,0.65), axis = "tlrb", scale = c(0.975,1))
panel_C
ggsave(filename = "disease_risk_main.eps",plot = panel_C, width = 19, height = 20.5, units = "cm")

ggsave(filename = "disease_risk_main.pdf",plot = panel_C, width = 19, height = 20.5, units = "cm")


# JOUBERT
joubert_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$JOUBERT, slot = "scale.data")
joubert_heatmap_data_t <- t(joubert_heatmap_data)

joubert_heatmap <- pheatmap(joubert_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 3, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

joubert_results <- GSEA_results %>% filter(Disease == "JOUBERT") 

level_order <- c('Ciliated cells','Granule cells-1','Granule cells-2','Purkinje neurons','Unknown','Endothelial','Choroid plexus','Roof plate','Glutamatergic DCN','Rhombic lip','Bergmann glia','GC precusors')
joubert_bar_plot <-  ggplot(joubert_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
joubert_bar_plot <- joubert_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
joubert_bar_plot <- joubert_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
joubert_bar_plot <- joubert_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
joubert_bar_plot <- joubert_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
joubert_bar_plot

bottom_row <- plot_grid(joubert_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.94)
top_row <- plot_grid(joubert_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
joubert_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
joubert_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Joubert syndrome",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
joubert_final_v2 <- plot_grid(
  title,joubert_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
joubert_final_v2

ggsave(filename = "joubert_celltype_trend_bonf_v2.eps",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")

ggsave(filename = "joubert_celltype_trend_bonf_v2.pdf",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")



# ALZ
ALZ_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$ALZNEG, slot = "scale.data")
ALZ_heatmap_data_t <- t(ALZ_heatmap_data)

ALZ_heatmap <- pheatmap(ALZ_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 4, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

ALZ_results <- GSEA_results %>% filter(Disease == "ALZNEG") 

level_order <- c('Unknown','Ciliated cells','Purkinje neurons','Roof plate','Choroid plexus','Endothelial','Granule cells-1','Granule cells-2','Glutamatergic DCN','Bergmann glia','GC precusors','Rhombic lip')
ALZ_bar_plot <-  ggplot(ALZ_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
ALZ_bar_plot <- ALZ_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
ALZ_bar_plot <- ALZ_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
ALZ_bar_plot <- ALZ_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
ALZ_bar_plot <- ALZ_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
ALZ_bar_plot

bottom_row <- plot_grid(ALZ_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.92)
top_row <- plot_grid(ALZ_bar_plot, nrow = 1, ncol = 1, scale = 0.98)
ALZ_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
ALZ_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Alzheimer's disease",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
ALZ_final_v2 <- plot_grid(
  title,ALZ_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
ALZ_final_v2

ggsave(filename = "ALZ_celltype_trend_bonf_v2.eps",plot = ALZ_final_v2,width = 6.5, height = 11.5, units = "cm")

ggsave(filename = "ALZ_celltype_trend_bonf_v2.pdf",plot = ALZ_final_v2,width = 6.5, height = 11.5, units = "cm")


# ARA
ARA_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$ARA, slot = "scale.data")
ARA_heatmap_data_t <- t(ARA_heatmap_data)

ARA_heatmap <- pheatmap(ARA_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 4, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

ARA_results <- GSEA_results %>% filter(Disease == "ARA") 

level_order <- c('Ciliated cells','Choroid plexus','Endothelial','Roof plate','Unknown','Glutamatergic DCN','GC precusors','Rhombic lip','Granule cells-1','Granule cells-2','Bergmann glia','Purkinje neurons')
ARA_bar_plot <-  ggplot(ARA_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
ARA_bar_plot <- ARA_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
ARA_bar_plot <- ARA_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
ARA_bar_plot <- ARA_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
ARA_bar_plot <- ARA_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
ARA_bar_plot

bottom_row <- plot_grid(ARA_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.92)
top_row <- plot_grid(ARA_bar_plot, nrow = 1, ncol = 1, scale = 0.98)
ARA_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
ARA_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "ARA",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
ARA_final_v2 <- plot_grid(
  title,ARA_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
ARA_final_v2

ggsave(filename = "ARA_celltype_trend_bonf_v2.eps",plot = ARA_final_v2,width = 6.5, height = 11.5, units = "cm")

ggsave(filename = "ARA_celltype_trend_bonf_v2.pdf",plot = ARA_final_v2,width = 6.5, height = 11.5, units = "cm")








panel_D <- plot_grid(joubert_final_v2,ALZ_final_v2,  align = "h", ncol = 2,rel_widths =c(0.5,0.5), axis = "tlrb")
panel_D


ggsave(filename = "supp_disease_risk_trend.eps",plot = panel_D,width = 14, height = 13, units = "cm")

ggsave(filename = "supp_disease_risk_trend.pdf",plot = panel_D,width = 14, height = 13, units = "cm")



# Alternative Gene set figures

GSEA_results <- read.table(file="./GSEA_27th_alternative/GSEA_disease_risk_organoid_celltypes_v2_alternative.txt", sep="\t",header=T)
GSEA_results <- GSEA_results[,1:9]
GSEA_results$p_adj <- p.adjust(GSEA_results$NOM.p.val, method = "BH")
GSEA_results$bon <- p.adjust(GSEA_results$NOM.p.val, method="bonferroni")
GSEA_results$log10_padj <- -log10(GSEA_results$p_adj)
GSEA_results$bon_log10_padj <- -log10(GSEA_results$bon)
GSEA_results

GSEA_results$direction <- ifelse(GSEA_results$NES > 0,"Up","Down")
#GSEA_results$Celltype <- gsub("_"," ",GSEA_results$Celltype)

write.table(GSEA_results, file="GSEA_disease_risk_organoid_celltypes_v2_alternative_final.txt", sep="\t",quote=F)
DefaultAssay(Avg_organoid) <- "SCT"
Autism_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes_v2$AUTISM, slot = "scale.data")
Autism_heatmap_data_t <- t(Autism_heatmap_data)

Autism_heatmap <- pheatmap(Autism_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 5, fontsize = 3, clustering_method = "ward.D2",border_color = NA)

Autism_results <- GSEA_results %>% filter(Disease == "AUTISM") 

level_order <- c('Glutamatergic DCN', 'Granule cells-2', 'Endothelial','Granule cells-1','GC precusors','Rhombic lip','Ciliated cells','Unknown','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
Autism_bar_plot <-  ggplot(Autism_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
Autism_bar_plot <- Autism_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
Autism_bar_plot <- Autism_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
Autism_bar_plot <- Autism_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
Autism_bar_plot <- Autism_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -1.0, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
Autism_bar_plot

bottom_row <- plot_grid(Autism_heatmap[[4]], nrow = 1, ncol = 1, scale = 0.95)
top_row <- plot_grid(Autism_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
Autism_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="tl", rel_heights = c(0.05,0.9))
Autism_final

# now add the title
#ast <- ggdraw() +draw_text(text = " *",x =0.1,y=-0.08) + theme(plot.margin = margin(0,0,0.01,0.01))
title <- ggdraw() + 
  draw_label(
    "Autism spectrum disorder",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
Autism_final_v2 <- plot_grid(
  title,Autism_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
Autism_final_v2
ggsave(filename = "Autism_celltype_trend_bonf_v3.eps",plot = Autism_final_v2,width = 10.5, height = 18.5, units = "cm")

ggsave(filename = "Autism_celltype_trend_bonf_v3.pdf",plot = Autism_final_v2,width = 10.5, height = 18.5, units = "cm")


# CBLM 
CBLM_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes_v2$CBLM, slot = "scale.data")
CBLM_heatmap_data_t <- t(CBLM_heatmap_data)

CBLM_heatmap <- pheatmap(CBLM_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 5, fontsize = 4, clustering_method = "ward.D2",legend = T,border_color = NA)

CBLM_results <- GSEA_results %>% filter(Disease == "CBLM") 

level_order <- c('Endothelial','Rhombic lip','Granule cells-1','Granule cells-2','Glutamatergic DCN','GC precusors','Ciliated cells','Unknown','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
CBLM_bar_plot <-  ggplot(CBLM_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
CBLM_bar_plot <- CBLM_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
CBLM_bar_plot <- CBLM_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
CBLM_bar_plot <- CBLM_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
CBLM_bar_plot <- CBLM_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
CBLM_bar_plot

bottom_row <- plot_grid(CBLM_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.945)
top_row <- plot_grid(CBLM_bar_plot, nrow = 1, ncol = 1, scale = 0.985)
CBLM_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
CBLM_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.34,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Cerebellar malformations",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
CBLM_final_v2 <- plot_grid(
  title,CBLM_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
CBLM_final_v2
ggsave(filename = "CBLM_celltype_trend_bonf_v3.eps",plot = CBLM_final_v2,width = 8.0, height = 11, units = "cm")

ggsave(filename = "CBLM_celltype_trend_bonf_v3.pdf",plot = CBLM_final_v2,width = 8.0, height = 11, units = "cm")





"HEREDATAXGREEN"
HAG_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes_v2$HEREDATAXGREEN, slot = "scale.data")
HAG_heatmap_data_t <- t(HAG_heatmap_data)

HAG_heatmap <- pheatmap(HAG_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 4, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

HAG_results <- GSEA_results %>% filter(Disease == "HEREDATAXGREEN") 

level_order <- c('Ciliated cells','Endothelial','Choroid plexus','Roof plate','Purkinje neurons','Unknown','Granule cells-1','Granule cells-2','Glutamatergic DCN','Bergmann glia','GC precusors','Rhombic lip')
HAG_bar_plot <-  ggplot(HAG_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
HAG_bar_plot <- HAG_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
HAG_bar_plot <- HAG_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
HAG_bar_plot <- HAG_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
HAG_bar_plot <- HAG_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
HAG_bar_plot

bottom_row <- plot_grid(HAG_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.92)
top_row <- plot_grid(HAG_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
HAG_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
HAG_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "HEREDATAXGREEN",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
HAG_final_v2 <- plot_grid(
  title,HAG_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
HAG_final_v2
ggsave(filename = "HEREDATAXGREEN_celltype_trend_bonf_v3.eps",plot = HAG_final_v2,width = 6.5, height = 19, units = "cm")

ggsave(filename = "HEREDATAXGREEN_celltype_trend_bonf_v3.pdf",plot = HAG_final_v2,width = 6.5, height = 19, units = "cm")


muultipanel_A <- plot_grid(CBLM_final_v2,SCA_final_v2, align = "v", nrow = 2, rel_heights = c(0.7,0.45) )
muultipanel_A

panel_B <- plot_grid(muultipanel_A,Autism_final_v2, align="h", ncol = 2, rel_widths = c(0.7,0.7))
panel_B

ggsave(filename = "disease_risk_main.eps",plot = panel_B, width = 18.5, height = 20, units = "cm")

ggsave(filename = "disease_risk_main.pdf",plot = panel_B, width = 18.5, height = 20, units = "cm")


# ID
ID_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes_v2$ID, slot = "scale.data")
ID_heatmap_data_t <- t(ID_heatmap_data)

ID_heatmap <- pheatmap(ID_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 2.0, cutree_rows = 5, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

ID_results <- GSEA_results %>% filter(Disease == "ID") 

level_order <- c('Granule cells-1','Granule cells-2','GC precusors','Glutamatergic DCN','Rhombic lip','Unknown','Endothelial','Ciliated cells','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
ID_bar_plot <-  ggplot(ID_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
ID_bar_plot <- ID_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
ID_bar_plot <- ID_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
ID_bar_plot <- ID_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
ID_bar_plot <- ID_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0, 0, -1.0, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
ID_bar_plot

bottom_row <- plot_grid(ID_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.945)
top_row <- plot_grid(ID_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
ID_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="tl", rel_heights = c(0.05,0.9))
ID_final

# now add the title
#ast <- ggdraw() +draw_text(text = "* *         *",x =0.35,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Intellectual disability",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
ID_final_v2 <- plot_grid(
  title,ID_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
ID_final_v2

ggsave(filename = "ID_celltype_trend_bonf_v3.eps",plot = ID_final_v2,width = 10.5, height = 19.5, units = "cm")

ggsave(filename = "ID_celltype_trend_bonf_v3.pdf",plot = ID_final_v2,width = 10.5, height = 19.5, units = "cm")





muultipanel_A <- plot_grid(CBLM_final_v2,SCA_final_v2, align = "v", nrow = 2, rel_heights = c(0.6,0.45) )
muultipanel_A

panel_B <- plot_grid(Autism_final_v2, ID_final_v2, align="h", ncol = 2, rel_widths = c(0.6,0.6), axis = "tlrb", scale = c(0.99,0.99))
panel_B


panel_C <- plot_grid(muultipanel_A,panel_B,align = "h", ncol = 2,rel_widths = c(0.35,0.65), axis = "tlrb", scale = c(0.975,1))
panel_C
ggsave(filename = "disease_risk_main.eps",plot = panel_C, width = 19, height = 20.5, units = "cm")

ggsave(filename = "disease_risk_main.pdf",plot = panel_C, width = 19, height = 20.5, units = "cm")


# JOUBERT
joubert_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes_v2$JOUBERT, slot = "scale.data")
joubert_heatmap_data_t <- t(joubert_heatmap_data)

joubert_heatmap <- pheatmap(joubert_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 3, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

joubert_results <- GSEA_results %>% filter(Disease == "JOUBERT") 

level_order <- c('Ciliated cells','Granule cells-1','Granule cells-2','Purkinje neurons','Unknown','Endothelial','Choroid plexus','Roof plate','Glutamatergic DCN','Rhombic lip','Bergmann glia','GC precusors')
joubert_bar_plot <-  ggplot(joubert_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
joubert_bar_plot <- joubert_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
joubert_bar_plot <- joubert_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
joubert_bar_plot <- joubert_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
joubert_bar_plot <- joubert_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
joubert_bar_plot

bottom_row <- plot_grid(joubert_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.94)
top_row <- plot_grid(joubert_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
joubert_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
joubert_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Joubert syndrome",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
joubert_final_v2 <- plot_grid(
  title,joubert_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
joubert_final_v2

ggsave(filename = "joubert_celltype_trend_bonf_v3.eps",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")

ggsave(filename = "joubert_celltype_trend_bonf_v3.pdf",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")



# ALZ
ALZ_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes_v2$ALZNEG, slot = "scale.data")
ALZ_heatmap_data_t <- t(ALZ_heatmap_data)

ALZ_heatmap <- pheatmap(ALZ_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 4, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

ALZ_results <- GSEA_results %>% filter(Disease == "ALZNEG") 

level_order <- c('Unknown','Ciliated cells','Purkinje neurons','Roof plate','Choroid plexus','Endothelial','Granule cells-1','Granule cells-2','Glutamatergic DCN','Bergmann glia','GC precusors','Rhombic lip')
ALZ_bar_plot <-  ggplot(ALZ_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
ALZ_bar_plot <- ALZ_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
ALZ_bar_plot <- ALZ_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(p.adj)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
ALZ_bar_plot <- ALZ_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
ALZ_bar_plot <- ALZ_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
ALZ_bar_plot

bottom_row <- plot_grid(ALZ_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.92)
top_row <- plot_grid(ALZ_bar_plot, nrow = 1, ncol = 1, scale = 0.98)
ALZ_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
ALZ_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Alzheimer's disease",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
ALZ_final_v2 <- plot_grid(
  title,ALZ_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
ALZ_final_v2

ggsave(filename = "ALZ_celltype_trend_bonf_v3.eps",plot = ALZ_final_v2,width = 6.5, height = 11.5, units = "cm")

ggsave(filename = "ALZ_celltype_trend_bonf_v3.pdf",plot = ALZ_final_v2,width = 6.5, height = 11.5, units = "cm")













Autism_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$AUTISM, slot = "scale.data")
Autism_heatmap_data_t <- t(Autism_heatmap_data)

Autism_heatmap <- pheatmap(Autism_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 5, fontsize = 3, clustering_method = "ward.D2",border_color = NA)

Autism_results <- GSEA_results %>% filter(Disease == "AUTISM") 

level_order <- c('Glutamatergic DCN', 'Granule cells-2', 'Endothelial','Granule cells-1','GC precusors','Rhombic lip','Ciliated cells','Unknown','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
Autism_bar_plot <-  ggplot(Autism_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
Autism_bar_plot <- Autism_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
Autism_bar_plot <- Autism_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(FDR)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
Autism_bar_plot <- Autism_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
Autism_bar_plot <- Autism_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -1.0, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
Autism_bar_plot

bottom_row <- plot_grid(Autism_heatmap[[4]], nrow = 1, ncol = 1, scale = 0.95)
top_row <- plot_grid(Autism_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
Autism_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="tl", rel_heights = c(0.05,0.9))
Autism_final

# now add the title
#ast <- ggdraw() +draw_text(text = " *",x =0.1,y=-0.08) + theme(plot.margin = margin(0,0,0.01,0.01))
title <- ggdraw() + 
  draw_label(
    "Autism spectrum disorder",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
Autism_final_v2 <- plot_grid(
  title,Autism_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
Autism_final_v2
ggsave(filename = "Autism_celltype_trend_bonf.eps",plot = Autism_final_v2,width = 10.5, height = 18.5, units = "cm")

ggsave(filename = "Autism_celltype_trend_bonf.pdf",plot = Autism_final_v2,width = 10.5, height = 18.5, units = "cm")


# CBLM 
CBLM_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$CBLM, slot = "scale.data")
CBLM_heatmap_data_t <- t(CBLM_heatmap_data)

CBLM_heatmap <- pheatmap(CBLM_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 5, fontsize = 4, clustering_method = "ward.D2",legend = T,border_color = NA)

CBLM_results <- GSEA_results %>% filter(Disease == "CBLM") 

level_order <- c('Endothelial','Rhombic lip','Granule cells-1','Granule cells-2','Glutamatergic DCN','GC precusors','Ciliated cells','Unknown','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
CBLM_bar_plot <-  ggplot(CBLM_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
CBLM_bar_plot <- CBLM_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
CBLM_bar_plot <- CBLM_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(FDR)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
CBLM_bar_plot <- CBLM_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
CBLM_bar_plot <- CBLM_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
CBLM_bar_plot

bottom_row <- plot_grid(CBLM_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.945)
top_row <- plot_grid(CBLM_bar_plot, nrow = 1, ncol = 1, scale = 1)
CBLM_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
CBLM_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.34,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Cerebellar malformations",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
CBLM_final_v2 <- plot_grid(
  title,CBLM_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
CBLM_final_v2
ggsave(filename = "CBLM_celltype_trend_bonf.eps",plot = CBLM_final_v2,width = 8.0, height = 11, units = "cm")

ggsave(filename = "CBLM_celltype_trend_bonf.pdf",plot = CBLM_final_v2,width = 8.0, height = 11, units = "cm")





# SCA
SCA_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$SCA, slot = "scale.data")
SCA_heatmap_data_t <- t(SCA_heatmap_data)

SCA_heatmap <- pheatmap(SCA_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 4, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

SCA_results <- GSEA_results %>% filter(Disease == "SCA") 

level_order <- c('Unknown','Bergmann glia','Purkinje neurons','Granule cells-1','Granule cells-2','GC precusors','Rhombic lip','Choroid plexus','Roof plate','Ciliated cells','Glutamatergic DCN','Endothelial')
SCA_bar_plot <-  ggplot(SCA_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
SCA_bar_plot <- SCA_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
SCA_bar_plot <- SCA_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(FDR)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
SCA_bar_plot <- SCA_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
SCA_bar_plot <- SCA_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
SCA_bar_plot

bottom_row <- plot_grid(SCA_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.945)
top_row <- plot_grid(SCA_bar_plot, nrow = 1, ncol = 1, scale = 1)
SCA_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
SCA_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Spinocerebellar ataxia",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
SCA_final_v2 <- plot_grid(
  title,SCA_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
SCA_final_v2
ggsave(filename = "SCA_celltype_trend.eps",plot = SCA_final_v2,width = 6.5, height = 10, units = "cm")

ggsave(filename = "SCA_celltype_trend.pdf",plot = SCA_final_v2,width = 6.5, height = 10, units = "cm")


muultipanel_A <- plot_grid(CBLM_final_v2,SCA_final_v2, align = "v", nrow = 2, rel_heights = c(0.7,0.45) )
muultipanel_A

panel_B <- plot_grid(muultipanel_A,Autism_final_v2, align="h", ncol = 2, rel_widths = c(0.7,0.7))
panel_B

ggsave(filename = "disease_risk_main.eps",plot = panel_B, width = 18.5, height = 20, units = "cm")

ggsave(filename = "disease_risk_main.pdf",plot = panel_B, width = 18.5, height = 20, units = "cm")


# ID
ID_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$ID, slot = "scale.data")
ID_heatmap_data_t <- t(ID_heatmap_data)

ID_heatmap <- pheatmap(ID_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 2.0, cutree_rows = 5, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

ID_results <- GSEA_results %>% filter(Disease == "ID") 

level_order <- c('Granule cells-1','Granule cells-2','GC precusors','Glutamatergic DCN','Rhombic lip','Unknown','Endothelial','Ciliated cells','Choroid plexus','Roof plate','Bergmann glia','Purkinje neurons')
ID_bar_plot <-  ggplot(ID_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
ID_bar_plot <- ID_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
ID_bar_plot <- ID_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(FDR)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
ID_bar_plot <- ID_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
ID_bar_plot <- ID_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0, 0, -1.0, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
ID_bar_plot

bottom_row <- plot_grid(ID_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.945)
top_row <- plot_grid(ID_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
ID_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="tl", rel_heights = c(0.05,0.9))
ID_final

# now add the title
#ast <- ggdraw() +draw_text(text = "* *         *",x =0.35,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Intellectual disability",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
ID_final_v2 <- plot_grid(
  title,ID_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
ID_final_v2

muultipanel_A <- plot_grid(CBLM_final_v2,SCA_final_v2, align = "v", nrow = 2, rel_heights = c(0.6,0.45) )
muultipanel_A

panel_B <- plot_grid(Autism_final_v2, ID_final_v2, align="h", ncol = 2, rel_widths = c(0.6,0.6), axis = "tlrb", scale = c(0.99,0.99))
panel_B


panel_C <- plot_grid(muultipanel_A,panel_B,align = "h", ncol = 2,rel_widths = c(0.35,0.65), axis = "tlrb",scale = c(0.975,1))
panel_C
ggsave(filename = "disease_risk_main_bonf.eps",plot = panel_C, width = 19, height = 20.5, units = "cm")

ggsave(filename = "disease_risk_main_bonf.pdf",plot = panel_C, width = 19, height = 20.5, units = "cm")


# JOUBERT
joubert_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$JOUBERT, slot = "scale.data")
joubert_heatmap_data_t <- t(joubert_heatmap_data)

joubert_heatmap <- pheatmap(joubert_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 3, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

joubert_results <- GSEA_results %>% filter(Disease == "JOUBERT") 

level_order <- c('Ciliated cells','Granule cells-1','Granule cells-2','Purkinje neurons','Unknown','Endothelial','Choroid plexus','Roof plate','Glutamatergic DCN','Rhombic lip','Bergmann glia','GC precusors')
joubert_bar_plot <-  ggplot(joubert_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
joubert_bar_plot <- joubert_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
joubert_bar_plot <- joubert_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(FDR)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
joubert_bar_plot <- joubert_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
joubert_bar_plot <- joubert_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
joubert_bar_plot

bottom_row <- plot_grid(joubert_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.94)
top_row <- plot_grid(joubert_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
joubert_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
joubert_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Joubert syndrome",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
joubert_final_v2 <- plot_grid(
  title,joubert_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
joubert_final_v2

ggsave(filename = "joubert_celltype_trend_bonf.eps",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")

ggsave(filename = "joubert_celltype_trend_bonf.pdf",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")



# ALZ
ALZ_heatmap_data <- FetchData(Avg_organoid,vars = disease_risk_genes$ALZ, slot = "scale.data")
ALZ_heatmap_data_t <- t(ALZ_heatmap_data)

ALZ_heatmap <- pheatmap(ALZ_heatmap_data_t, treeheight_row = 0,treeheight_col = 0, fontsize_col = 4,fontsize_row = 3, cutree_rows = 4, fontsize = 3, clustering_method = "ward.D2",legend = T,border_color = NA)

ALZ_results <- GSEA_results %>% filter(Disease == "ALZ") 

level_order <- c('Unknown','Ciliated cells','Purkinje neurons','Roof plate','Choroid plexus','Endothelial','Granule cells-1','Granule cells-2','Glutamatergic DCN','Bergmann glia','GC precusors','Rhombic lip')
ALZ_bar_plot <-  ggplot(ALZ_results ,aes(x= factor(Celltype, levels = level_order), y= bon_log10_padj, fill=direction)) + geom_bar(stat="identity", width = 0.5) + theme_bw() 
ALZ_bar_plot <- ALZ_bar_plot + xlab("") + theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank()) 
ALZ_bar_plot <- ALZ_bar_plot + geom_hline(yintercept = -log10(0.05), linetype= "dashed",colour="black",size=0.5) + 
  theme(legend.position = "right") + scale_fill_manual(values=c("Up"="red","Down"= "dodgerblue4")) + ylab("-log10(FDR)") + theme(axis.title.y  = element_text(size=4), axis.text.y = element_text(size=3))
ALZ_bar_plot <- ALZ_bar_plot+ theme(legend.text = element_text(size=3), legend.title = element_blank(), panel.border  = element_blank()) 
ALZ_bar_plot <- ALZ_bar_plot +theme(axis.line = element_line(color = 'black', size=0.3)) +theme(plot.margin = unit(c(0.0, 0, -0.8, 0), "cm")) + theme(legend.key.size = unit(0.3, "cm"))
ALZ_bar_plot

bottom_row <- plot_grid(ALZ_heatmap[[4]], nrow = 1, ncol = 1,scale = 0.93)
top_row <- plot_grid(ALZ_bar_plot, nrow = 1, ncol = 1, scale = 0.99)
ALZ_final <- plot_grid(top_row,bottom_row, align = "v",nrow = 2, axis="l", rel_heights = c(0.05,0.9))
ALZ_final

# now add the title
#ast <- ggdraw() +draw_text(text = "*",x =0.725,y=-0.08) + theme(plot.margin = margin(0.0,0,0,0.01))
title <- ggdraw() + 
  draw_label(
    "Alzheimer's disease",
    fontface = 'bold', size=6,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
ALZ_final_v2 <- plot_grid(
  title,ALZ_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.05, 1)
)
ALZ_final_v2

ggsave(filename = "ALZ_celltype_trend_bonf.eps",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")

ggsave(filename = "ALZ_celltype_trend_bonf.pdf",plot = joubert_final_v2,width = 6.5, height = 10, units = "cm")



panel_D <- plot_grid(joubert_final_v2,ALZ_final_v2,  align = "h", ncol = 2,rel_widths =c(0.5,0.5), axis = "tlrb")
panel_D


ggsave(filename = "supp_disease_risk_trend_bonf.eps",plot = panel_D,width = 14, height = 13, units = "cm")

ggsave(filename = "supp_disease_risk_trend_bonf.pdf",plot = panel_D,width = 14, height = 13, units = "cm")

