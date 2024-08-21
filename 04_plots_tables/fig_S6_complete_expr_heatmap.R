###############################################################
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
###############################################################


###############################################################
tile_height <- unit(base_text_size, "pt")
tile_width <- unit(0.6, "pt")
ht_opt("simple_anno_size" = tile_height)
ht_opt("heatmap_row_names_gp" = text_setting_small_it)
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"
names(cluster_colors) <- c(1:4)
###############################################################


###############################################################
gene_order <- c("POSTN", "NRP2", "UTP20", "CEACAM5", "MRAP2", "SERPINB2", "FAM83D", "SEC14L1", "TLE1", "TFF3", "DPYSL3", "STOM", "SPAG5", "MFSD2A", "AADAC", "KRT6A", "ABO", "DHX35", "GATA2", "KBTBD4", "FGFR1", "SNCA", "TMEM200A", "GALNT5", "CCL26", "EIF4E2", "GAD1", "ADM", "RPP25", "IGF2BP3", "AKAP12", "TXNL4A", "DOK1", "GNPNAT1", "UPK3A", "PDE10A", "RECQL4", "ECT2", "BFSP1", "PPTC7", "EGFL6", "NTRK2", "ELOVL5", "CLCA1", "CPT1A", "MTNR1B", "HOMER2", "RAMP1", "PAWR", "NDUFA8", "GTPBP8", "ZBTB7C", "B3GNT6", "SLC37A1", "KPNA3", "SYAP1", "DCAF12", "MUC12", "NEK6", "ABHD2", "PXDN", "SPRR1B", "UBE2T", "MACC1", "MELK", "FANCI", "CD44", "ANLN", "TSPAN5", "CEP72", "CKAP2", "NEDD9", "KLF4", "HYLS1", "DEFB1", "SPATS2", "IL18R1", "BIRC5", "PHLDB2", "AMMECR1", "NOP16", "HRASLS", "HSD11B2", "DDX49", "SERPINB13", "S100A10", "FA2H", "KRT23", "LEO1", "ST6GAL1", "NACC1", "FKBP5", "HIGD2A", "TSPAN3", "SLC39A8", "NF2", "MS4A2", "UBE2J1", "POC1B", "CTBP2", "TXNDC17", "RUNX2", "ARNTL2", "UPK1B", "NCAPG", "SURF2", "C14orf166", "EFHD2", "SIGLEC6", "TNPO3", "SIGLEC8", "GRK5", "TCN1", "ZMAT4", "NAV1", "TMEM64", "CTSG", "REXO2", "FAM120A", "FAM110C", "ZNF467", "KCNE3", "DENND4A", "PNPO", "PTGS1", "VSIG2", "SH3RF2", "MOCS1", "RNF10", "KRT10", "GALNT1", "H2AFX", "CLC", "LRRC31", "AGPS", "PLCB4", "FBXL7", "ANO7", "CCNB2", "CPSF3", "CD109", "SLC18A2", "TRAK1", "TRIM59", "LRRC8A", "GSR", "KRT6B", "HPGDS", "CYP2R1", "DQX1", "SLC24A3", "CEP152", "CPA3", "SLC44A1", "GPATCH4", "DLG4", "PTTG1", "TM4SF1", "FBN2", "CDC42EP5", "CDH26", "TBC1D7", "TOP2A", "ALOX15", "MXD1", "NOP10", "GDF2", "CST2", "HPCAL1", "KLF9", "SLC6A8", "PIM1", "AHSG", "C12orf49", "MYADM", "C12orf57", "PMCH", "LRRC8D", "CST1", "ZMAT2", "ITLN1", "FHL1", "S100A16", "ANGPT4", "GSN", "KIT", "EME1", "TBC1D16", "UCHL3", "SERPINB5", "MYLK3", "TAOK3", "PTPRH", "XRCC4", "ETV1", "PKP2", "BDNF", "SCGB3A1", "C3", "SCGB1A1", "ATP13A4", "SLC5A1", "CYP2A13", "SLC19A2", "ST8SIA4", "ZC3H6", "VEGFA", "MUC5B", "SLC39A10", "ZNF418", "CAPRIN2", "EPB41L3", "FABP6", "SGCE", "ZNF274", "FOLR1", "CACHD1", "WNT5B", "SLC6A16", "CPNE8", "UBE2E2", "ZBTB46", "LRRK2", "TMEM45A", "DCDC2", "RPRML", "CYP2J2", "EZH1", "MYO9A", "SCNN1B", "CDK14", "SEC14L3", "KCNMB2", "CXCL1", "SLC4A4", "RAB3IP", "CWH43", "SCN4B", "NPNT", "HNMT", "FLNB", "MAOB", "SEC63", "CNTN3", "FXYD6", "DAPK1", "ADIPOR2", "SLC22A4", "PTPRM", "SHISA2", "STXBP1", "EFHD1", "C16orf89", "SLC27A2", "BAZ2B", "ITGB5", "RPRM", "LIMCH1", "COL9A2", "GLS", "SLITRK6", "APCDD1", "C1orf52", "PIP", "CYP4X1", "FMN2", "INPP5B", "FAM46C", "CLDN8", "SLAMF7", "RBPJ", "RSAD2", "GPR155", "MGP", "DPYSL2", "RFTN1", "EPB41L2", "ZNF382", "ANTXR2", "PKD2", "STEAP2", "SPAG6", "C2CD4B", "ATL2", "PRKAR2B", "SULF1", "HEY1", "CD47", "ZNF331", "CNTD1", "GLI3", "DPY19L3", "SUSD4", "HEG1", "PTGFR", "FHDC1", "KLK10", "ADAM28", "CCDC81", "CHST9", "ABCA13", "OSBPL6", "INTS10", "INO80", "HYI", "NCOA7", "NOTCH2", "SLC13A3", "PPP1R3B", "RTN4RL1", "EFNB2", "IRF2BP2", "HTRA1", "RGS17", "ABCC4", "TIGD7", "GSTA4", "RAB6B", "SLC41A1", "CMTM8", "AFAP1L1", "IRS2", "CHGB", "GPX3", "CAPN13", "WNK4", "AKAP11", "USP3", "RTN1", "FAM117A", "RTN4", "EPM2AIP1", "KIAA0319L", "ZHX2", "ILDR1", "EPHX3", "FHOD3", "LYPD6B", "CA8", "P4HA3", "TRIM55", "ZIK1", "ANXA6", "ZMAT1", "FNIP2", "C6", "GAB2", "ITPR1", "SCNN1G", "GPD1L", "TGFBR3", "TM2D3", "STEAP4", "SLC15A2", "GRP", "HGSNAT", "PLAG1", "FRAS1", "KCNA1")
combined_pheno$Severity <- ordered(combined_pheno$Severity, levels = c("Healthy", "Mild/moderate", "Severe", "Unspecified"))
combined_pheno$Sex <- ordered(combined_pheno$Sex, levels = c("female", "male"))
combined_pheno <- combined_pheno[with(combined_pheno, order(Cluster, Severity, Sex, Dataset)), ]
combined_expr <- combined_expr[gene_order, rownames(combined_pheno)]
combined_expr <- as.matrix(t(scale(t(combined_expr))))
combined_expr[combined_expr > 2.5] <- 2.5
combined_expr[combined_expr < -2.5] <- -2.5
###############################################################


###############################################################
min_expr <- -2.5
max_expr <- 2.5
expr_matrix_colour_pts <- c(min_expr, min_expr/2, 0, max_expr/2, max_expr)
#expr_cols <- colorRamp2(expr_matrix_colour_pts, rev(c("#762A83", "#C2A5CF", "#F7F7F7", "#ACD39E", "#1B7837")))
expr_cols <- colorRamp2(expr_matrix_colour_pts, scico(5, palette = "vik"))
expr_matrix_legend_pts <- c(min_expr, max_expr)
expr_matrix_legend_labels <- c(round(min_expr, 2), round(max_expr, 2))
expr_legend <- Legend(col_fun = expr_cols,
  at = expr_matrix_legend_pts,
  labels = expr_matrix_legend_labels,
  title = "Row-scaled\nexpression", title_position = "topcenter",
  grid_height = unit(base_text_size / 3, "pt"),
  legend_width = 4 * unit(base_text_size, "pt"),
  background = "transparent",
  title_gp = text_setting_small,
  title_gap = gap_size, gap = gap_size,
  labels_gp = text_setting_sub_small,
  direction = "horizontal")

create_discrete_legend <- function(legend_colours, legend_title, num_row = 1)
{
  discrete_legend <- Legend(labels = names(legend_colours),
    legend_gp = gpar(fill = legend_colours),
    title = legend_title, title_position = "leftcenter",
    grid_height = unit(base_text_size, "pt"),
    grid_width = unit(base_text_size, "pt"),
    background = "transparent",
    title_gp = text_setting_small,
    title_gap = gap_size, gap = gap_size,
    labels_gp = text_setting_sub_small,
    nrow = num_row)
}

cluster_legend <- create_discrete_legend(cluster_colors, "Cluster")
severity_legend <- create_discrete_legend(severity_colours, "Severity")
dataset_legend <- create_discrete_legend(dataset_colours, "Dataset", 2)
sex_legend <- create_discrete_legend(sex_colors, "Sex", 2)

ht_mp_annotation <- HeatmapAnnotation(
  Cluster = combined_pheno$Cluster, Severity = combined_pheno$Severity, Sex = combined_pheno$Sex, Dataset = combined_pheno$Dataset,
  col = list(Cluster = cluster_colors, Severity = severity_colours, Dataset = dataset_colours, Sex = sex_colors),
  show_legend = FALSE,
  annotation_name_gp = text_setting_small,
  annotation_name_side = "left")

ht_mp_expr <- Heatmap(combined_expr,
  col = expr_cols, na_col = "#FFEE99",
  show_column_names = FALSE,
  width = ncol(combined_expr) * tile_width,
  name = "expr_matrix", column_title = "Sample",
  column_title_side = "bottom",# column_title_gp = text_setting_small,
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  #column_split = as.factor(combined_pheno$Cluster),
  #show_column_dend = FALSE,
  gap = gap_size, top_annotation = ht_mp_annotation)
###############################################################


###############################################################
cairo_pdf(paste0(figures_folder, "fig_S6_conorm_expr_heatmap.pdf"), width = 7.2, height = 51.5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 1000, ncol = 100)))

pushViewport(viewport(layout.pos.row = 15:999, layout.pos.col = 2:99))
draw(ht_mp_expr, height = length(gene_order) * tile_height, gap = gap_size, background = "transparent", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 2:12, layout.pos.col = 2:12))
draw(expr_legend)
upViewport(1)

pushViewport(viewport(layout.pos.row = 2:6, layout.pos.col = 14:38))
draw(cluster_legend)
upViewport(1)

pushViewport(viewport(layout.pos.row = 2:6, layout.pos.col = 44:95))
draw(severity_legend)
upViewport(1)

pushViewport(viewport(layout.pos.row = 8:13, layout.pos.col = 14:81))
draw(dataset_legend)
upViewport(1)

pushViewport(viewport(layout.pos.row = 8:13, layout.pos.col = 83:98))
draw(sex_legend)
upViewport(1)

dev.off()
###############################################################

