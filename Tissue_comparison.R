library(WGCNA)
library(ggplot2)

options(stringsAsFactors = FALSE)
disableWGCNAThreads()

load("GSE160306_processed.RData")

# Filtrar muestras Control (Macula)
control_macula_samples <- rownames(metadata[metadata$`disease_group:ch1` == "Control" & 
                                              rownames(metadata) %in% macula_samples, ])
expr_macula_control <- expr_macula[, control_macula_samples]

dim(expr_macula_control)

# Filtrar muestras Control (Periferia)
control_periphery_samples <- rownames(metadata_periphery_clean[metadata_periphery_clean$`disease_group:ch1` == "Control", ])
expr_periphery_control <- expr_periphery_clean[, control_periphery_samples]

dim(expr_periphery_control)

# Soft Thresholding Power
powers_control <- c(1:30)
# Macula Control
sft_macula_control <- pickSoftThreshold(t(expr_macula_control),
                                        powerVector = powers_control,
                                        verbose = 5,
                                        networkType = "signed")

png(file.path(img_dir, "soft_threshold_macula_control.png"), width = 900, height = 500)
par(mfrow = c(1,2))
plot(sft_macula_control$fitIndices[,1], 
     -sign(sft_macula_control$fitIndices[,3])*sft_macula_control$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale independence - Macula Control",
     type = "n")
text(sft_macula_control$fitIndices[,1], 
     -sign(sft_macula_control$fitIndices[,3])*sft_macula_control$fitIndices[,2],
     labels = powers_control, col = "red")
abline(h = 0.85, col = "blue")
abline(h = 0.80, col = "orange", lty = 2)

plot(sft_macula_control$fitIndices[,1], 
     sft_macula_control$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity - Macula Control",
     type = "n")
text(sft_macula_control$fitIndices[,1], 
     sft_macula_control$fitIndices[,5],
     labels = powers_control, col = "red")
dev.off()

# Periferia Control
sft_periphery_control <- pickSoftThreshold(t(expr_periphery_control),
                                           powerVector = powers_control,
                                           verbose = 5,
                                           networkType = "signed")

# Construcción de red (Macula Control)
net_macula_control <- blockwiseModules(
  t(expr_macula_control),
  power = 28,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

table(net_macula_control$colors)

# Construcción de red (Periferia Control)
net_periphery_control <- blockwiseModules(
  t(expr_periphery_control),
  power = 28,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

table(net_periphery_control$colors)

save(expr_macula, expr_periphery_clean,
     metadata, metadata_periphery_clean,
     periphery_samples_clean, macula_samples,
     sft_macula, sft_periphery2,
     net_macula, net_periphery,
     MEsmacula, MEsperiphery,
     module_trait_cor_macula, module_trait_pval_macula,
     module_trait_cor_periphery, module_trait_pval_periphery,
     go_green, go_purple_macula, go_tan_macula, go_salmon_macula,
     go_turquoise_periphery, go_grey60_periphery, go_saddlebrown_periphery,
     kegg_green_macula, kegg_purple_macula, kegg_tan_macula, kegg_salmon_macula,
     kegg_turquoise_periphery, kegg_grey60_periphery, kegg_saddlebrown_periphery,
     kME_macula, kME_periphery,
     hub_genes_macula, hub_genes_periphery,
     mp, preservation_df,
     control_macula_samples, expr_macula_control,
     control_periphery_samples, expr_periphery_control,
     sft_macula_control, sft_periphery_control,
     net_macula_control, net_periphery_control,
     file = "GSE160306_processed.RData")

