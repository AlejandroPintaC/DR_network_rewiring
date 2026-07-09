pacman::p_load("limma", "WGCNA", "ggplot2", "clusterProfiler", "org.Hs.eg.db")

options(stringsAsFactors = FALSE)
disableWGCNAThreads()

load("GSE160306_processed.RData")

#  Comparación de  PDR Macula vs PDR Periferia


# Filtramos muestras PDR (Macula)
pdr_macula_samples <- rownames(metadata[metadata$`disease_group:ch1` == "PDR" &
                                          rownames(metadata) %in% macula_samples, ])
expr_macula_pdr <- expr_macula[, pdr_macula_samples]

# Filtrar muestras PDR (Periferia)
pdr_periphery_samples <- rownames(metadata_periphery_clean[metadata_periphery_clean$`disease_group:ch1` == "PDR", ])
expr_periphery_pdr <- expr_periphery_clean[, pdr_periphery_samples]

dim(expr_macula_pdr)
dim(expr_periphery_pdr)

# Comparación entre Diferencia vs Diferencia (Control → PDR) 

# Calcular expresión promedio por grupo en macula
control_macula_mean <- rowMeans(expr_macula[, control_macula_samples])
pdr_macula_mean <- rowMeans(expr_macula[, pdr_macula_samples])

# Delta macula: cambio promedio Control → PDR
delta_macula <- pdr_macula_mean - control_macula_mean

# Calcular expresión promedio por grupo en periferia
control_periphery_mean <- rowMeans(expr_periphery_clean[, control_periphery_samples])

# Ver distribución del delta
summary(delta_macula)

# Delta periferia: cambio promedio Control → PDR
pdr_periphery_mean <- as.numeric(expr_periphery_clean[, pdr_periphery_samples])
names(pdr_periphery_mean) <- rownames(expr_periphery_clean)

delta_periphery <- pdr_periphery_mean - control_periphery_mean

summary(delta_periphery)

# Correlación entre ambos deltas (¿van en la misma dirección?)
cor_deltas <- cor(delta_macula, delta_periphery, method = "pearson")
cor_deltas

# Scatter plot delta macula vs delta periferia
png(file.path(img_dir, "delta_macula_vs_periphery.png"), width = 800, height = 700)
plot(delta_macula, delta_periphery,
     xlab = "Delta Macula (PDR - Control, log2-CPM)",
     ylab = "Delta Periferia (PDR - Control, log2-CPM)",
     main = paste0("Comparación D: Cambio Control→PDR\nMacula vs Periferia (r=", round(cor_deltas, 3), ")"),
     pch = 19, cex = 0.4, col = "steelblue", alpha = 0.5)
abline(h = 0, v = 0, col = "grey", lty = 2)
abline(lm(delta_periphery ~ delta_macula), col = "red", lwd = 2)
dev.off()

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
     partition_macula_factor, partition_periphery_factor,
     nmi_control, ari_control,
     go_macula_control, go_periphery_control,
     similarity_matrix, similarity_long,
     metadata_macula, design_macula, fit_macula_eb,
     metadata_periphery_dge, design_periphery, fit_periphery_eb,
     pdr_macula_samples, expr_macula_pdr,
     delta_macula, delta_periphery, cor_deltas,
     file = "GSE160306_processed.RData")

