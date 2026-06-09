pacman::p_load("WGCNA")

options(stringsAsFactors = FALSE)
disableWGCNAThreads()

load("GSE160306_processed.RData")

# Preparamos los datos para preservación de módulos
# De referencia utilizamos Macula y como test la periferia 
multiExpr <- list(
  Macula = list(data = t(expr_macula)),
  Periphery = list(data = t(expr_periphery_clean))
)

multiColor <- list(
  Macula = net_macula$colors
)

# Calcular preservación de módulos
set.seed(42)
mp <- modulePreservation(
  multiData = multiExpr,
  multiColor = multiColor,
  dataIsExpr = TRUE,
  referenceNetworks = 1,
  nPermutations = 200,
  randomSeed = 42,
  quickCor = 0,
  verbose = 3
)

# Extraer estadísticos principales: Zsummary y medianRank
module_names <- rownames(mp$preservation$Z$ref.Macula$inColumnsAlsoPresentIn.Periphery)
Zsummary <- mp$preservation$Z$ref.Macula$inColumnsAlsoPresentIn.Periphery[, "Zsummary.pres"]
medianRank <- mp$preservation$observed$ref.Macula$inColumnsAlsoPresentIn.Periphery[, "medianRank.pres"]

# Crear tabla de resultados ordenada por Zsummary
preservation_df <- data.frame(
  Module = module_names,
  Zsummary = Zsummary,
  medianRank = medianRank,
  ModuleSize = as.numeric(table(net_macula$colors)[module_names])
)
preservation_df <- preservation_df[order(-preservation_df$Zsummary), ]

# Clasificar módulos por preservación
# Zsummary > 10: preservado | 2-10: débil | < 2: no preservado
preservation_df$Status <- ifelse(preservation_df$Zsummary > 10, "Preserved",
                                 ifelse(preservation_df$Zsummary > 2, "Weak", "Not preserved"))

table(preservation_df$Status)
preservation_df[preservation_df$Status == "Weak", ]

# Visualización
png(file.path(img_dir, "module_preservation.png"), width = 800, height = 600)
plot(preservation_df$ModuleSize, preservation_df$Zsummary,
     col = preservation_df$Module,
     pch = 19, cex = 1.5,
     xlab = "Module size",
     ylab = "Zsummary",
     main = "Module preservation: Macula vs Periphery",
     log = "x")
abline(h = 10, col = "red", lty = 2)
abline(h = 2, col = "blue", lty = 2)
text(preservation_df$ModuleSize, preservation_df$Zsummary,
     labels = preservation_df$Module,
     cex = 0.6, pos = 3)
dev.off()

# Guardar resultados
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
     file = "GSE160306_processed.RData")
