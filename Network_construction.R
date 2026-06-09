pacman::p_load("WGCNA",
               "ggplot2")
options(stringsAsFactors = FALSE)
disableWGCNAThreads()

load("GSE160306_processed.RData")

# Ahora construimos la red (Mácula)
net_macula <- blockwiseModules(
  t(expr_macula),
  power = 10,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

# Verificar modulos encontrados
table(net_macula$colors)

# Dendograma de modulos
for(i in 1:4){
  png(file.path(img_dir, paste0("dendrogram_macula_block", i, ".png")), 
      width = 1200, height = 600)
  plotDendroAndColors(
    net_macula$dendrograms[[i]],
    net_macula$colors[net_macula$blockGenes[[i]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("Block", i, "- Macula")
  )
  dev.off()
}

# Ahora construimos la siguiente red (Periferia)
net_periphery <- blockwiseModules(
  t(expr_periphery_clean),
  power = 18,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

# Dendogramas periferia 
for(i in 1:length(net_periphery$dendrograms)){
  png(file.path(img_dir, paste0("dendrogram_periphery_block", i, ".png")), 
      width = 1200, height = 600)
  plotDendroAndColors(
    net_periphery$dendrograms[[i]],
    net_periphery$colors[net_periphery$blockGenes[[i]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("Block", i, "- Periphery")
  )
  dev.off()
}

# Tabla de modulos con colores y numero de genes para ambos tejidos
# Tabla de módulos (Macula)
module_table_macula <- data.frame(
  Modulo = names(table(net_macula$colors)),
  N_genes = as.numeric(table(net_macula$colors)),
  Tejido = "Macula"
)

# Tabla de módulos (Periferia)
module_table_periphery <- data.frame(
  Modulo = names(table(net_periphery$colors)),
  N_genes = as.numeric(table(net_periphery$colors)),
  Tejido = "Periphery"
)

# Combinar
module_table <- rbind(module_table_macula, module_table_periphery)

# Visualizar ordenado por número de genes
module_table[order(module_table$Tejido, -module_table$N_genes), ]

# Visualizacion en gráfica de colores de cada modulo para cada tejido 
# Gráfica de barras de módulos (Macula)
module_macula_plot <- module_table_macula[order(-module_table_macula$N_genes), ]
module_macula_plot <- module_macula_plot[module_macula_plot$Modulo != "grey", ]

png(file.path(img_dir, "modules_macula.png"), width = 1000, height = 600)
barplot(module_macula_plot$N_genes,
        names.arg = module_macula_plot$Modulo,
        col = module_macula_plot$Modulo,
        las = 2,
        main = "Módulos WGCNA - Macula",
        ylab = "Número de genes",
        cex.names = 0.7)
dev.off()

# Gráfica de barras de módulos (Periferia)
module_periphery_plot <- module_table_periphery[order(-module_table_periphery$N_genes), ]
module_periphery_plot <- module_periphery_plot[module_periphery_plot$Modulo != "grey", ]

png(file.path(img_dir, "modules_periphery.png"), width = 1000, height = 600)
barplot(module_periphery_plot$N_genes,
        names.arg = module_periphery_plot$Modulo,
        col = module_periphery_plot$Modulo,
        las = 2,
        main = "Módulos WGCNA - Periphery",
        ylab = "Número de genes",
        cex.names = 0.7)
dev.off()

# Contar módulos reales (excluyendo grey)
length(unique(net_macula$colors)) - 1
length(unique(net_periphery$colors)) - 1

# Calcular eigengenes (Macula)
MEsmacula <- moduleEigengenes(t(expr_macula), net_macula$colors)$eigengenes

# Calcular eigengenes (Periferia)
MEsperiphery <- moduleEigengenes(t(expr_periphery_clean), net_periphery$colors)$eigengenes

# Ver dimensiones
dim(MEsmacula)
dim(MEsperiphery)

# Convertimos "disease_group" a variable numérica ordinal (Macula)
disease_numeric <- as.numeric(factor(metadata[macula_samples, "disease_group:ch1"],
                                     levels = c("Control", "Diabetic", "NPDR", "PDR")))

table(disease_numeric)

# Ahora calculamos la correlacion entre eigengenes y disease group
# Correlación eigengenes vs disease_group (Macula)
module_trait_cor_macula <- cor(MEsmacula, disease_numeric, use = "p")
module_trait_pval_macula <- corPvalueStudent(module_trait_cor_macula, nrow(MEsmacula))

# Ver los módulos más correlacionados
module_trait_results_macula <- data.frame(
  Modulo = rownames(module_trait_cor_macula),
  Correlacion = module_trait_cor_macula[, 1],
  Pvalor = module_trait_pval_macula[, 1]
)

# Ordenar por correlación absoluta
module_trait_results_macula <- module_trait_results_macula[
  order(-abs(module_trait_results_macula$Correlacion)), ]

head(module_trait_results_macula, 10)

#Visualización con heatmap de correlaciones
png(file.path(img_dir, "module_trait_heatmap_macula.png"), width = 700, height = 900)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix = module_trait_cor_macula,
  xLabels = "Disease progression",
  yLabels = rownames(module_trait_cor_macula),
  ySymbols = rownames(module_trait_cor_macula),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = paste0(round(module_trait_cor_macula, 2), 
                      "\n(", round(module_trait_pval_macula, 3), ")"),
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Module-trait correlation - Macula"
)
dev.off()

# Convertimos "disease_group" a variable numérica ordinal (Periferia)
disease_numeric_periphery <- as.numeric(factor(
  metadata_periphery_clean$`disease_group:ch1`,
  levels = c("Control", "Diabetic", "NPDR", "PDR")))

# Correlación eigengenes vs disease_group (Periferia)
module_trait_cor_periphery <- cor(MEsperiphery, disease_numeric_periphery, use = "p")
module_trait_pval_periphery <- corPvalueStudent(module_trait_cor_periphery, nrow(MEsperiphery))

# Top módulos correlacionados
module_trait_results_periphery <- data.frame(
  Modulo = rownames(module_trait_cor_periphery),
  Correlacion = module_trait_cor_periphery[, 1],
  Pvalor = module_trait_pval_periphery[, 1]
)

module_trait_results_periphery <- module_trait_results_periphery[
  order(-abs(module_trait_results_periphery$Correlacion)), ]

head(module_trait_results_periphery, 10)

#Heatmap de periferia
png(file.path(img_dir, "module_trait_heatmap_periphery.png"), width = 700, height = 900)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix = module_trait_cor_periphery,
  xLabels = "Disease progression",
  yLabels = rownames(module_trait_cor_periphery),
  ySymbols = rownames(module_trait_cor_periphery),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = paste0(round(module_trait_cor_periphery, 2), 
                      "\n(", round(module_trait_pval_periphery, 3), ")"),
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Module-trait correlation - Periphery"
)
dev.off()

# Actualizar RData con eigengenes y correlaciones
save(expr_macula, expr_periphery_clean,
     metadata, metadata_periphery_clean,
     periphery_samples_clean, macula_samples,
     sft_macula, sft_periphery2,
     net_macula, net_periphery,
     MEsmacula, MEsperiphery,
     module_trait_cor_macula, module_trait_pval_macula,
     module_trait_cor_periphery, module_trait_pval_periphery,
     file = "GSE160306_processed.RData")
