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
