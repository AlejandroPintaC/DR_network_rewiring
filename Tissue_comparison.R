
pacman::p_load("WGCNA", "ggplot2", "igraph", "aricode", "clusterProfiler", "org.Hs.eg.db", "reshape2")

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

# Tenemos que comparar las particiones de módulos control de macula y periferia

partition_macula <- net_macula_control$colors
partition_periphery <- net_periphery_control$colors

#Checar que los genes estén en el mismo orden
all(names(partition_macula) == names(partition_periphery))

partition_macula_factor <- as.factor(partition_macula)
partition_periphery_factor <- as.factor(partition_periphery)

# Verificar que no hay NAs
sum(is.na(partition_macula_factor))
sum(is.na(partition_periphery_factor))

# NMI entre dos particiones 
nmi_control <- NMI(partition_macula_factor, partition_periphery_factor)
nmi_control

# Adjusted Rand Index (ARI) para segundo punto de comparasión
ari_control <- ARI(partition_macula_factor, partition_periphery_factor)
ari_control

# Función para obtener enriquecimiento GO en módulos
get_all_go <- function(net_colors, module_list){
  results <- list()
  for(mod in module_list){
    if(mod == "grey") next  # excluir módulo grey
    genes <- names(net_colors[net_colors == mod])
    entrez <- tryCatch(
      bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) NULL
    )
    if(is.null(entrez) || nrow(entrez) < 5) next
    
    go_result <- tryCatch(
      enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05),
      error = function(e) NULL
    )
    
    if(!is.null(go_result) && nrow(go_result) > 0){
      results[[mod]] <- go_result@result$Description
    }
  }
  return(results)
}

# Obtener GO de todos los módulos (Macula control)
modules_macula_control <- unique(net_macula_control$colors)
go_macula_control <- get_all_go(net_macula_control$colors, modules_macula_control)

length(go_macula_control)

# Obtener GO de modulos (Periferia control)
modules_periphery_control <- unique(net_periphery_control$colors)
go_periphery_control <- get_all_go(net_periphery_control$colors, modules_periphery_control)

length(go_periphery_control)

# Similitud funcional a nivel de red completa (Jaccard)
all_functions_macula <- unique(unlist(go_macula_control))
all_functions_periphery <- unique(unlist(go_periphery_control))

jaccard_functional <- length(intersect(all_functions_macula, all_functions_periphery)) / 
  length(union(all_functions_macula, all_functions_periphery))

jaccard_functional

# Cuántas funciones únicas tiene cada tejido y cuántas comparten
length(all_functions_macula)
length(all_functions_periphery)
length(intersect(all_functions_macula, all_functions_periphery))

# Revisar si go_macula_control y go_periphery_control son realmente diferentes
identical(go_macula_control, go_periphery_control)

# Ver las primeras funciones de cada uno
head(all_functions_macula, 5)
head(all_functions_periphery, 5)

# Matriz de similitud funcional módulo a módulo (Jaccard)
modules_macula_with_go <- names(go_macula_control)
modules_periphery_with_go <- names(go_periphery_control)

similarity_matrix <- matrix(0, 
                            nrow = length(modules_macula_with_go), 
                            ncol = length(modules_periphery_with_go),
                            dimnames = list(modules_macula_with_go, modules_periphery_with_go))

for(i in modules_macula_with_go){
  for(j in modules_periphery_with_go){
    funcs_i <- go_macula_control[[i]]
    funcs_j <- go_periphery_control[[j]]
    jaccard_ij <- length(intersect(funcs_i, funcs_j)) / length(union(funcs_i, funcs_j))
    similarity_matrix[i, j] <- jaccard_ij
  }
}

dim(similarity_matrix)

# Encontrar los pares con mayor similitud funcional
similarity_long <- melt(similarity_matrix)
colnames(similarity_long) <- c("Modulo_Macula", "Modulo_Periferia", "Jaccard")

# Ordenar de mayor a menor similitud
similarity_long <- similarity_long[order(-similarity_long$Jaccard), ]

head(similarity_long, 15)

# Ver funciones específicas que comparten turquoise-turquoise
intersect(go_macula_control[["turquoise"]], go_periphery_control[["turquoise"]])[1:15]

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
     file = "GSE160306_processed.RData")

