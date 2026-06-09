pacman::p_load("WGCNA", "ggplot2", "org.Hs.db")

options(stringsAsFactors = FALSE)
disableWGCNAThreads()

load("GSE160306_processed.RData")

## Macula ##
# Calcular Module Membership (kME) de Macula
# kME mide que tan correlacionado está cada gen con el eigengene de su modulo
# Genes con kME alto son hub genes 
kME_macula <- signedKME(t(expr_macula), MEsmacula)

# Ver dimensiones
dim(kME_macula)
head(kME_macula[, 1:5])


# Función para extraer hub genes de un módulo (kME_threshold = 0.7)
get_hub_genes <- function(kME_matrix, module_colors, module_name, kME_threshold = 0.7){
  col_name <- paste0("kME", module_name)
  module_genes <- names(module_colors[module_colors == module_name])
  kME_values <- kME_matrix[module_genes, col_name]
  hub_idx <- kME_values > kME_threshold
  
  if(sum(hub_idx) == 0) return(data.frame(Gene = character(), kME = numeric(), Module = character()))
  
  return(data.frame(
    Gene = module_genes[hub_idx],
    kME = kME_values[hub_idx],
    Module = module_name
  ))
}

# Hub genes módulos significativos - Macula
modules_sig_macula <- c("green", "purple", "tan", "salmon")

hub_genes_macula <- data.frame()
for(mod in modules_sig_macula){
  result <- get_hub_genes(kME_macula, net_macula$colors, mod)
  if(nrow(result) > 0){
    hub_genes_macula <- rbind(hub_genes_macula, result)
  }
}

table(hub_genes_macula$Module)

# Convertir Ensembl IDs a símbolos de genes
hub_symbols <- bitr(hub_genes_macula$Gene,
                    fromType = "ENSEMBL",
                    toType = "SYMBOL",
                    OrgDb = org.Hs.eg.db)

# Unir con tabla de hub genes
hub_genes_macula <- merge(hub_genes_macula, hub_symbols,
                          by.x = "Gene", by.y = "ENSEMBL")

# Ordenar por módulo y kME descendente
hub_genes_macula <- hub_genes_macula[order(hub_genes_macula$Module, 
                                           -hub_genes_macula$kME), ]

# Ver los top 5 por módulo
do.call(rbind, lapply(modules_sig_macula, function(mod){
  head(hub_genes_macula[hub_genes_macula$Module == mod, ], 5)
}))

## Periferia ##
# Calcular kME - Periferia
kME_periphery <- signedKME(t(expr_periphery_clean), MEsperiphery)

# Hub genes módulos significativos - Periferia
modules_sig_periphery <- c("turquoise", "grey60", "saddlebrown")

hub_genes_periphery <- data.frame()
for(mod in modules_sig_periphery){
  result <- get_hub_genes(kME_periphery, net_periphery$colors, mod)
  if(nrow(result) > 0){
    hub_genes_periphery <- rbind(hub_genes_periphery, result)
  }
}

table(hub_genes_periphery$Module)

# Convertir Ensembl IDs a símbolos - Periferia
hub_symbols_periphery <- bitr(hub_genes_periphery$Gene,
                              fromType = "ENSEMBL",
                              toType = "SYMBOL",
                              OrgDb = org.Hs.eg.db)

# Unir con tabla de hub genes
hub_genes_periphery <- merge(hub_genes_periphery, hub_symbols_periphery,
                             by.x = "Gene", by.y = "ENSEMBL")

# Ordenar por módulo y kME descendente
hub_genes_periphery <- hub_genes_periphery[order(hub_genes_periphery$Module,
                                                 -hub_genes_periphery$kME), ]

# Ver top 5 por módulo
do.call(rbind, lapply(modules_sig_periphery, function(mod){
  head(hub_genes_periphery[hub_genes_periphery$Module == mod, ], 5)
}))

# KEGG pathway enrichment de los módulos significativos de Macula
modules_sig_macula <- c("green", "purple", "tan", "salmon")

for(mod in modules_sig_macula){
  genes <- names(net_macula$colors[net_macula$colors == mod])
  entrez <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  kegg_result <- enrichKEGG(gene = entrez$ENTREZID,
                            organism = "hsa",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
  
  if(!is.null(kegg_result) && nrow(kegg_result) > 0){
    png(file.path(img_dir, paste0("KEGG_", mod, "_macula.png")), width = 900, height = 700)
    print(dotplot(kegg_result, showCategory = 15,
                  title = paste0("KEGG - ME", mod, " Macula")))
    dev.off()
  }
  
  assign(paste0("kegg_", mod, "_macula"), kegg_result)
}

# KEGG pathway enrichment de los  módulos significativos de Periferia
modules_sig_periphery <- c("turquoise", "grey60", "saddlebrown")

for(mod in modules_sig_periphery){
  genes <- names(net_periphery$colors[net_periphery$colors == mod])
  entrez <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  kegg_result <- enrichKEGG(gene = entrez$ENTREZID,
                            organism = "hsa",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
  
  if(!is.null(kegg_result) && nrow(kegg_result) > 0){
    png(file.path(img_dir, paste0("KEGG_", mod, "_periphery.png")), width = 900, height = 700)
    print(dotplot(kegg_result, showCategory = 15,
                  title = paste0("KEGG - ME", mod, " Periphery")))
    dev.off()
  }
  
  assign(paste0("kegg_", mod, "_periphery"), kegg_result)
}

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
     file = "GSE160306_processed.RData")
