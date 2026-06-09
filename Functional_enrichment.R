BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

pacman::p_load("WGCNA", "clusterProfiler", "org.Hs.eg.db", "ggplot2" )

options(stringsAsFactors = FALSE)
load("GSE160306_processed.RData")

# Obtenemos los genes del modulo green (Macula)
green_genes <- names(net_macula$colors[net_macula$colors == "green"])

# Convertir Ensembl IDs a Entrez IDs
green_entrez <- bitr(green_genes, 
                     fromType = "ENSEMBL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

head(green_entrez)
nrow(green_entrez)

# Gene Ontology enrichment (MEgreen Macula)
go_green <- enrichGO(gene = green_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",  # Biological Process
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)

# Ver resultados
head(summary(go_green)[, c("Description", "GeneRatio", "pvalue", "p.adjust")], 15)

# Visualizacion del enriquecimiento (Macula)
png(file.path(img_dir, "GO_green_macula.png"), width = 900, height = 700)
dotplot(go_green, showCategory = 15, title = "GO Biological Process - MEgreen Macula")
dev.off()

# Módulos significativos
modules_macula_sig <- c("purple", "tan", "salmon")
modules_periphery_sig <- c("turquoise", "grey60", "saddlebrown")

# Loop para macula
for(mod in modules_macula_sig){
  genes <- names(net_macula$colors[net_macula$colors == mod])
  entrez <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  go_result <- enrichGO(gene = entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  
  png(file.path(img_dir, paste0("GO_", mod, "_macula.png")), width = 900, height = 700)
  print(dotplot(go_result, showCategory = 15, 
                title = paste0("GO BP - ME", mod, " Macula")))
  dev.off()
  
  assign(paste0("go_", mod, "_macula"), go_result)
}

# Loop para periferia 
for(mod in modules_periphery_sig){
  genes <- names(net_periphery$colors[net_periphery$colors == mod])
  entrez <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  go_result <- enrichGO(gene = entrez$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  
  png(file.path(img_dir, paste0("GO_", mod, "_periphery.png")), width = 900, height = 700)
  print(dotplot(go_result, showCategory = 15,
                title = paste0("GO BP - ME", mod, " Periphery")))
  dev.off()
  
  assign(paste0("go_", mod, "_periphery"), go_result)
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
     file = "GSE160306_processed.RData")
