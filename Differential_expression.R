pacman::p_load("limma", "WGCNA", "ggplot2")

options(stringsAsFactors = FALSE)

load("GSE160306_processed.RData")

# Preparamos la  metadata de macula con disease_group como factor
metadata_macula <- metadata[macula_samples, ]
metadata_macula$disease_group <- factor(metadata_macula$`disease_group:ch1`,
                                        levels = c("Control", "Diabetic", "NPDR", "PDR"))

table(metadata_macula$disease_group)

# Matriz de diseño
design_macula <- model.matrix(~0 + disease_group, data = metadata_macula)
colnames(design_macula) <- levels(metadata_macula$disease_group)

head(design_macula)

# Ajustamos el modelo lineal
fit_macula <- lmFit(expr_macula, design_macula)

# Definimos contrastes de interés o comparaciones entre grupos
contrast_matrix <- makeContrasts(
  Diabetic_vs_Control = Diabetic - Control,
  NPDR_vs_Control = NPDR - Control,
  PDR_vs_Control = PDR - Control,
  NPDR_vs_Diabetic = NPDR - Diabetic,
  PDR_vs_NPDR = PDR - NPDR,
  levels = design_macula
)

# Aplicamos los contrastes y calcular estadísticos
fit_macula_contrasts <- contrasts.fit(fit_macula, contrast_matrix)
fit_macula_eb <- eBayes(fit_macula_contrasts)

# Resumen de genes diferencialmente expresados por contraste
summary(decideTests(fit_macula_eb, p.value = 0.05))

design_periphery <- model.matrix(~0 + disease_group, data = metadata_periphery_dge)
colnames(design_periphery) <- levels(metadata_periphery_dge$disease_group)

head(design_periphery)

# Preparamos metadata de periferia
metadata_periphery_dge <- metadata_periphery_clean
metadata_periphery_dge$disease_group <- factor(metadata_periphery_dge$`disease_group:ch1`,
                                               levels = c("Control", "Diabetic", "NPDR", "PDR"))

table(metadata_periphery_dge$disease_group)

# Construcción de matriz de diseño para el modelo lineal de expresión diferencial
# ~0 + disease_group: sin intercepto, una columna por grupo, permitiendo el contraste explícito entre grupos
design_periphery <- model.matrix(~0 + disease_group, data = metadata_periphery_dge)
colnames(design_periphery) <- levels(metadata_periphery_dge$disease_group)

head(design_periphery)

# Ajustamos modelo lineal
fit_periphery <- lmFit(expr_periphery_clean, design_periphery)

# Definimos contrastes
contrast_matrix_periphery <- makeContrasts(
  Diabetic_vs_Control = Diabetic - Control,
  NPDR_vs_Control = NPDR - Control,
  PDR_vs_Control = PDR - Control,
  NPDR_vs_Diabetic = NPDR - Diabetic,
  PDR_vs_NPDR = PDR - NPDR,
  levels = design_periphery
)

# Aplicamos contrastes y calculamos estadísticos
fit_periphery_contrasts <- contrasts.fit(fit_periphery, contrast_matrix_periphery)
fit_periphery_eb <- eBayes(fit_periphery_contrasts)

# Resumen de DEGs por contraste
summary(decideTests(fit_periphery_eb, p.value = 0.05))

# IMPORTANTE!!!1: Los contrastes PDR_vs_Control y PDR_vs_NPDR en periferia usan n=1 para el grupo PDR (tras exclusión del outlier sample_78).
# limma no puede estimar varianza intra-grupo con n=1, por lo que estos resultados (864 y 837 DEGs respectivamente) son artefactos estadísticos
# y NO deben interpretarse como hallazgos biológicos reales.
# Se reportan únicamente como limitación metodológica del dataset!!!!.

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
     file = "GSE160306_processed.RData")



