BiocManager::install("WGCNA")
# Carga de librerias
library(WGCNA)
library(GEOquery)
library(ggplot2)

# Configuración recomendada para WGCNA
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
disableWGCNAThreads()
# Carga de datos
load("GSE160306_processed.RData")
ls()

# Verificar que no haya genes con demasiados valores faltantes o varianza cero
# Verificar calidad de datos para WGCNA (genes macula)
# WGCNA necesita muestras en filas y genes en columnas
gsg_macula <- goodSamplesGenes(t(expr_macula), verbose = 3)
gsg_macula$allOK

# Same, pero genes de periferia 
gsg_periphery <- goodSamplesGenes(t(expr_periphery_clean), verbose = 3)
gsg_periphery$allOK

# Soft Thresholding Power (parametro clave que determina como se pondera las correlaciones entre los genes)
# Selección de soft thresholding power - Macula
powers <- c(1:20)

sft_macula <- pickSoftThreshold(t(expr_macula), 
                                powerVector = powers,
                                verbose = 5,
                                networkType = "signed")


#Tenemos que verificar que el slope sea negativo (red scale-free)
png(file.path(img_dir, "soft_threshold_macula.png"), width = 900, height = 500)
par(mfrow = c(1,2))

# R^2
plot(sft_macula$fitIndices[,1], 
     -sign(sft_macula$fitIndices[,3])*sft_macula$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale independence - Macula",
     type = "n")
text(sft_macula$fitIndices[,1], 
     -sign(sft_macula$fitIndices[,3])*sft_macula$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.90, col = "blue")

# Conectividad media
plot(sft_macula$fitIndices[,1], 
     sft_macula$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity - Macula",
     type = "n")
text(sft_macula$fitIndices[,1], 
     sft_macula$fitIndices[,5],
     labels = powers, col = "red")
dev.off()
# Ahora STP para periferia
powers2 <- c(1:30)
sft_periphery2 <- pickSoftThreshold(t(expr_periphery_clean), 
                                    powerVector = powers2,
                                    verbose = 5,
                                    networkType = "signed")
png(file.path(img_dir, "soft_threshold_periphery.png"), width = 900, height = 500)
par(mfrow = c(1,2))

plot(sft_periphery2$fitIndices[,1], 
     -sign(sft_periphery2$fitIndices[,3])*sft_periphery2$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale independence - Periphery",
     type = "n")
text(sft_periphery2$fitIndices[,1], 
     -sign(sft_periphery2$fitIndices[,3])*sft_periphery2$fitIndices[,2],
     labels = powers2, col = "red")
abline(h = 0.90, col = "blue")

plot(sft_periphery2$fitIndices[,1], 
     sft_periphery2$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity - Periphery",
     type = "n")
text(sft_periphery2$fitIndices[,1], 
     sft_periphery2$fitIndices[,5],
     labels = powers2, col = "red")
dev.off()

save(expr_macula, expr_periphery_clean,
     metadata, metadata_periphery_clean,
     periphery_samples_clean, macula_samples,
     sft_macula, sft_periphery2,
     file = "GSE160306_processed.RData")

# Poderes seleccionados:
# Macula: power = 10 (R^2 = 0.903, mean connectivity = 154)
# Periferia: power = 18 (R^2 = 0.887, mean connectivity = 12.2)
# Nota: para periferia se usó umbral R^2 >= 0.85 por n < 50