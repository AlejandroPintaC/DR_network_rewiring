library(GEOquery)
gse160306 <- getGEO("GSE160306", GSEMatrix = TRUE)

#Veamos estructura general 
str(gse160306)

#Ver metadatos de las muestras 
pData(gse160306[[1]])

#Ver numero de muestras por grupo
metadata <- pData(gse160306[[1]])
table(metadata$`disease_group:ch1`)
table(metadata$`disease_group:ch1`, metadata$`sample_site:ch1`)

expr_matrix <- exprs(gse160306[[1]])
dim(expr_matrix)
head(expr_matrix[, 1:5])

# Descarga de archivo suplementarios
getGEOSuppFiles("GSE160306")

#Lectura del archivo suplementario
# Leer el archivo
expr_data <- read.table("GSE160306/GSE160306_human_retina_DR_totalRNA_normalized_cpm.txt.gz", 
                        header = TRUE, 
                        sep = ",", 
                        row.names = 1)

dim(expr_data)
head(expr_data[, 1:5])

# Leer las primeras líneas crudas para ver el formato
readLines("GSE160306/GSE160306_human_retina_DR_totalRNA_normalized_cpm.txt.gz", n = 3)

# Extraer metadatos relevantes
head(metadata$`sampleID:ch1`)

# Ordenar metadatos por sampleID para que coincidan con las columnas de la matriz
rownames(metadata) <- metadata$`sampleID:ch1`
metadata <- metadata[colnames(expr_data), ]

# Verificar que están alineados
all(rownames(metadata) == colnames(expr_data))

# Resumen general de la matriz
summary(expr_data[, 1:3])

# Distribución de grupos
table(metadata$`disease_group:ch1`)

#Ver visualmente como se separan los grupos antes de entrar a WGCNA (Weighted Gene Co-expression Network analysis)
BiocManager::install("ggplot2")
install.packages("ggplot2")
library(ggplot2)

# Principal Component Analysis (PCA)
pca <- prcomp(t(expr_data), scale. = TRUE)

# Preparar datos para graficar
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  disease_group = metadata$`disease_group:ch1`,
  sample_site = metadata$`sample_site:ch1`
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = disease_group, shape = sample_site)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA - GSE160306")
#Guardar PCA
png(file.path(img_dir, "PCA_GSE160306.png"), width = 800, height = 600)
ggplot(pca_df, aes(x = PC1, y = PC2, color = disease_group, shape = sample_site)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA - GSE160306")
dev.off()

# Hay que verla varianza explicada de cada componente
summary(pca)$importance[2, 1:10]

# verificar si PC1 se asocia a una variable clinica(disease_group o sample_site)
boxplot(pca$x[,1] ~ metadata$`disease_group:ch1`, 
        main = "PC1 por grupo de enfermedad",
        xlab = "Grupo", ylab = "PC1")
png(file.path(img_dir, "PC1_disease_group.png"), width = 800, height = 600)
boxplot(pca$x[,1] ~ metadata$`disease_group:ch1`,
        main = "PC1 por grupo de enfermedad",
        xlab = "Grupo", ylab = "PC1")
dev.off()

#Boxplot pero por sitio de muestreo
png(file.path(img_dir, "PC1_sample_site.png"), width = 800, height = 600)
boxplot(pca$x[,1] ~ metadata$`sample_site:ch1`,
        main = "PC1 por sitio de muestreo",
        xlab = "Sitio", ylab = "PC1")
dev.off()

# Separar muestras por sitio
macula_samples <- rownames(metadata[metadata$`sample_site:ch1` == "Macula", ])
periphery_samples <- rownames(metadata[metadata$`sample_site:ch1` == "Periphery", ])

# Subsets de expresión
expr_macula <- expr_data[, macula_samples]
expr_periphery <- expr_data[, periphery_samples]

# Verificar dimensiones
dim(expr_macula)
dim(expr_periphery)

#Verificar cuantas muestras hay por grupo de enfermedad tenemos en cada subset
table(metadata[macula_samples, "disease_group:ch1"])
table(metadata[periphery_samples, "disease_group:ch1"])

# PCA de cada subset para confirmar que la separar por sitio ya no haya efecto dominante 
# PCA macula
pca_macula <- prcomp(t(expr_macula), scale. = TRUE)

pca_macula_df <- data.frame(
  PC1 = pca_macula$x[, 1],
  PC2 = pca_macula$x[, 2],
  disease_group = metadata[macula_samples, "disease_group:ch1"]
)

png(file.path(img_dir, "PCA_macula.png"), width = 800, height = 600)
ggplot(pca_macula_df, aes(x = PC1, y = PC2, color = disease_group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA - Macula")
dev.off()

# PCA periphery
pca_periphery <- prcomp(t(expr_periphery), scale. = TRUE)

pca_periphery_df <- data.frame(
  PC1 = pca_periphery$x[, 1],
  PC2 = pca_periphery$x[, 2],
  disease_group = metadata[periphery_samples, "disease_group:ch1"]
)

png(file.path(img_dir, "PCA_periphery.png"), width = 800, height = 600)
ggplot(pca_periphery_df, aes(x = PC1, y = PC2, color = disease_group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA - Periphery")
dev.off()

#Identificacion del ouline en la periferia
which(pca_periphery$x[, 1] > 150)

metadata["sample_78", c("disease_group:ch1", "sample_site:ch1", "rin:ch1", "age:ch1", "comorbidities:ch1")] # Esto para verficar que muestra es y a que grupo pertenece

metadata["sample_78", c("disease_group:ch1", "rin:ch1", "age:ch1", "post_mortem_interval:ch1")] # Ver resto de metadatos

# Excluir outlier de periferia
periphery_samples_clean <- periphery_samples[periphery_samples != "sample_78"]
expr_periphery_clean <- expr_data[, periphery_samples_clean]
metadata_periphery_clean <- metadata[periphery_samples_clean, ]

# Verificar
dim(expr_periphery_clean)
table(metadata_periphery_clean$`disease_group:ch1`)
pca_periphery_clean <- prcomp(t(expr_periphery_clean), scale. = TRUE)

#PCA rapido de periferia limpia para confirmar que el outlier ya no esta
pca_periphery_clean_df <- data.frame(
  PC1 = pca_periphery_clean$x[, 1],
  PC2 = pca_periphery_clean$x[, 2],
  disease_group = metadata_periphery_clean$`disease_group:ch1`
)

png(file.path(img_dir, "PCA_periphery_clean.png"), width = 800, height = 600)
ggplot(pca_periphery_clean_df, aes(x = PC1, y = PC2, color = disease_group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA - Periphery (clean)")
dev.off()
