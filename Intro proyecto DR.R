setwd("/Users/alejandropintacastro/Desktop/Diabetic retinopathy") # Este va a ser mi directory
update.packages(ask = FALSE) # Updatear todos los paquetes que ya tenemos
BiocManager::install(version = "3.22")
BiocManager::install("GEOquery")
img_dir <- "/Users/alejandropintacastro/Desktop/Diabetic retinopathy/Imagenes"
