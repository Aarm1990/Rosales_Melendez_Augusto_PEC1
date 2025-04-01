#Codigo para PEC 1

#Importar dataset de Caquexia.
caquexia<-read.csv("human_cachexia.csv")

#Exploración del dataset
head(caquexia)

#Instalación de paquete summarized experiments
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install("SummarizedExperiment", force = TRUE)



#Cargar paquete summarized experiment
library(BiocManager)
library(SummarizedExperiment)

#Extraer metadatos (Id de paciente y condicion del paciente)

metadatos<-caquexia[,1:2]
colnames(metadatos)<-c("Id_paciente", "Condicion")

#Extraer resultados de experimento (Metabolitos).
metabolitos<-log(caquexia[,-c(1:2)]) # aplicar transformacion log a los datos
# Se aplica transformacion log porque para lograr normalizar los datos.


# Transponer la tabla de metabolitos 
metabolitos <- t(as.matrix(metabolitos)) 
colnames(metabolitos)<-c(metadatos$Id_paciente)# asignar nombres a las columnas ( tomado de id del paciente.)



#Creación de objeto Summarized experiments

objeto_caquexia <- SummarizedExperiment(
  assays = list(counts = as.matrix(metabolitos)), # asignar resultados
  colData = DataFrame(metadatos) # asignar los metadatos.
)


#Exportando objeto a tipo Rda
save(objeto_caquexia, file="objeto_caquexia.Rda")


#Análisis exploratorio de los datos


# Boxplot de los metabolitos 
par(mar = c(10, 4, 4, 2)) # Configurar el limite inferior del grafico.
"Hacer un boxplot de los metabolitos"
boxplot(t(assay(objeto_caquexia)), # Graficar los valores de los metabolitos
        outline = FALSE, 
        las = 2, # rotar nombres en el eje y.
        main = "Boxplot por metabolito", 
        ylab = "Intensidad",
        cex.axis = 0.7  
        )

"Al inspeccionar el boxplot de los metabolitos se determina que despues de la transfo
mación logaritmica siguen una distribución normal"

" Los metabolitos con mayor concentracion son la creatinina y el de menor concentracion 
el fumarato"

# se obtendra el sumario estadistico de la creatinina y el fumarato para el analisis 
# descriptivo

summary(log(caquexia$Creatine)) # Medidas de dispersion de la creatina
sd(log(caquexia$Creatine)) # Desviacion estandar
summary(log(caquexia$Fumarate)) # Medidas de dispersion del fumarato.
sd(log(caquexia$Fumarate))# Desviacion estandar del fumarato


"Realización del analisis de compoentes principales"

#Transponer matriz para la PCA
matrix<-t(assay(objeto_caquexia))
#Estandarizar datos y aplicar PCA
pca<-prcomp(matrix,scale. = TRUE)
summary(pca)

# Extraer la condicion de los metadatos
condicion <- colData(objeto_caquexia)$Condicion

# Crear gráfico con resultados del PCA
plot(pca$x[,1], pca$x[,2], # Seleccionar las dos columnas que corresponden al pca1 y pca 2.
     col = ifelse(condicion == "cachexic", "red", "blue"),# colorear los perfiles metabolicos segun su condicion.
     pch = 19,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"), # asignar nombres basados en la importancia de los metabolitos del PCA 1
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)"), # Mismo paso que arriba.
     main = "PCA de pacientes según perfiles metabólicos")
legend("topright", legend = c("Cachexic", "Control"), col = c("red", "blue"), pch = 19)


# Representar los componentes del PCA
componentes <- pca$rotation

# Seleccionar top 10 metabolitos que más contribuyen a PC1
importancia_PC1 <- sort(abs(componentes[,1]), decreasing = TRUE)[1:10]

# Gráfico de barras de los 10 metabolitos que mas contribuyen al PCA 1
barplot(importancia_PC1,
        las = 2,
        col = "steelblue",
        main = "Contribución al PC1 (Top 10 metabolitos)",
        ylab = "Valor absoluto del loading")


"Análisis comparativo: Diferencia de medias de metabolitos entre grupos"

# Extraer los metadatos y valores de los metabolitos.
metabolitos <- assay(objeto_caquexia)         # Matriz de metabolitos
grupo <- colData(objeto_caquexia)$Condicion   # Vector de grupo (cachexic / control)

# Crear dataframe para guardar resultados de la T de student.
resultados <- data.frame(
  Metabolito = rownames(metabolitos),
  p_value = NA, # Llenar con valores nuelos
  media_cachexic = NA,
  media_control = NA
)

# Loop por cada metabolito para comparar medias entre grupos 
for (i in 1:nrow(metabolitos)) { #Bucle que por cada metabolito extrae los valores y los compara segun grupos.
  valores <- metabolitos[i, ] # Iterar por cada metabolito
  grupo1 <- valores[grupo == "cachexic"] # asignar cada valor si tiene caquexia.
  grupo2 <- valores[grupo == "control"]
  
  prueba <- t.test(grupo1, grupo2) # Hacer la T de student.
  
  resultados$p_value[i] <- prueba$p.value # Guardar valor de p en el dataframe resultados
  resultados$media_cachexic[i] <- mean(grupo1, na.rm = TRUE) # Guardar las medias del grupo de caquexia
  resultados$media_control[i] <- mean(grupo2, na.rm = TRUE) # Guardar las medias del grupo control.
}

# Ajustar valor de p con comparaciones multiples FDR
resultados$FDR <- p.adjust(resultados$p_value, method = "fdr") 

# Ordenar por p-valor en orden de significancia estadistica.
resultados <- resultados[order(resultados$p_value), ]

# Ver los top 10 de valores de p.
head(resultados, 10)




