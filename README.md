# **Análisis expresión diferencial rna-seq**
### 0 horas vs 72 horas
### Basado en el estudio: Diferenciación de células madre embrionarias a endodermo
### Experimento: E-MTAB-9194 (Expression Atlas)
### Comparación: embryonic stem cell (0h) vs definitive endoderm cell (72h)

## 0) Instalacion y carga de paquetes necesarios 
```r
library(BiocManager)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(EnhancedVolcano)
#install.packages("pheatmap")
library(pheatmap)      # Para heatmaps de correlación
library(RColorBrewer)  # Paletas de colores
#los que no descargué ya los había descargado antes 
```
### Set working directory
```r
setwd("C:\\Users\\Laura Gonzalez\\Desktop\\Proyecto_DEseq2_GonzalezLaura")
getwd()
```
## 1) CARGA DE DATOS DESDE EXPRESSION ATLAS

### Se usa read.delim porque los archivos están en formato tabulado (.tsv)
```r
countss <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-9194/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
metadatas <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-9194/resources/ExperimentDesignFile.RnaSeq/experiment-design")
```
## 2) SELECCIÓN DE MUESTRAS DE INTERÉS (0h y 72h), ya que sin esta selección inicial se hacia muy dificil el tratamiento de los datos nulos 

### - ERR4235451, ERR4235464, ERR4235465: 0h (réplicas 1,  2, 3)
### - ERR4235461, ERR4235462, ERR4235463: 72h (réplicas 1, 2, 3)
```r
counts <- countss[, c("Gene.ID", "Gene.Name",
                      "ERR4235451", "ERR4235461", "ERR4235462",
                      "ERR4235463", "ERR4235464", "ERR4235465")]
head(counts)

metadata <- metadatas[metadatas$Run %in% c("ERR4235451", "ERR4235461", 
                                           "ERR4235462", "ERR4235463",
                                           "ERR4235464", "ERR4235465"), ]
head(metadata)
```
## 3) ACOMODAR LOS DATOS AL FORMATO QUE DESeq2 ESPERA

# 3.1) Poner IDs de gen como rownames
rownames(counts) <- counts$Gene.ID

# 3.2) Guardar info de genes
genes <- counts[, c("Gene.ID", "Gene.Name")]

# 3.3) Dejar solo columnas numéricas de conteos
counts <- counts[, -c(1, 2)]

# 3.4) Preparar metadata
rownames(metadata) <- metadata$Run
metadata <- metadata[, "Factor.Value.time.", drop = FALSE]
colnames(metadata) <- "tiempo"

# 3.5) Limpiar etiquetas
metadata$tiempo[metadata$tiempo == "0 hour"] <- "0"
metadata$tiempo[metadata$tiempo == "72 hour"] <- "72"

# 3.6) Declarar factor con orden
metadata$tiempo <- factor(metadata$tiempo, levels = c("0", "72"))

## 4) VERIFICACIÓN: gen LZTS1

gene_id <- genes$Gene.ID[genes$Gene.Name == "LZTS1"]
gene_counts <- counts[gene_id, ]
gene_data <- cbind(metadata, counts = as.numeric(gene_counts))

ggplot(gene_data, aes(x = tiempo, y = counts, fill = tiempo)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(title = "Conteos crudos de LZTS1 por tiempo") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

## 5) CREAR OBJETO DESEQ2 Y CORRER ANÁLISIS

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ tiempo)

# Filtrar genes con baja expresión
dds <- dds[rowSums(counts(dds)) > 10, ]

# Ejecutar análisis
dds <- DESeq(dds)

# Contraste 72 vs 0
res <- results(dds, contrast = c("tiempo", "72", "0"), alpha = 0.01)

## 6) PCA PARA CONTROL DE CALIDAD

vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "tiempo", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = tiempo, shape = tiempo)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = name), vjust = -1.5, size = 3.5, fontface = "bold") +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  theme_minimal()

# 3) Acomodar los datos al formato que DESeq2 espera
#    - Filas de 'counts' = genes (con sus IDs en rownames)
#    - Columnas de 'counts' = muestras (los nombres deben coincidir con metadata)
#    - Rownames de 'metadata' = IDs de muestra
# 3.1) Poner IDs de gen como nombres de fila
rownames(counts) <- counts$Gene.ID
head(counts)
# 3.2) Guardar información de genes para referencias posteriores
genes <- counts[, c("Gene.ID", "Gene.Name")]
head(genes)
# 3.3) Eliminar columnas no numéricas de counts (dejar solo conteos)
counts <- counts[, -c(1, 2)]
head(counts)
# 3.4) Preparar metadata
rownames(metadata) <- metadata$Run
metadata <- metadata[, "Factor.Value.time.", drop = FALSE]
colnames(metadata) <- "tiempo"
head(metadata)
# 3.5) Limpiar etiquetas de tiempo
metadata$tiempo[metadata$tiempo == "0 hour"] <- "0"
metadata$tiempo[metadata$tiempo == "72 hour"] <- "72"
head(metadata)
# 3.6) Declarar 'tiempo' como factor y fijar el orden (referencia primero)
#      - Esto define cómo se interpretan los contrastes (comparaciones)
metadata$tiempo <- factor(metadata$tiempo, levels = c("0", "72"))
metadata$tiempo

# 4) VERIFICACIÓN RÁPIDA: Gen marcador LZTS1
gene_id <- genes$Gene.ID[genes$Gene.Name == "LZTS1"]
gene_counts <- counts[gene_id, ]
gene_data <- cbind(metadata, counts = as.numeric(gene_counts))

# Boxplot de expresión cruda
p1 <- ggplot(gene_data, aes(x = tiempo, y = counts, fill = tiempo)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(title = "Conteos crudos de LZTS1 por tiempo",
       x = "Tiempo", y = "Conteos (sin normalizar)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

print(p1)

# 5) Crear el objeto DESeq y correr el análisis
#    - design = ~ tiempo indica que queremos probar diferencias por tiempo
#    - Filtramos genes muy poco expresados (ruido) para evitar falsos positivos
# 5.1) Crear objeto DESeq
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ tiempo)
# Filtrar genes con baja suma de conteos (umbral > 10)
dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

# 5.4) Extraer resultados del contraste
res <- results(dds, contrast = c("tiempo", "72", "0"), alpha = 0.01)
head(res)

# 6) CONTROL DE CALIDAD: PCA
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "tiempo", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p2 <- ggplot(pca_data, aes(PC1, PC2, color = tiempo, shape = tiempo)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text(aes(label = name), vjust = -1.5, size = 3.5, fontface = "bold") +
  xlab(paste0("PC1: ", percentVar[1], "% varianza (Diferenciación celular)")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza (Variación residual)")) +
  labs(title = "PCA: Separación entre 0h y 72h",
       subtitle = "PC1 captura toda la variabilidad biológica",
       caption = "Datos transformados con VST | n=500 genes más variables") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    legend.position = "right"
  ) +
  scale_color_manual(values = c("0" = "#E41A1C", "72" = "#377EB8"),
                     name = "Tiempo",
                     labels = c("0h (ESC)", "72h (Endodermo)")) +
  scale_shape_manual(values = c("0" = 16, "72" = 17),
                     name = "Tiempo",
                     labels = c("0h (ESC)", "72h (Endodermo)"))

print(p2)

# 7) CONTROL DE CALIDAD: CORRELACIÓN DE SPEARMAN Y CLUSTERING

normalized_counts <- counts(dds, normalized = TRUE)

mean_expression <- rowMeans(normalized_counts)
high_expr_genes <- names(mean_expression[mean_expression > 10])

cat(length(high_expr_genes))

cor_matrix <- cor(normalized_counts[high_expr_genes, ], method = "spearman")
print(round(cor_matrix, 3))

annotation_col <- data.frame(
  Tiempo = metadata$tiempo,
  row.names = colnames(normalized_counts)
)

pheatmap(cor_matrix,
         annotation_col = annotation_col,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("white", "lightblue", "darkblue"))(50),
         main = "Correlación de Spearman entre réplicas\n(Clustering Jerárquico)",
         fontsize = 10,
         display_numbers = TRUE,
         number_format = "%.2f")

# 8) FILTRADO ESTRICTO SEGÚN CRITERIOS DEL PAPER

res_df_all <- as.data.frame(res)
res_df_all <- merge(res_df_all, genes, by.x = "row.names", by.y = "Gene.ID")
colnames(res_df_all)[1] <- "Gene.ID"

res_df_clean <- res_df_all[!is.na(res_df_all$padj) & 
                             !is.na(res_df_all$log2FoldChange), ]

res_df_filtered <- res_df_clean[
  abs(res_df_clean$log2FoldChange) >= log2(1.5) & 
    res_df_clean$padj < 0.01, 
]

res_df_filtered <- res_df_filtered[order(res_df_filtered$padj), ]

print(res_df_filtered[1:20, c("Gene.Name", "log2FoldChange", "padj")])

# 9) VERIFICACIÓN DE GENES MARCADORES
genes_a_verificar <- c("EOMES", "SOX17", "GATA4", "LZTS1")
marcadores <- res_df_filtered[res_df_filtered$Gene.Name %in% genes_a_verificar, ]

if(nrow(marcadores) > 0) {
  
  print(marcadores[, c("Gene.Name", "log2FoldChange", "padj", "baseMean")])
} else {
  
  marcadores_all <- res_df_clean[res_df_clean$Gene.Name %in% genes_a_verificar, ]
  print(marcadores_all[, c("Gene.Name", "log2FoldChange", "padj", "baseMean")])
}

# 10) VISUALIZACIONES

plotMA(res, alpha = 0.01, 
       main = "MA plot: Expresión media vs log2 Fold Change",
       ylim = c(-10, 10))

p3 <- EnhancedVolcano(res,
                      lab = genes$Gene.Name[match(rownames(res), genes$Gene.ID)],
                      x = "log2FoldChange",
                      y = "padj",
                      title = "Volcano plot: 72h vs 0h",
                      subtitle = "Genes diferencialmente expresados (FC≥1.5, padj<0.01)",
                      pCutoff = 0.01,
                      FCcutoff = log2(1.5),
                      pointSize = 2.0,
                      labSize = 4.0,
                      legendPosition = 'right',
                      legendLabSize = 12,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.3,
                      xlim = c(-14, 14),
                      ylim = c(0, 200))
print(p3)

# 11) ANOTACIÓN GENÓMICA CON biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

atributos <- c("ensembl_gene_id", "chromosome_name", 
               "start_position", "end_position")

todos_los_genes <- getBM(attributes = atributos, 
                         values = list(ensembl_gene_id = c()), 
                         mart = ensembl)

colnames(todos_los_genes)[1] <- "Gene.ID"

# 12) GUARDAR RESULTADOS
datos_unidos_filtrados <- merge(todos_los_genes, res_df_filtered, by = "Gene.ID")
datos_unidos_filtrados$chromosome_name <- paste0("chr", 
                                                 datos_unidos_filtrados$chromosome_name)

datos_unidos_completos <- merge(todos_los_genes, res_df_clean, by = "Gene.ID")
datos_unidos_completos$chromosome_name <- paste0("chr", 
                                                 datos_unidos_completos$chromosome_name)

subset_marcadores <- datos_unidos_completos[
  datos_unidos_completos$Gene.Name %in% genes_a_verificar, 
]

write.csv(datos_unidos_filtrados, "deseq_filtrado_FC1.5_padj0.01.csv", 
          row.names = FALSE)
write.csv(datos_unidos_completos, "deseq_completo.csv", 
          row.names = FALSE)
write.csv(subset_marcadores, "deseq_marcadores.csv", 
          row.names = FALSE)

write.csv(normalized_counts, "normalized_counts.csv", 
          row.names = TRUE)

write.csv(cor_matrix, "correlation_matrix_spearman.csv", 
          row.names = TRUE)

cat("- CSVs con resultados (ver lista arriba)\n")
cat("- Gráficos: PCA, MA plot, Volcano plot, Heatmap de correlación\n\n")
