# **Análisis expresión diferencial RNA-SEQ**

### 0 horas vs 72 horas

### Basado en el estudio: Diferenciación de células madre embrionarias a endodermo

### Experimento: E-MTAB-9194 (Expression Atlas)

### Comparación: embryonic stem cell (0h) vs definitive endoderm cell (72h)

### 0) Instalación y carga de paquetes necesarios
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

Set working directory
```r
setwd("C:\\Users\\Laura Gonzalez\\Desktop\\Proyecto_DEseq2_GonzalezLaura")
getwd()
```

### 1) Carga de datos desde Expression Atlas

Se usa read.delim porque los archivos están en formato tabulado (.tsv)
```r
countss <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-9194/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
metadatas <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-9194/resources/ExperimentDesignFile.RnaSeq/experiment-design")
```

### 2) Selección de muestras de interés (0h y 72h)

Ya que sin esta selección inicial se hacía muy difícil el tratamiento de los datos nulos

- ERR4235451, ERR4235464, ERR4235465: 0h (réplicas 1, 2, 3)
- ERR4235461, ERR4235462, ERR4235463: 72h (réplicas 1, 2, 3)
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

### 3) Preparación de datos para DESeq2

Acomodar los datos al formato que DESeq2 espera:
- Filas de 'counts' = genes (con sus IDs en rownames)
- Columnas de 'counts' = muestras (los nombres deben coincidir con metadata)
- Rownames de 'metadata' = IDs de muestra

3.1) Asignar IDs de genes como nombres de fila
```r
rownames(counts) <- counts$Gene.ID
head(counts)
```

3.2) Guardar información de genes para referencias posteriores
```r
genes <- counts[, c("Gene.ID", "Gene.Name")]
head(genes)
```

3.3) Eliminar columnas no numéricas de counts (dejar solo conteos)
```r
counts <- counts[, -c(1, 2)]
head(counts)
```

3.4) Preparar metadata: asignar IDs como rownames y extraer columna de tiempo
```r
rownames(metadata) <- metadata$Run
metadata <- metadata[, "Factor.Value.time.", drop = FALSE]
colnames(metadata) <- "tiempo"
head(metadata)
```

3.5) Limpiar etiquetas de tiempo para simplificar
```r
metadata$tiempo[metadata$tiempo == "0 hour"] <- "0"
metadata$tiempo[metadata$tiempo == "72 hour"] <- "72"
head(metadata)
```

3.6) Declarar 'tiempo' como factor y fijar el orden (referencia primero)

Esto define cómo se interpretan los contrastes en DESeq2
```r
metadata$tiempo <- factor(metadata$tiempo, levels = c("0", "72"))
metadata$tiempo
```

### 4) Verificación rápida: Gen marcador LZTS1
```r
gene_id <- genes$Gene.ID[genes$Gene.Name == "LZTS1"]
gene_counts <- counts[gene_id, ]
gene_data <- cbind(metadata, counts = as.numeric(gene_counts))
```

Boxplot de expresión cruda
```r
p1 <- ggplot(gene_data, aes(x = tiempo, y = counts, fill = tiempo)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(title = "Conteos crudos de LZTS1 por tiempo",
       x = "Tiempo", y = "Conteos (sin normalizar)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

print(p1)
```

### 5) Crear el objeto DESeq y correr el análisis

design = ~ tiempo indica que queremos probar diferencias por tiempo. Filtramos genes muy poco expresados (ruido) para evitar falsos positivos

5.1) Crear objeto DESeq
```r
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ tiempo)
```

Filtrar genes con baja suma de conteos (aquí umbral > 10 es un ejemplo simple)
```r
dds <- dds[rowSums(counts(dds)) > 10, ]
```

Ejecuta la estimación de tamaños, dispersión y el modelo (Wald por defecto)
```r
dds <- DESeq(dds)
```

5.4) Extraer resultados del contraste de interés

- contrast = c("tiempo", "72", "0") = 72 vs 0
- alpha = 0.01 define el umbral FDR (padj) para marcar significancia
- log2FoldChange se interpreta como log2(72 / 0)
```r
res <- results(dds, contrast = c("tiempo", "72", "0"), alpha = 0.01)
head(res)
```

### 6) Control de calidad: PCA

Transformación VST (Variance Stabilizing Transformation) para PCA
```r
vsd <- vst(dds, blind = FALSE)
```

PCA plot
```r
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
```

### 7) Control de calidad: Correlación de Spearman y clustering

7.1) Obtener conteos normalizados de DESeq2
```r
normalized_counts <- counts(dds, normalized = TRUE)
```

7.2) Calcular matriz de correlación de Spearman, filtramos genes con expresión media normalizada > 10
```r
mean_expression <- rowMeans(normalized_counts)
high_expr_genes <- names(mean_expression[mean_expression > 10])

cat(length(high_expr_genes))
```

Calcular correlación solo con genes altamente expresados
```r
cor_matrix <- cor(normalized_counts[high_expr_genes, ], method = "spearman")

print(round(cor_matrix, 3))
```

7.3) Heatmap de correlación con clustering jerárquico
```r
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
```

### 8) Filtrado estricto según criterios del paper

Criterios: |log2FC| ≥ log2(1.5) = 0.585 AND padj < 0.01

8.1) Convertir a data.frame y agregar nombres de genes
```r
res_df_all <- as.data.frame(res)
head(res_df_all)
head(genes)

res_df_all <- merge(res_df_all, genes, by.x = "row.names", by.y = "Gene.ID")
colnames(res_df_all)[1] <- "Gene.ID"
```

8.2) Eliminar genes con NA
```r
res_df_clean <- res_df_all[!is.na(res_df_all$padj) & 
                             !is.na(res_df_all$log2FoldChange), ]
```

8.3) Aplicar filtros estrictos según paper original
```r
res_df_filtered <- res_df_clean[
  abs(res_df_clean$log2FoldChange) >= log2(1.5) & 
    res_df_clean$padj < 0.01, 
]
```

8.4) Ordenar por significancia
```r
res_df_filtered <- res_df_filtered[order(res_df_filtered$padj), ]

print(res_df_filtered[1:20, c("Gene.Name", "log2FoldChange", "padj")]) #los 20 más significativos
```

### 9) Verificación de genes marcadores EOMES, SOX17, LZTS1 y GATA4
```r
genes_a_verificar <- c("EOMES", "SOX17", "GATA4", "LZTS1")
marcadores <- res_df_filtered[res_df_filtered$Gene.Name %in% genes_a_verificar, ]

if(nrow(marcadores) > 0) {
  print(marcadores[, c("Gene.Name", "log2FoldChange", "padj", "baseMean")])
} else {
  marcadores_all <- res_df_clean[res_df_clean$Gene.Name %in% genes_a_verificar, ]
  print(marcadores_all[, c("Gene.Name", "log2FoldChange", "padj", "baseMean")])
}
```

### 10) Visualizaciones

10.1) MA plot: muestra relación entre abundancia media y cambio de expresión
```r
plotMA(res, alpha = 0.01, 
       main = "MA plot: Expresión media vs log2 Fold Change",
       ylim = c(-10, 10)) # puntos rojos suelen ser genes significativos
```

10.2) Volcano plot
```r
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
```

### 11) Anotación genómica con biomaRt
```r
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

atributos <- c("ensembl_gene_id", "chromosome_name", 
               "start_position", "end_position")

todos_los_genes <- getBM(attributes = atributos, 
                         values = list(ensembl_gene_id = c()), 
                         mart = ensembl)

colnames(todos_los_genes)[1] <- "Gene.ID"
head(todos_los_genes)
```

### 12) Guardar resultados

12.1) Unir resultados filtrados con coordenadas genómicas
```r
datos_unidos_filtrados <- merge(todos_los_genes, res_df_filtered, by = "Gene.ID")
datos_unidos_filtrados$chromosome_name <- paste0("chr", 
                                                 datos_unidos_filtrados$chromosome_name)
head(datos_unidos_filtrados)
```

12.2) Unir TODOS los resultados con coordenadas (para referencia)
```r
datos_unidos_completos <- merge(todos_los_genes, res_df_clean, by = "Gene.ID")
datos_unidos_completos$chromosome_name <- paste0("chr", 
                                                 datos_unidos_completos$chromosome_name)
head(datos_unidos_completos)
```

12.3) Subset de genes marcadores
```r
subset_marcadores <- datos_unidos_completos[
  datos_unidos_completos$Gene.Name %in% genes_a_verificar, 
]
```

12.4) Guardar archivos CSV
```r
write.csv(datos_unidos_filtrados, "deseq_filtrado_FC1.5_padj0.01.csv", 
          row.names = FALSE)
write.csv(datos_unidos_completos, "deseq_completo.csv", 
          row.names = FALSE)
write.csv(subset_marcadores, "deseq_marcadores.csv", 
          row.names = FALSE)
```

12.5) Guardar matriz de conteos normalizados
```r
write.csv(normalized_counts, "normalized_counts.csv", 
          row.names = TRUE)
```

12.6) Guardar matriz de correlación
```r
write.csv(cor_matrix, "correlation_matrix_spearman.csv", 
          row.names = TRUE)

cat("- CSVs con resultados (ver lista arriba)\n")
cat("- Gráficos: PCA, MA plot, Volcano plot, Heatmap de correlación\n\n")
```
