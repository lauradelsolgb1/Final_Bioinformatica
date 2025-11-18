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

### 4) VERIFICACIÓN RÁPIDA: Genes marcadores (EOMES, SOX17, LZTS1, GATA4)
```r
genes_a_verificar <- c("EOMES", "SOX17", "LZTS1", "GATA4")
gene_ids <- genes$Gene.ID[genes$Gene.Name %in% genes_a_verificar]
names(gene_ids) <- genes$Gene.Name[genes$Gene.Name %in% genes_a_verificar]
```
4.2) Extraer conteos crudos para estos genes
```r
genes_counts <- counts[gene_ids, ]


genes_data <- as.data.frame(t(genes_counts))
genes_data$Sample <- rownames(genes_data)
genes_data$tiempo <- metadata[rownames(genes_data), "tiempo"]
colnames(genes_data)[1:length(gene_ids)] <- names(gene_ids)

genes_long <- melt(genes_data, 
                   id.vars = c("Sample", "tiempo"),
                   variable.name = "Gene",
                   value.name = "Counts")
```
4.4) Crear boxplot con los 4 genes
```r
p1 <- ggplot(genes_long, aes(x = tiempo, y = Counts, fill = tiempo)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  facet_wrap(~Gene, scales = "free_y", ncol = 2) +
  labs(title = "Conteos crudos de genes marcadores por tiempo",
       subtitle = "EOMES, SOX17, LZTS1, GATA4",
       x = "Tiempo", 
       y = "Conteos (sin normalizar)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(strip.text = element_text(size = 11, face = "bold"),
        legend.position = "top")

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
log2 fold change (MLE): tiempo 72 vs 0 
Wald test p-value: tiempo 72 vs 0 
DataFrame with 6 rows and 6 columns
                  baseMean log2FoldChange     lfcSE      stat      pvalue        padj
                 <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
ENSG00000000003 7103.25156     -0.2256343 0.0375516  -6.00864 1.87082e-09 5.50784e-09
ENSG00000000005  243.95564     -2.8230870 0.1677082 -16.83333 1.39066e-63 1.16495e-62
ENSG00000000419 6073.28534     -0.0800895 0.0362277  -2.21073 2.70548e-02 4.60104e-02
ENSG00000000457  853.23085      0.2192106 0.0713694   3.07149 2.12990e-03 4.23214e-03
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
22724 
```

Calcular correlación solo con genes altamente expresados
```r
cor_matrix <- cor(normalized_counts[high_expr_genes, ], method = "spearman")

print(round(cor_matrix, 3))
           ERR4235451 ERR4235461 ERR4235462 ERR4235463 ERR4235464 ERR4235465
ERR4235451      1.000      0.919      0.995      0.919      0.918      0.995
ERR4235461      0.919      1.000      0.919      0.995      0.996      0.919
ERR4235462      0.995      0.919      1.000      0.919      0.918      0.995
ERR4235463      0.919      0.995      0.919      1.000      0.996      0.919
ERR4235464      0.918      0.996      0.918      0.996      1.000      0.918
ERR4235465      0.995      0.919      0.995      0.919      0.918      1.000
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
    Gene.Name log2FoldChange padj
42    NDUFAF7      -1.444831    0
106     ITGA3       1.832659    0
116      YBX2      -1.804203    0
128    MAP3K9       1.277301    0
131     KDM7A       2.077349    0
157     PROM1      -3.585410    0
168     SCN4A      -3.583746    0
175  CACNA2D2      -1.741931    0
180     TEAD3       2.622174    0
188    JARID2      -1.039358    0
223      PAX7       4.837460    0
240       CD9      -1.985595    0
265      MRC2       2.454198    0
328    MAMLD1       5.297447    0
358   SLC38A5      -3.322193    0
367    ATP1A2      -3.921468    0
433       VIM       1.370777    0
468  ARHGAP31       3.751753    0
483      GAB2       1.767040    0
490    TMSB10       1.519033    0
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
      Gene.Name log2FoldChange padj  baseMean
791       LZTS1       4.684407    0 32976.732
6841      GATA4       4.826014    0  6038.701
10235     EOMES       5.791150    0 54252.276
10593     SOX17       5.456209    0 15988.506
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

