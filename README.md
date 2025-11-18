ANÁLISIS DE EXPRESIÓN DIFERENCIAL RNA-SEQ: 0 horas vs 72 horas
Basado en el estudio: Diferenciación de células madre embrionarias a endodermo
Experimento: E-MTAB-9194 (Expression Atlas)
Comparación: embryonic stem cell (0h) vs definitive endoderm cell (72h)
Criterios: FC ≥ 1.5, p-adj < 0.01 (según metodología del paper original)

0) Instalación y carga de paquetes necesarios

library(BiocManager)
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(EnhancedVolcano)
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
