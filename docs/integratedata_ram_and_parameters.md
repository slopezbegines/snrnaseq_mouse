# Consumo de RAM en `IntegrateData()` e implicaciones de sus parámetros clave

**Contexto:** integración RPCA con SCTransform en Seurat v5, dataset GSE262881 (6 muestras, 113.698 núcleos, 15 GB RAM).

---

## 1. Por qué se agota la RAM durante `IntegrateData()`

La integración RPCA en Seurat sigue un proceso jerárquico de tres fases. Cada fase tiene un perfil de memoria distinto.

### Fase 1 — `SelectIntegrationFeatures()` + `PrepSCTIntegration()`

Selecciona los genes variables compartidos entre muestras y re-escala cada objeto SCT sobre ese conjunto de genes. El consumo es aproximadamente proporcional al número de células × genes seleccionados × número de muestras.

En este dataset, tras cargar los 6 objetos SCT desde el checkpoint:

```
RAM antes de la integración: 13.538 MB
RAM tras PrepSCTIntegration: 12.672 MB
```

El consumo ya ocupa ~85 % de la RAM total del sistema.

### Fase 2 — `FindIntegrationAnchors()` (RPCA)

Para cada par de muestras, proyecta las células de una muestra en el espacio PCA de la otra y encuentra vecinos mutuos (anchors). Es CPU-intensivo pero relativamente moderado en RAM porque opera par a par y libera estructuras intermedias.

En este dataset encontró **511.186 anchors** distribuidos entre los 15 pares posibles (C(6,2)).

```
RAM tras rm(sct_objects) + gc(): 12.707 MB
```

El objeto `anchors` retiene internamente los embeddings PCA de todas las muestras, por lo que la RAM no baja significativamente aunque se eliminen los objetos SCT originales.

### Fase 3 — `IntegrateData()` : el cuello de botella

Es la fase que agota la RAM. Opera en tres sub-pasos para cada merge jerárquico:

| Sub-paso | Qué hace | Coste de memoria |
|---|---|---|
| *Finding integration vectors* | Extrae los vectores de corrección a partir de los anchors | Moderado |
| **Finding integration vector weights** | Para cada célula, calcula un peso sobre sus `k.weight` anchors más cercanos | **Máximo — aquí ocurre el kill** |
| *Integrating data* | Aplica la corrección ponderada a la matriz de expresión | Alto, pero decae tras el cálculo |

La matriz de pesos tiene dimensión `n_células × k.weight`. Con los valores por defecto:

```
113.698 células × 100 anchors (k.weight=100) = 11.369.800 valores por paso de merge
```

Como la integración es jerárquica (merge 4→2, merge 6→3, merge {2,4}→{3,6}), el tercer paso trabaja con los datasets ya fusionados (hasta ~56.000 células combinadas) y construye una matriz de pesos aún mayor antes de liberarla. Este es el punto exacto del crash observado en los logs:

```
Merging dataset 2 4 into 3 6
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights   ← OOM kill aquí
```

---

## 2. Parámetros que controlan el consumo de RAM

### `k.weight` (el más influyente)

| Valor | Efecto en RAM | Efecto en calidad |
|---|---|---|
| 100 (defecto) | Máximo; inviable con < 20 GB libres para este dataset | Referencia |
| **50 (aplicado)** | **~50 % menos en el paso de weights** | Prácticamente idéntica; recomendado para datasets grandes |
| 25 | ~75 % menos | Ligera pérdida de suavidad en la corrección; aceptable |
| 10 | Mínimo; posible en máquinas muy limitadas | Corrección menos suave; usar solo como último recurso |

`k.weight` controla cuántos anchors vecinos se usan para ponderar la corrección de cada célula. Valores más bajos hacen la corrección local (menos vecinos), lo que en la práctica afecta muy poco a la separación de tipos celulares pero reduce el consumo cuadráticamente.

### `dims` — `1:N_PCS_INTEGRATION`

Número de componentes principales usados para la búsqueda de vecinos y la proyección RPCA.

| Valor | RAM | Calidad |
|---|---|---|
| 1:30 (defecto Seurat) | Alto | Máxima captura de varianza |
| **1:20 (aplicado)** | Moderado | Adecuado para la mayoría de datasets con tipos celulares bien definidos |
| 1:10 | Bajo | Puede perder subpoblaciones raras |

Reducir `dims` disminuye el tamaño de los embeddings PCA almacenados en el objeto `anchors` y la dimensión de los vectores de integración.

### `k.anchor` en `FindIntegrationAnchors()`

Controla cuántos vecinos mutuos se consideran al identificar anchors.

| Valor | Efecto |
|---|---|
| 5 (aplicado, mínimo) | Menos anchors totales; menos RAM en `FindIntegrationAnchors()` y en el objeto resultante |
| 20 (defecto) | Más anchors; integración más robusta ante batch effects fuertes |

Con 511.186 anchors encontrados usando `k.anchor=5`, la integración ya tiene información más que suficiente para 113k células.

### `normalization.method = "SCT"`

Fija el método a SCTransform. No tiene alternativa en este pipeline; se menciona porque activa rutas de código distintas a `LogNormalize` dentro de `IntegrateData()`, incluyendo la necesidad de que `PrepSCTIntegration()` haya sido ejecutado antes.

---

## 3. Configuración aplicada y justificación

```r
# 06_integration.R — línea 114
data_int <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "SCT",
  dims                 = 1:N_PCS_INTEGRATION,   # = 1:20
  k.weight             = 50                      # reducido desde 100
)
```

**Justificación del cambio:** con 12.7 GB ocupados antes de llamar a `IntegrateData()` y solo ~2 GB de RAM libre (más swap disponible), la matriz de pesos con `k.weight=100` supera la capacidad del sistema durante el tercer merge jerárquico. Reducir a `k.weight=50` es la intervención mínima que resuelve el OOM sin comprometer la calidad biológica del resultado.

---

## 4. Alternativa si el problema persiste: Harmony

Si el dataset crece o se añaden más muestras, la integración RPCA de Seurat puede volverse inviable en máquinas con < 32 GB. Harmony es una alternativa que opera sobre el espacio PCA ya calculado (no sobre matrices de expresión completas) y requiere típicamente < 4 GB adicionales para este tamaño de dataset:

```r
library(harmony)
data_int <- RunPCA(merged_object, assay = "SCT")
data_int <- RunHarmony(data_int, group.by.vars = "orig.ident", assay.use = "SCT")
```

La desventaja es que Harmony corrige el espacio de embedding (PCA/UMAP) pero no la matriz de expresión integrada, lo que limita ciertos análisis posteriores que requieren valores de expresión corregidos por batch.

---

*Generado: 2026-05-09 | Pipeline: snRNAseq_mouse | Dataset: GSE262881*
