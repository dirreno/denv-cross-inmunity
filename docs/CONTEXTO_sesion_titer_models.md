# Contexto de sesión — modelos de títulos antigénicos (augur titers)

Documento de traspaso. Resume qué se hizo, por qué, qué se descubrió y qué queda
pendiente. Escrito para que otra sesión retome el trabajo sin repetir el camino.

---

## 1. Pregunta del proyecto

Se dispone de secuencias colombianas de DENV sin datos serológicos, y de las mediciones
de neutralización de Katzelnick et al. (2015), que sí tienen títulos pero para otros
virus. El objetivo es **estimar títulos para los aislados colombianos**.

La duda inicial del usuario fue si `augur titers` sirve para eso. La respuesta corta,
tras investigarlo a fondo: **sí predice fuera de muestra, pero con una limitación
estructural que resultó determinante en este caso.**

---

## 2. Decisiones conceptuales (y por qué)

### Qué predice realmente `augur titers`

El modelo (Neher et al., 2016) es:

```
log2(título) = potencia_suero − avidez_virus − Σ(caídas)
```

La validación cruzada del paper original oculta **celdas** de la matriz de títulos,
conservando el virus y el suero en el entrenamiento: es *matrix completion*. El caso de
este proyecto es distinto — una **fila entera nueva**, un virus nunca medido. Su avidez
no es identificable.

**Mitigante importante encontrado en el código:** la avidez es el término más
regularizado del modelo (`lam_avi=3.0` frente a `lam_drop=1.0` y `lam_pot=0.5`). Se
contrae agresivamente hacia cero ya en el entrenamiento, así que fijarla en 0 para un
virus nuevo no es un salto grande. Esto debilitó una objeción inicial.

### tree vs sub

- `tree` localiza el cambio antigénico en la **filogenia** (`dTiter`/`cTiter` por nodo).
  Solo transfiere señal a linajes atravesados por pares medidos; las ramas sin datos
  reciben 0.
- `sub` lo localiza en la **secuencia** (`dTiterSub`/`cTiterSub`, más una tabla de
  efectos por sustitución). Transfiere a cualquier secuencia portadora, medida o no.

**`sub` es el motor para este problema; `tree` es el contraste.** Correr siempre **por
serotipo**, nunca conjunto: el modelo asume un escenario tipo influenza (un serotipo,
distancia aditiva sobre el árbol).

### Precedente

`blab/dengue-antigenic-dynamics` (Bell et al., 2019) es exactamente este modelo aplicado
a dengue con los datos de Katzelnick. Ese repo es **Python 2** y su propio README
recomienda usar la implementación modular de `augur`. Sirve como referencia de
configuración, no como código a correr.

---

## 3. Hechos verificados sobre el código de augur

Comprobados leyendo la fuente instalada (`augur 33.0.0`), no de memoria. **No volver a
cuestionarlos sin releer la fuente.**

| Hecho | Detalle |
|---|---|
| Formato de `--titers` | TSV de **5 columnas sin encabezado**: virus, virus_referencia, id_suero, id_fuente, título |
| Censura | La maneja augur: `<10` → 5.0 (umbral/2), `>5120` → 10240. **No pre-convertir.** `*` se descarta silenciosamente |
| Salida `tree` | `titers`, `potency`, `avidity`, `nodes{name: dTiter, cTiter}` |
| Salida `sub` | `titers`, `potency`, `avidity`, `substitution`, y `nodes{name: dTiterSub, cTiterSub}` solo si se pasa `--tree` |
| `titers` en el JSON | **Son MEDICIONES, no predicciones**: `[valor_normalizado_observado, mediana_cruda_observada]`. Confirmado en `compile_titers`, que itera `self.titers.titers_normalized` |
| `potency` | `ref → {serum_id: valor, 'mean_potency': valor}`. **Excluir la clave `mean_potency`**, no es un suero |
| `avidity` | Plano, `virus → efecto`. Solo virus medidos |
| `substitution` | Plano, `"E:A123B" → efecto`, posición **1-based**. Solo efectos `> 1e-4`, redondeados a 4 decimales |
| Virus sin mediciones | `predict_titer` **sí funciona** si está en el árbol/alineamiento; `avidity` cae a 0.0 por defecto. Devuelve `None` solo si el virus no está |
| Ajuste | `nnl1reg`: mínimos cuadrados **no negativos** con L1. El modelo no puede representar una sustitución que *suba* el título |
| Validación por virus | `make_training_set(subset_strains=True)` parte por virus de prueba en vez de por celda |
| El CLI no emite títulos predichos | Nunca llama `predict_titer` para pares no medidos |

---

## 4. Qué se implementó en `notebooks/AugurPipeline.ipynb`

Celdas añadidas, en orden de ejecución:

1. **Rutas de títulos** — añadidas a la celda *Merged Paths*
   (`TITERS_MERGE`, `SERA_MAP_MERGE`, `TITER_TREE_MODEL_MERGE`, `TITER_SUB_MODEL_MERGE`, `RAW_TITERS`).
2. **`augur_titers_tree()` / `augur_titers_sub()`** — wrappers (los escribió el usuario a
   partir de la doc oficial). Ojo: `--titers` y `--alignment` son nargs múltiples;
   `--gene-names` debe ir en el mismo orden que `--alignment`.
3. **Regenerar GFF temporal de merged** — ver sección 5.
4. **`augur translate` para merged** — no existía; `sub` lo necesita.
5. **Preparación del insumo de títulos** — conversión + resolución de nombres de suero.
6. **Ejecución de `tree` y `sub`** por serotipo.
7. **Análisis exploratorio** — 5 bloques + figura (`figures/titer_model_eda.png`).
8. **Títulos predichos para colombianas** — reconstrucción de `predict_titer`
   (→ `output/predicted_titers_colombia_sub.csv`).

---

## 5. Problemas encontrados y resueltos

### GFF de merged con múltiples seqids (bloqueaba `augur translate`)

La celda de merge llamaba `merge_gffs(TMP_GFF_COL, TMP_GFF_EXP, TMP_GFF_MERGE)`,
concatenando dos GFF con seqids distintos. `augur translate` exige **uno solo**.

Fusionar GFFs no tiene sentido: el GFF solo declara las coordenadas del gen E
**relativas al alineamiento**, idénticas para todas las secuencias ya alineadas. El
seqid es una etiqueta.

**Resuelto:** se comentó la llamada a `merge_gffs` y se añadió una celda que regenera el
GFF desde el alineamiento merged, con `assert` de un único seqid.
**DENV3 mide 1479 nt, no 1485** — el largo se toma por serotipo, no fijo.

### Resolución de nombres de suero

Los nombres de suero (`DEN1_BURMA_2005_61117`) no coinciden con las puntas
(`DENV1/MYANMAR/2005/61117`). **No emparejar por país** (Burma/Myanmar) **ni confiar en
el año** (el del suero puede ser el de extracción, no el del virus: `DEN2_SENEGAL_2003_SENDAK_HD_0674`
→ `DENV2/SENEGAL/1970/SENDAK-HD-10674`).

Criterio adoptado: serotipo + año + accession normalizado, aceptando solo coincidencias
unívocas. Resultado: 24 automáticos + 5 manuales (`SERA_MANUAL`) = **29 de 36**.

Los **7 excluidos** (`SERA_ABSENT`) lo están porque su virus de referencia no existe en
`full-genome-only-sequences.fasta`. Verificado — no es el filtro de recorte.
Recuperarlos exige bajar secuencias de GenBank.

### Bug preexistente, NO corregido

En las celdas de `augur_translate` de Colombia y experimental se pasa
`ref_seq=TMP_GFF_COL.format(serotype)` pero el bucle itera sobre `sero`. `serotype` es
una variable colgada del bucle anterior y vale 4 en las cuatro iteraciones: **los cuatro
serotipos se traducen contra el GFF de DENV4.** La celda merged añadida usa `sero`
correctamente.

---

## 6. Resultados

### Cobertura del insumo

| | mediciones | sueros al insumo | sueros retenidos | censuradas |
|---|---|---|---|---|
| DENV1 | 534 | 8 | 6 | 1,7 % |
| DENV2 | 348 | 9 | 7 | 12,9 % |
| DENV3 | 418 | 6 | 4 | 1,4 % |
| DENV4 | 318 | 6 | 4 | 28,3 % |

Virus: 404/412 coinciden con puntas. Al restringir a intraserotípico se pierden **todas**
las mediciones entre serotipos — estos modelos no dicen nada de reactividad cruzada.

### Ajuste

RMSE sobre pares medidos: **0,48 / 0,64 / 0,53 / 0,74** log2. Para DENV1 coincide con la
métrica que augur imprime al entrenar (`rms deviation on training set = 0.47599`). Está
**por debajo del ruido del ensayo** (~1 log2): el modelo ajusta bien.

### Hallazgo principal (negativo)

El modelo seleccionó un vocabulario diminuto:

| | subs. con efecto | presentes en col. | secuencias col. | proteínas E únicas | **valores predichos únicos** | predicciones nulas |
|---|---|---|---|---|---|---|
| DENV1 | 12 | 7 | 304 | 73 | **4** | 17 % |
| DENV2 | 10 | 7 | 474 | 141 | **6** | 23 % |
| DENV3 | 10 | 5 | 263 | 64 | **3** | 25 % |
| DENV4 | 6 | **0** | 46 | 17 | **1** | **100 %** |

**El modelo `sub` entrenado sobre el panel de Katzelnick no resuelve la variación
antigénica intraserotípica de Colombia.** En DENV4 ninguna sustitución con efecto está
presente en las secuencias colombianas: todas las predicciones son el modelo nulo.

Causa: dispersión inducida por la regularización L1 sobre un conjunto serológico
pequeño. No es un defecto del ajuste. Para DENV1, augur además colapsó la matriz de
diseño de 158 a 141 columnas agrupando sustituciones colineales.

**Implicación estratégica:** `sub` no puede corroborar la proyección cartográfica que el
usuario ya tiene (`notebooks/AntigenicProjection.ipynb`, `output/projected_coordinates_*.csv`),
porque carece de poder de resolución aquí. Esa proyección usa la distancia genética
completa, no solo las sustituciones con efecto ajustado, y por eso conserva resolución.
**Queda como instrumento principal.**

### Incertidumbre

SD de las avideces ajustadas: 0,39 / 0,40 / 0,34 / 0,41 log2 (≈ 1,3× en título). Es la
cota inferior del error de cada predicción colombiana, donde la avidez se fija en 0.

---

## 7. Verificación hecha

- **Reconstrucción de `predict_titer`:** se reajustó el modelo con la API de Python
  (`SubstitutionModel` + `prepare()` + `train()`) y se comparó par por par contra la
  reconstrucción desde el JSON. DENV1, 524 pares: **máx 1,235e-4**, RMSE 4,5e-5. Esa
  discrepancia es exactamente el umbral de truncamiento con que augur escribe los
  efectos. **La reconstrucción es exacta salvo redondeo.**
- **Anclaje de escala:** SD residual ~1e-5 log2 → la relación
  `log2(crudo) = ancla − normalizado` se sostiene; `pred_titer` es interpretable.
- **Celdas probadas antes de entregar**, con datos sintéticos en scratchpad: la de
  preparación, la exploratoria y la de predicción (esta última con round-trip exacto,
  RMSE 6e-16).

---

## 8. Pendiente

1. **Validación leave-one-VIRUS-out.** Es el nivel que corresponde a la pregunta
   científica y no se ha hecho. Usar `make_training_set(subset_strains=True)`. Reportar
   contra dos referencias: el modelo nulo (solo potencia) y el oráculo (coordenadas
   verdaderas). La CV por celdas de Neher et al. (2016) **no** aplica a este caso.
2. **Celda de comparación contra `predict_titer`** para los cuatro serotipos (hecha solo
   para DENV1, fuera del notebook).
3. **Modelo `tree` como contraste** — el bloque 4 del análisis exploratorio calcula el
   Spearman `cTiter` vs `cTiterSub`, pero no se ha discutido el resultado.
4. **Decidir si se recuperan los 7 sueros excluidos** bajando sus secuencias de GenBank.
5. **Cruzar las sustituciones con efecto contra el alineamiento colombiano** para
   identificar si alguna con peso es exclusiva de linajes colombianos.

---

## 9. Advertencias de método acordadas

- No aflojar la regularización para fabricar resolución que los datos no sostienen.
- Los sueros de Katzelnick son de **infección primaria** (monotípicos). Si el objetivo
  final es inmunidad cruzada poblacional en Colombia —donde domina la infección
  secundaria— hay una limitación de validez externa que debe declararse.
- El piso de ruido del ensayo (~1 log2) acota cualquier afirmación sobre variación
  intraserotípica.
- `cTiter` / `cTiterSub` son **posición antigénica relativa**, no títulos. No
  presentarlos como títulos.

---

## 10. Archivos relevantes

| Ruta | Contenido |
|---|---|
| `notebooks/AugurPipeline.ipynb` | Pipeline completo; celdas de títulos al final |
| `notebooks/AntigenicProjection.ipynb` | Proyección cartográfica (ruta principal) |
| `resources/experimental_data/neutralization-titers.txt` | Mediciones de Katzelnick |
| `resources/experimental_data/2d-antigenic-map.txt` | Mapa con coordenadas AG **y SR** |
| `resources/merged_data/titers/` | TSV de insumo, mapa de sueros, JSON de modelos |
| `output/predicted_titers_colombia_sub.csv` | Títulos predichos |
| `docs/metodos_modelos_titulos.md` | Redacción metodológica con referencias |

---

## 11. Referencias

Bell, S. M., Katzelnick, L., & Bedford, T. (2019). Dengue genetic divergence generates
within-serotype antigenic variation, but serotypes dominate evolutionary dynamics.
*eLife*, *8*, e42496. https://doi.org/10.7554/eLife.42496

Huddleston, J., Hadfield, J., Sibley, T. R., Lee, J., Fay, K., Ilcisin, M., Harkins, E.,
Bedford, T., Neher, R. A., & Hodcroft, E. B. (2021). Augur: A bioinformatics toolkit for
phylogenetic analyses of human pathogens. *Journal of Open Source Software*, *6*(57),
2906. https://doi.org/10.21105/joss.02906

Katzelnick, L. C., Fonville, J. M., Gromowski, G. D., Bustos Arriaga, J., Green, A.,
James, S. L., Lau, L., Montoya, M., Wang, C., VanBlargan, L. A., Russell, C. A., Thu,
H. M., Pierson, T. C., Buchy, P., Aaskov, J. G., Muñoz-Jordán, J. L., Vasilakis, N.,
Gibbons, R. V., Tesh, R. B., … Harris, E. (2015). Dengue viruses cluster antigenically
but not as discrete serotypes. *Science*, *349*(6254), 1338–1343.
https://doi.org/10.1126/science.aac5017

Neher, R. A., Bedford, T., Daniels, R. S., Russell, C. A., & Shraiman, B. I. (2016).
Prediction, dynamics, and visualization of antigenic phenotypes of seasonal influenza
viruses. *Proceedings of the National Academy of Sciences*, *113*(12), E1701–E1709.
https://doi.org/10.1073/pnas.1525578113

Smith, D. J., Lapedes, A. S., de Jong, J. C., Bestebroer, T. M., Rimmelzwaan, G. F.,
Osterhaus, A. D. M. E., & Fouchier, R. A. M. (2004). Mapping the antigenic and genetic
evolution of influenza virus. *Science*, *305*(5682), 371–376.
https://doi.org/10.1126/science.1097211
