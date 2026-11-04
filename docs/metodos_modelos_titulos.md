# Modelos de títulos antigénicos: métodos

## Datos serológicos y preparación del insumo

Se partió de las mediciones de neutralización publicadas por Katzelnick et al. (2015),
que comprenden 412 virus de prueba y 36 sueros de infección primaria. Estas mediciones
se reformatearon al esquema de entrada de `augur titers` (Huddleston et al., 2021):
un archivo delimitado por tabulaciones y sin encabezado, con las columnas virus de
prueba, virus de referencia del suero, identificador de suero, identificador de fuente
y título.

La correspondencia entre los nombres de cepa del archivo serológico y las puntas del
árbol filogenético requirió resolución explícita, ya que ambos conjuntos emplean
convenciones de nomenclatura distintas. Los nombres de virus de prueba coincidieron
directamente en 404 de 412 casos. Los nombres de suero, en cambio, no admiten
emparejamiento por país —el archivo serológico emplea topónimos históricos (por
ejemplo, *Burma* donde el árbol registra *Myanmar*)— ni por año, dado que el año
consignado en el identificador de suero corresponde en algunos casos a la extracción
del suero y no al aislamiento del virus. El emparejamiento se realizó, por tanto, sobre
serotipo y número de acceso normalizado, aceptando únicamente coincidencias unívocas.
Este procedimiento resolvió 24 de los 36 sueros; otros 5 se asignaron manualmente tras
inspección individual. Los 7 sueros restantes se excluyeron porque su virus de
referencia no está presente en el conjunto de secuencias de genoma completo disponible
y, en consecuencia, no puede figurar como punta del árbol.

Dado que los árboles se construyeron por serotipo, se conservaron únicamente los pares
intraserotípicos, en los que tanto el virus de prueba como el virus de referencia del
suero son puntas del mismo árbol. Esta restricción excluye la totalidad de las
mediciones interserotípicas y, por consiguiente, los modelos aquí descritos no informan
sobre reactividad cruzada entre serotipos. Tiene como contrapartida que la censura a la
izquierda, predominante en las mediciones interserotípicas, resulta marginal en el
conjunto retenido (1,7 %, 12,9 %, 1,4 % y 28,3 % de los títulos para DENV1–DENV4,
respectivamente). Los valores censurados se transmitieron sin transformación previa,
delegando en `augur` su convención interna, que asigna a un título reportado como `<X`
el valor X/2. Las celdas sin medición se descartaron.

El conjunto final comprendió 534, 348, 418 y 318 mediciones para DENV1–DENV4,
correspondientes a 8, 9, 6 y 6 sueros.

## Ajuste de los modelos

Sobre cada serotipo se ajustaron de forma independiente los dos modelos de títulos
implementados en `augur` (versión 33.0.0), siguiendo la formulación de Neher et al.
(2016). Ambos descomponen el título normalizado en escala log2 como

    log2(título) = potencia del suero − avidez del virus − Σ(caídas de título)

y difieren en la atribución del término de caídas. El **modelo de árbol** (`augur
titers tree`) asigna una caída a cada rama de la filogenia y anota en cada nodo la
caída propia (`dTiter`) y la acumulada desde la raíz (`cTiter`). El **modelo de
sustituciones** (`augur titers sub`) atribuye las caídas a sustituciones individuales de
aminoácido, anotando los atributos homólogos `dTiterSub` y `cTiterSub`. El primero
localiza el cambio antigénico en la filogenia; el segundo, en la secuencia. La
diferencia es consecuente para este trabajo: el modelo de sustituciones transfiere el
efecto estimado a cualquier secuencia portadora, haya sido o no medida, mientras que el
modelo de árbol solo transfiere señal a los linajes atravesados por pares medidos.

El modelo de árbol se ajustó sobre los árboles datados obtenidos con `augur refine`; el
de sustituciones, sobre los alineamientos de aminoácidos de la proteína de envoltura (E)
generados por `augur translate`. En ambos casos la estimación emplea mínimos cuadrados
no negativos con regularización L1 (`nnl1reg`), con los coeficientes de penalización
predeterminados del paquete: 1,0 para las caídas, 0,5 para la potencia del suero y 3,0
para la avidez del virus. La restricción de no negatividad implica que el modelo no
puede representar sustituciones que incrementen el título.

De los 29 sueros incorporados al insumo, el ajuste retuvo 21 (6, 7, 4 y 4 para
DENV1–DENV4). Para DENV1, el procedimiento colapsó la matriz de diseño de 158 a 141
columnas al agrupar sustituciones colineales de co-ocurrencia invariable, que los datos
no permiten separar.

## Reconstrucción de títulos predichos

La interfaz de línea de comandos de `augur titers` no emite títulos predichos para pares
virus–suero no medidos: su salida comprende las mediciones compiladas, los parámetros
ajustados (potencias, avideces y efectos de sustitución) y las anotaciones por nodo
(`cTiter`, `cTiterSub`), que constituyen una posición antigénica relativa y no un
título. Para obtener títulos predichos para los aislados colombianos se reimplementó la
función `predict_titer` del modelo de sustituciones a partir de los parámetros
almacenados, sumando el efecto estimado de cada sustitución de aminoácido que separa la
cepa de referencia de cada suero de la secuencia colombiana evaluada.

La avidez es un parámetro específico de virus que no puede estimarse para una cepa sin
mediciones propias; siguiendo el comportamiento de `augur`, se fijó en cero. La
desviación estándar de las avideces ajustadas sobre las cepas medidas se reporta como
cota inferior de la incertidumbre de cada predicción, y resultó de 0,39, 0,40, 0,34 y
0,41 unidades log2 para DENV1–DENV4, equivalentes a un factor aproximado de 1,3 en el
título.

La escala se ancló empíricamente en lugar de asumir la forma de la normalización
interna. La salida de `augur` conserva, para cada par medido, el título normalizado y la
mediana del título crudo; sobre esos pares se estimó un ancla por suero y se verificó su
consistencia interna, cuya desviación residual fue del orden de 10⁻⁵ unidades log2.

## Validación

Se distinguieron tres niveles de verificación.

**Fidelidad de la reconstrucción.** Para confirmar que la reimplementación reproduce la
función de `augur`, se reajustó el modelo de sustituciones mediante la interfaz de
Python del paquete y se compararon sus predicciones con las reconstruidas, par por par.
Sobre los 524 pares medidos de DENV1 la discrepancia máxima absoluta fue de 1,2 × 10⁻⁴
unidades log2, valor que corresponde al umbral de truncamiento con que `augur` escribe
los efectos de sustitución en su salida. La reconstrucción es, por tanto, exacta salvo
redondeo.

**Ajuste sobre las mediciones.** El error cuadrático medio de las predicciones del
modelo frente a los títulos medidos fue de 0,48, 0,64, 0,53 y 0,74 unidades log2 para
DENV1–DENV4. Para DENV1 este valor coincide con la métrica que el propio paquete reporta
durante el entrenamiento. Dado que la reproducibilidad de los ensayos de neutralización
se sitúa en torno a 1 unidad log2, estos residuos se encuentran por debajo del umbral de
ruido experimental. Cabe subrayar que un ajuste perfecto no constituye el objetivo del
procedimiento: la regularización L1 contrae deliberadamente los coeficientes en favor de
la capacidad de generalización.

**Capacidad predictiva fuera de muestra.** La validación cruzada de Neher et al. (2016)
omite celdas individuales de la matriz de títulos conservando el virus y el suero en el
conjunto de entrenamiento, de modo que los términos de potencia y avidez siguen siendo
estimables. Ese esquema no corresponde a la pregunta aquí planteada, que requiere
predecir la totalidad de la fila de un virus nunca medido. El esquema pertinente es la
omisión por virus, disponible en el paquete mediante el parámetro `subset_strains`.
*(Nota para los autores: este tercer nivel no se ha ejecutado todavía; no debe
reportarse como realizado.)*

## Resultados del ajuste relevantes para la interpretación

El modelo de sustituciones seleccionó un conjunto reducido de posiciones: 12, 10, 10 y 6
sustituciones con efecto no nulo para DENV1–DENV4. De ellas, 7, 7, 5 y 0 segregan entre
las secuencias colombianas analizadas.

En consecuencia, la resolución del modelo sobre el conjunto colombiano es limitada. Las
304, 474, 263 y 46 secuencias evaluadas por serotipo —que comprenden 73, 141, 64 y 17
proteínas E distintas— reciben únicamente 4, 6, 3 y 1 valores de título predicho
distintos. Para DENV4, ninguna de las sustituciones con efecto estimado está presente en
las secuencias colombianas, de modo que la totalidad de sus predicciones se reduce al
modelo nulo, es decir, a la potencia del suero. La proporción de predicciones sin
contribución de sustituciones fue del 17 %, 23 %, 25 % y 100 % para DENV1–DENV4.

Este resultado es atribuible a la dispersión inducida por la regularización L1 sobre un
conjunto serológico de tamaño reducido, y no a un defecto del ajuste: el modelo describe
adecuadamente las mediciones disponibles, pero el vocabulario antigénico que identifica
no cubre la variación presente en los linajes colombianos.

## Programas empleados

Los análisis se realizaron con `augur` 33.0.0 (Huddleston et al., 2021), componente del
ecosistema Nextstrain (Hadfield et al., 2018), sobre alineamientos generados con MAFFT
7.526 y filogenias inferidas con IQ-TREE 3.0.1. El procesamiento se efectuó en Python
con Biopython 1.86, NumPy 2.4.2, SciPy 1.17.0 y pandas 2.3.3.

## Referencias

Bell, S. M., Katzelnick, L., & Bedford, T. (2019). Dengue genetic divergence generates
within-serotype antigenic variation, but serotypes dominate evolutionary dynamics.
*eLife*, *8*, e42496. https://doi.org/10.7554/eLife.42496

Hadfield, J., Megill, C., Bell, S. M., Huddleston, J., Potter, B., Callender, C.,
Sagulenko, P., Bedford, T., & Neher, R. A. (2018). Nextstrain: Real-time tracking of
pathogen evolution. *Bioinformatics*, *34*(23), 4121–4123.
https://doi.org/10.1093/bioinformatics/bty407

Huddleston, J., Hadfield, J., Sibley, T. R., Lee, J., Fay, K., Ilcisin, M., Harkins, E.,
Bedford, T., Neher, R. A., & Hodcroft, E. B. (2021). Augur: A bioinformatics toolkit for
phylogenetic analyses of human pathogens. *Journal of Open Source Software*, *6*(57),
2906. https://doi.org/10.21105/joss.02906

Katzelnick, L. C., Coello Escoto, A., Huang, A. T., Garcia-Carreras, B., Chowdhury, N.,
Maljkovic Berry, I., Chavez, C., Buchy, P., Duong, V., Dussart, P., Gromowski, G.,
Macareo, L., Thaisomboonsuk, B., Fernandez, S., Smith, D. J., Jarman, R., Whitehead,
S. S., Salje, H., & Cummings, D. A. T. (2021). Antigenic evolution of dengue viruses
over 20 years. *Science*, *374*(6570), 999–1004. https://doi.org/10.1126/science.abk0058

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
