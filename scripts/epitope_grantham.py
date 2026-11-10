# -*- coding: utf-8 -*-
"""Distancia antigénica p-all-epitope de Grantham sobre el gen E.

Las definiciones —matriz de Grantham, sitios epitópicos, numeración común y
normalización— son las mismas que usa `notebooks/DistancesCalculation.ipynb`
para construir `output/distances_epitope_grantham_DENV<s>.csv`. Están aquí
para poder calcular además matrices **dentro de un conjunto** (colombianas
contra colombianas), que aquel notebook no produce: solo calcula colombianas
contra las referencias del mapa experimental.

`antigenic_distance` es la implementación de referencia, par a par. `matriz`
hace lo mismo vectorizado sobre un conjunto entero, y `autoprueba()` comprueba
que las dos coinciden; si se toca una, hay que volver a correr la prueba:

    python scripts/epitope_grantham.py <ruta.fasta>
"""
import sys

import numpy as np

E_LEN_STD = 495   # longitud del gen E en DENV1, DENV2 y DENV4
DENV3_GAP = 157   # indice 0-based donde DENV3 pierde los residuos 158-159

AA = "ARNDCQEGHILKMFPSTWYV"

GRANTHAM = {
    'A': {'A':0,'R':112,'N':111,'D':126,'C':195,'Q':91,'E':107,'G':60,'H':86,'I':94,'L':96,'K':106,'M':84,'F':113,'P':27,'S':99,'T':58,'W':148,'Y':112,'V':64},
    'R': {'A':112,'R':0,'N':86,'D':96,'C':180,'Q':43,'E':54,'G':125,'H':29,'I':97,'L':102,'K':26,'M':91,'F':97,'P':103,'S':110,'T':71,'W':101,'Y':77,'V':96},
    'N': {'A':111,'R':86,'N':0,'D':23,'C':139,'Q':46,'E':42,'G':80,'H':68,'I':149,'L':153,'K':94,'M':142,'F':158,'P':91,'S':46,'T':65,'W':174,'Y':143,'V':133},
    'D': {'A':126,'R':96,'N':23,'D':0,'C':154,'Q':61,'E':45,'G':94,'H':81,'I':168,'L':172,'K':101,'M':160,'F':177,'P':108,'S':65,'T':85,'W':181,'Y':160,'V':152},
    'C': {'A':195,'R':180,'N':139,'D':154,'C':0,'Q':154,'E':170,'G':159,'H':174,'I':198,'L':198,'K':202,'M':196,'F':205,'P':169,'S':112,'T':149,'W':215,'Y':194,'V':192},
    'Q': {'A':91,'R':43,'N':46,'D':61,'C':154,'Q':0,'E':29,'G':87,'H':24,'I':109,'L':113,'K':53,'M':101,'F':116,'P':76,'S':68,'T':42,'W':130,'Y':99,'V':96},
    'E': {'A':107,'R':54,'N':42,'D':45,'C':170,'Q':29,'E':0,'G':98,'H':40,'I':134,'L':138,'K':56,'M':126,'F':140,'P':93,'S':80,'T':65,'W':152,'Y':122,'V':121},
    'G': {'A':60,'R':125,'N':80,'D':94,'C':159,'Q':87,'E':98,'G':0,'H':98,'I':135,'L':138,'K':127,'M':127,'F':153,'P':42,'S':56,'T':59,'W':184,'Y':147,'V':109},
    'H': {'A':86,'R':29,'N':68,'D':81,'C':174,'Q':24,'E':40,'G':98,'H':0,'I':94,'L':99,'K':32,'M':87,'F':100,'P':77,'S':89,'T':47,'W':115,'Y':83,'V':84},
    'I': {'A':94,'R':97,'N':149,'D':168,'C':198,'Q':109,'E':134,'G':135,'H':94,'I':0,'L':5,'K':102,'M':10,'F':21,'P':95,'S':142,'T':89,'W':61,'Y':33,'V':29},
    'L': {'A':96,'R':102,'N':153,'D':172,'C':198,'Q':113,'E':138,'G':138,'H':99,'I':5,'L':0,'K':107,'M':15,'F':22,'P':98,'S':145,'T':92,'W':61,'Y':36,'V':32},
    'K': {'A':106,'R':26,'N':94,'D':101,'C':202,'Q':53,'E':56,'G':127,'H':32,'I':102,'L':107,'K':0,'M':95,'F':102,'P':103,'S':121,'T':78,'W':110,'Y':85,'V':97},
    'M': {'A':84,'R':91,'N':142,'D':160,'C':196,'Q':101,'E':126,'G':127,'H':87,'I':10,'L':15,'K':95,'M':0,'F':28,'P':87,'S':135,'T':81,'W':67,'Y':36,'V':21},
    'F': {'A':113,'R':97,'N':158,'D':177,'C':205,'Q':116,'E':140,'G':153,'H':100,'I':21,'L':22,'K':102,'M':28,'F':0,'P':114,'S':155,'T':103,'W':40,'Y':22,'V':50},
    'P': {'A':27,'R':103,'N':91,'D':108,'C':169,'Q':76,'E':93,'G':42,'H':77,'I':95,'L':98,'K':103,'M':87,'F':114,'P':0,'S':74,'T':38,'W':147,'Y':110,'V':68},
    'S': {'A':99,'R':110,'N':46,'D':65,'C':112,'Q':68,'E':80,'G':56,'H':89,'I':142,'L':145,'K':121,'M':135,'F':155,'P':74,'S':0,'T':58,'W':177,'Y':144,'V':124},
    'T': {'A':58,'R':71,'N':65,'D':85,'C':149,'Q':42,'E':65,'G':59,'H':47,'I':89,'L':92,'K':78,'M':81,'F':103,'P':38,'S':58,'T':0,'W':128,'Y':92,'V':69},
    'W': {'A':148,'R':101,'N':174,'D':181,'C':215,'Q':130,'E':152,'G':184,'H':115,'I':61,'L':61,'K':110,'M':67,'F':40,'P':147,'S':177,'T':128,'W':0,'Y':37,'V':88},
    'Y': {'A':112,'R':77,'N':143,'D':160,'C':194,'Q':99,'E':122,'G':147,'H':83,'I':33,'L':36,'K':85,'M':36,'F':22,'P':110,'S':144,'T':92,'W':37,'Y':0,'V':55},
    'V': {'A':64,'R':96,'N':133,'D':152,'C':192,'Q':96,'E':121,'G':109,'H':84,'I':29,'L':32,'K':97,'M':21,'F':50,'P':68,'S':124,'T':69,'W':88,'Y':55,'V':0},
}

# Posiciones en numeracion 1-based de DENV1/2/4 (convencion de la literatura)
LOC_EDI = [[47, 47], [123, 126], [148, 165], [167, 167],
           [170, 177], [180, 180], [279, 279], [293, 293]]
LOC_EDII = [[67, 74], [76, 79], [88, 88], [97, 118],
            [215, 215], [219, 219], [233, 233], [237, 237],
            [244, 249], [256, 256], [266, 268], [272, 276]]
LOC_EDIII = [[282, 282], [298, 311], [314, 370],
             [372, 372], [382, 385], [387, 390]]
DOMINIOS = (LOC_EDI, LOC_EDII, LOC_EDIII)

W_GRANTHAM = 215      # maximo de la matriz, para normalizar (Ec. 2)


def to_std_numbering(seq):
    """Lleva una secuencia a la numeracion comun de DENV1/2/4.

    DENV3 tiene dos residuos menos por una delecion en 158-159, asi que se
    insertan dos huecos y la posicion k significa lo mismo en cualquier par.
    """
    if len(seq) == E_LEN_STD - 2:
        return seq[:DENV3_GAP] + "--" + seq[DENV3_GAP:]
    return seq


def grantham_distance(s1, s2):
    """Suma de Grantham sobre las columnas comparables.

    Los huecos y los stop se saltan. Un aminoacido que no este en la matriz
    (X, B, Z) tampoco suma, pero OJO: en antigenic_distance si cuenta en el
    denominador, porque alli el filtro solo descarta '-' y '*'.
    """
    if len(s1) != len(s2):
        raise ValueError("Sequences must be of the same length")
    distancia = 0
    for a, b in zip(s1, s2):
        if a in "-*" or b in "-*":
            continue
        if a not in GRANTHAM or b not in GRANTHAM[a]:
            continue
        distancia += GRANTHAM[a][b]
    return distancia


def antigenic_distance(s1, s2):
    """Distancia p-all-epitope entre dos secuencias del gen E.

    Ec. 2: E_i = (suma de G) / (n_i * 215) * 20 por dominio.
    Ec. 3: la distancia es el promedio de los tres dominios.
    """
    s1, s2 = to_std_numbering(s1), to_std_numbering(s2)

    def por_dominio(locations):
        total, n = 0, 0
        for ini, fin in locations:
            a, b = s1[ini - 1:fin], s2[ini - 1:fin]
            pares = [(x, y) for x, y in zip(a, b) if x not in "-*" and y not in "-*"]
            if not pares:
                continue
            total += grantham_distance("".join(x for x, _ in pares),
                                       "".join(y for _, y in pares))
            n += len(pares)
        return total, n

    es = []
    for loc in DOMINIOS:
        total, n = por_dominio(loc)
        es.append((total / (n * W_GRANTHAM)) * 20)
    return sum(es) / 3


# ─────────────────────── version vectorizada ───────────────────────
#
# La fila y la columna 20 son el comodin de los aminoacidos que no estan en la
# matriz: distancia 0 contra todo, pero cuentan en el denominador, igual que en
# antigenic_distance. El -1 marca hueco o stop, que no cuentan en ninguno.
_GM = np.zeros((21, 21), dtype=np.int32)
for i, a in enumerate(AA):
    for j, b in enumerate(AA):
        _GM[i, j] = GRANTHAM[a][b]

_COD = np.full(128, 20, dtype=np.int8)
for i, a in enumerate(AA):
    _COD[ord(a)] = i
_COD[ord("-")] = -1
_COD[ord("*")] = -1


def _columnas(locations):
    cols = []
    for ini, fin in locations:
        cols.extend(range(ini - 1, fin))
    return np.array(cols, dtype=int)


def matriz(ids, seqs):
    """Matriz simetrica de distancia p-all-epitope dentro de un conjunto.

    Devuelve un array (n, n). Equivale a llamar antigenic_distance en cada par,
    y autoprueba() lo verifica.
    """
    seqs = [to_std_numbering(s) for s in seqs]
    largo = {len(s) for s in seqs}
    if len(largo) != 1:
        raise ValueError("secuencias de distinta longitud tras la numeracion "
                         "comun: %s" % sorted(largo))
    cod = np.array([[_COD[ord(c)] if ord(c) < 128 else 20 for c in s]
                    for s in seqs], dtype=np.int8)

    n = len(seqs)
    acum = np.zeros((n, n), dtype=float)
    for loc in DOMINIOS:
        A = cod[:, _columnas(loc)]
        val = A >= 0
        amb = val[:, None, :] & val[None, :, :]
        idx = np.maximum(A, 0).astype(np.int16)
        g = _GM[idx[:, None, :], idx[None, :, :]]
        total = np.where(amb, g, 0).sum(-1)
        cuenta = amb.sum(-1)
        with np.errstate(divide="ignore", invalid="ignore"):
            acum += np.where(cuenta > 0,
                             total / (cuenta * W_GRANTHAM) * 20, np.nan)
    return acum / 3


def autoprueba(fasta, n_pares=200, semilla=0):
    """Compara la version vectorizada con la de referencia en pares al azar."""
    from Bio import SeqIO
    regs = [r for r in SeqIO.parse(fasta, "fasta")
            if not r.id.startswith("NODE")][:120]
    ids = [r.id for r in regs]
    seqs = [str(r.seq) for r in regs]
    M = matriz(ids, seqs)

    rng = np.random.default_rng(semilla)
    peor = 0.0
    for _ in range(n_pares):
        i, j = rng.integers(0, len(ids), 2)
        ref = antigenic_distance(seqs[i], seqs[j])
        peor = max(peor, abs(ref - M[i, j]))
    print("%s: %d secuencias, %d pares comprobados, diferencia maxima %.3g"
          % (fasta, len(ids), n_pares, peor))
    return peor


if __name__ == "__main__":
    rutas = sys.argv[1:] or [
        "resources/colombia_data/AA_sequences/E_DENV2_aa_E.fasta",
        "resources/colombia_data/AA_sequences/E_DENV3_aa_E.fasta",
    ]
    peor = max(autoprueba(r) for r in rutas)
    print("OK" if peor < 1e-9 else "DIFIEREN")
