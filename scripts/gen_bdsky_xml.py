# -*- coding: utf-8 -*-
"""Genera los XML de BEAST2/BDSKY, uno por unidad de analisis.

La unidad es el linaje de la metadata, salvo cuando un linaje contiene dos
episodios separados por un hueco de muestreo: entonces se usa un clado del
arbol datado de augur, nombrado por su nodo (clave `clado`).
"""
import csv, json, math, os, re

B   = "D:/EQUIPO/Documents HDD/Tesis/denv-cross-inmunity"
OUT = B + "/resources/BEAST/lineage"
os.makedirs(OUT, exist_ok=True)


def leer_fasta(path, corta_pipe=False):
    seqs, n, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if n:
                seqs[n] = "".join(buf)
            h = line[1:].split()[0]
            n, buf = (h.split("|")[0] if corta_pipe else h), []
        else:
            buf.append(line.upper())
    if n:
        seqs[n] = "".join(buf)
    return seqs


def decimal(fecha):
    """'2023-04-17' / '2023-04' / '2023' -> anio decimal."""
    p = fecha.split("-")
    y = int(p[0])
    bis = (y % 4 == 0 and y % 100 != 0) or y % 400 == 0
    dias = 366 if bis else 365
    acum = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    if bis:
        acum = [d + (1 if i >= 2 else 0) for i, d in enumerate(acum)]
    if len(p) == 1:
        return y + 0.5
    if len(p) == 2:
        return y + (acum[int(p[1]) - 1] + 15) / dias
    return y + (acum[int(p[1]) - 1] + int(p[2]) - 0.5) / dias


def xml_safe(s):
    return re.sub(r"[^A-Za-z0-9._/-]", "_", s)


def puntas_del_clado(sero, nodo):
    """Nombres de punta bajo un nodo del arbol datado de augur refine."""
    s = open(B + "/resources/colombia_data/refined/E_DENV%d_tree_time.nwk" % sero).read().strip()
    pos = [0]

    def leer():
        if s[pos[0]] == "(":
            pos[0] += 1
            hijos = [leer()]
            while s[pos[0]] == ",":
                pos[0] += 1
                hijos.append(leer())
            pos[0] += 1                                   # cierra ')'
            m = re.match(r"[^,\)\;]*", s[pos[0]:])
            etq = m.group(0).split(":")[0].strip("'\"")
            pos[0] += m.end()
            return {"n": etq, "h": hijos}
        m = re.match(r"[^,\)\;]*", s[pos[0]:])
        etq = m.group(0).split(":")[0].strip("'\"")
        pos[0] += m.end()
        return {"n": etq, "h": []}

    def hojas(n):
        if not n["h"]:
            return [n["n"]]
        out = []
        for h in n["h"]:
            out += hojas(h)
        return out

    def buscar(n):
        if n["n"] == nodo:
            return n
        for h in n["h"]:
            r = buscar(h)
            if r:
                return r
        return None

    raiz = leer()
    objetivo = buscar(raiz)
    if objetivo is None:
        raise SystemExit("nodo %s no esta en el arbol de DENV%d" % (nodo, sero))
    return set(hojas(objetivo))


# re_rev / sp_rev: tiempos de cambio HACIA ATRAS desde la muestra mas reciente.
# Deben incluir el 0 y tener tantas entradas como dimension del parametro.
ANALISIS = [
    dict(tag="DENV2_2III_D2", sero=2, linaje="2III_D.2", desde=2021, aln="aligned",
         re_rev=[0, 0.6, 1.2, 1.8, 2.4, 3.0, 3.6], sp_rev=[0, 3.5],
         clock_M=math.log(9.0e-4), clock_S=0.35, org=18.5),
    dict(tag="DENV2_2II_F1112", sero=2, linaje="2II_F.1.1.2", desde=None, aln="aligned",
         re_rev=[0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8], sp_rev=[0, 2.5],
         clock_M=math.log(9.0e-4), clock_S=0.35, org=9.0),
    # 3III_C.1 NO se analiza como clado completo. Sus 155 puntas son dos
    # episodios con seis anios sin una sola secuencia entre ellos (2017-2022):
    # el bloque historico 2001-2016 y un clado monofiletico de 76 puntas de
    # 2023-24 (NODE_0000187, TMRCA 2015.35) cuyas hermanas son colombianas de
    # 2013-15. Un solo proceso de nacimiento-muerte sobre ambos obliga al
    # modelo a sostener la cadena a traves del hueco, que es la via por la que
    # Re vuelve a colapsar a 1. Se separan en dos analisis.
    dict(tag="DENV3_3III_C1_hist", sero=3, linaje="3III_C.1", hasta=2023, aln="cropped",
         re_rev=[0, 3, 6, 9, 12, 15], sp_rev=[0, 15.1],
         clock_M=math.log(7.64e-4), clock_S=0.25, org=19.0),
    dict(tag="DENV3_3III_C1_2023", sero=3, clado="NODE_0000187", aln="cropped",
         re_rev=[0, 0.3, 0.6, 0.9, 1.2, 1.5], sp_rev=[0, 1.45],
         clock_M=math.log(7.64e-4), clock_S=0.30, org=11.0),
    # Reloj recentrado: con LogN(log 7.64e-4, 0.25) el 68% de la masa posterior
    # quedaba por encima del percentil 97.5 del prior (post 1.33e-3, prior95
    # hasta 1.25e-3), con ESS>1200 y las dos semillas de acuerdo: el prior era
    # restriccion activa y eso era una mala especificacion que habia que
    # corregir. OJO: el analisis de sensibilidad posterior mostro que NO era la
    # causa del Re plano. Con el prior ancho el reloj subio a 1.50e-3 y los seis
    # Re salieron identicos al tercer decimal. El Re plano de C.2 viene de las
    # 64 fechas de solo anio sobre 109, no del prior del reloj.
    dict(tag="DENV3_3III_C2", sero=3, linaje="3III_C.2", desde=None, aln="cropped",
         re_rev=[0, 2.3, 4.6, 6.9, 9.2, 11.5], sp_rev=[0, 13.0],
         clock_M=math.log(1.3e-3), clock_S=0.40, org=15.0),
    # Sensibilidad del prior de delta, identica a DENV3_3III_C1_hist salvo S.
    # Con LogN(log 26.07, 0.5) el 50% de la masa posterior de delta quedaba por
    # DEBAJO del percentil 2.5 del prior (post 10.41/anio = 35 dias, prior2.5
    # 9.78/anio), con ESS>2600 y las dos semillas de acuerdo: el prior es
    # restriccion activa en el borde inferior. Se mantiene la mediana de 14 dias
    # y solo se ensancha S a 1.0 (95% [3.67, 185]/anio) para dejar que el dato
    # mueva delta si quiere, y se mira si Re se mueve con el.
    dict(tag="DENV3_3III_C1_hist_delta", sero=3, linaje="3III_C.1", hasta=2023,
         aln="cropped", re_rev=[0, 3, 6, 9, 12, 15], sp_rev=[0, 15.1],
         clock_M=math.log(7.64e-4), clock_S=0.25, org=19.0, delta_S=1.0),
]

CADENA = 100000000
DELTA_M = math.log(365.0 / 14.0)   # 14 dias -> 26.07/anio (Faria 2017)
DELTA_S = 0.5

PLANTILLA = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "gen_bdsky_plantilla.xml"), encoding="utf-8").read()

for A in ANALISIS:
    s, tag = A["sero"], A["tag"]
    corta = (A["aln"] == "aligned")
    S = leer_fasta(B + "/resources/colombia_data/%s/align_DENV%d.fasta" % (A["aln"], s),
                   corta_pipe=corta)
    lin, fec = {}, {}
    with open(B + "/resources/colombia_data/metadata/E_DENV%d_meta.tsv" % s,
              encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lin[r["strain"]] = r.get("lineage", "")
            fec[r["strain"]] = r.get("date", "")

    if A.get("clado"):
        sel = puntas_del_clado(s, A["clado"])
        ids = [k for k in S if k in sel and fec.get(k)]
    else:
        ids = [k for k in S if lin.get(k) == A["linaje"] and fec.get(k)]
    if A.get("desde"):
        ids = [k for k in ids if int(fec[k][:4]) >= A["desde"]]
    if A.get("hasta"):
        ids = [k for k in ids if int(fec[k][:4]) < A["hasta"]]
    ids = [k for k in ids if sum(1 for c in S[k] if c in "ACGT") / len(S[k]) >= 0.30]
    ids.sort()
    fechas = {k: decimal(fec[k]) for k in ids}
    reciente = max(fechas.values())
    primera = min(fechas.values())
    altura = reciente - primera
    L = len(S[ids[0]])

    seqs = "\n".join(
        '        <sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>'
        % (xml_safe(k), xml_safe(k), S[k]) for k in ids)
    trait = ",".join("%s=%.5f" % (xml_safe(k), fechas[k]) for k in ids)

    nre, nsp = len(A["re_rev"]), len(A["sp_rev"])
    campos = {
        "TAG": tag,
        "LINAJE": A.get("linaje") or ("clado " + A["clado"]),
        "N": len(ids), "L": L,
        "PRIMERA": "%.3f" % primera, "RECIENTE": "%.3f" % reciente,
        "ALTURA": "%.2f" % altura,
        "SEQS": seqs, "TRAIT": trait, "CADENA": CADENA,
        "NRE": nre, "NSP": nsp,
        "RET": " ".join("%g" % t for t in A["re_rev"]),
        "SPT": " ".join("%g" % t for t in A["sp_rev"]),
        "SPVAL": " ".join(["1.0E-10"] + ["0.001"] * (nsp - 1)),
        "SPIND": " ".join(["false"] + ["true"] * (nsp - 1)),
        # debe superar la altura del arbol inicial de RandomTree, no solo el prior
        "ORG0": "%.2f" % (max(A["org"], altura + 2.0) * 2.0),
        "ORGM": "%.4f" % math.log(A["org"]),
        "CLOCK0": "%.3e" % math.exp(A["clock_M"]),
        "CLOCKM": "%.4f" % A["clock_M"], "CLOCKS": A["clock_S"],
        "DELTAM": "%.4f" % A.get("delta_M", DELTA_M),
        "DELTAS": A.get("delta_S", DELTA_S),
        "NCAT": 2 * len(ids) - 2,
        # popSize pequeno -> RandomTree poco mas alto que el rango de fechas
        "POPSIZE": "1.0",
        "SLIDE": "%.3f" % max(0.1, altura / 20),
    }
    xml = PLANTILLA
    for k, v in campos.items():
        xml = xml.replace("{{%s}}" % k, str(v))

    ruta = OUT + "/BDSKY_%s.xml" % tag
    open(ruta, "w", encoding="utf-8").write(xml)
    print("%-18s n=%-4d %5d nt  altura=%5.2f a  Re=%d  s=%d  -> %s (%.1f MB)"
          % (tag, len(ids), L, altura, nre, nsp, os.path.basename(ruta),
             os.path.getsize(ruta) / 1e6))
