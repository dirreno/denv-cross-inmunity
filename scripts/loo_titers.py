#!/usr/bin/env python3
"""Validacion leave-one-VIRUS-out del modelo de sustituciones de augur.

Responde la pregunta que ni el RMSE de ajuste ni la validacion cruzada de
Neher et al. (2016) responden: si se borra un virus ENTERO de la tabla de
titulos, cuanto se acerca la prediccion a sus titulos reales?

Es el nivel que corresponde al caso colombiano. Al apartar la fila completa
de un virus su avidez deja de ser estimable y cae a 0, exactamente como en una
secuencia sin medir. Ocultar celdas sueltas (matrix completion) conserva la
avidez y mide un problema mas facil que el nuestro.

POR QUE NO SE USA make_training_set(subset_strains=True) TAL CUAL
-----------------------------------------------------------------
Ese metodo hace una particion ALEATORIA por fraccion, no leave-one-out:

    training_strains = sample(sorted(tmp), int(training_fraction*len(tmp)))

Reparte por virus de prueba (key[0]), que es la granularidad correcta, pero no
permite elegir cual se aparta. Ademas excluye del sorteo a los virus de
referencia y los reinserta todos en entrenamiento, asi que nunca se validan.

Aqui se sustituye ese metodo en la instancia por uno que aparta UN virus
concreto, y se deja que prepare() siga su curso: la matriz de diseno y el
ajuste salen de augur sin tocar nada mas.

Se replica una decision de augur a proposito: los virus de referencia de los
sueros NO se apartan. Sin ellos el suero correspondiente se queda sin anclaje.

La CLI de augur no tiene modo leave-one-out, asi que esto usa la API de Python.
Solo necesita que 'import augur' funcione en el mismo entorno que la CLI.

USO
---
    python3 scripts/loo_titers.py --inspect          # solo comprobar la API
    python3 scripts/loo_titers.py --sero 1 --limit 5 # prueba rapida
    python3 scripts/loo_titers.py                    # los 4 serotipos
    python3 scripts/loo_titers.py --repo /ruta/al/repo

VERSION
-------
Lo que sigue se escribio leyendo la fuente de augur en 'master'. El documento
de contexto del proyecto dice 33.0.0; la documentacion que se estaba
consultando era de 10.0.4. --inspect imprime la version real y la fuente
efectiva de los metodos que importan, para poder contrastarlo antes de creer
ningun numero.
"""
import argparse
import csv
import inspect
import os
import sys
import types

import numpy as np

TITERS = "resources/merged_data/titers/titers_DENV{}.tsv"
ALIGN = "resources/merged_data/AA_sequences/E_DENV{}_aa_E.fasta"
TREE = "resources/merged_data/refined/E_DENV{}_tree_time.nwk"
OUT_PAIRS = "output/loo_titers_pairs.csv"
OUT_SUMMARY = "output/loo_titers_summary.csv"
GENE = "E"


# --------------------------------------------------------------- comprobacion
def cargar_titer_model():
    """Devuelve el modulo titer_model, sea cual sea la version de augur."""
    import importlib
    for ruta in ("augur.titer_model", "augur.titer_models", "augur.titers"):
        try:
            return importlib.import_module(ruta)
        except ImportError:
            continue
    raise SystemExit(
        "No encuentro el modulo de modelos de titulos dentro de augur.\n"
        "Prueba:  python3 -c 'import augur,pkgutil;"
        "print([m.name for m in pkgutil.iter_modules(augur.__path__)])'")


def inspeccionar():
    """Imprime la API real antes de usarla.

    No es ceremonia: la fuente que se leyo al escribir esto venia de 'master',
    y este script va a correr contra una version que puede ser 10.x o 33.x.
    Si las firmas no coinciden, se ve aqui y no en un numero silenciosamente
    equivocado.
    """
    import augur
    tm = cargar_titer_model()

    ver = getattr(augur, "__version__", "?")
    print(f"augur {ver}  ->  {tm.__file__}\n")
    try:
        mayor = int(str(ver).split(".")[0])
        if mayor < 20:
            print(f"AVISO: version {ver}. Los hechos verificados en "
                  "docs/CONTEXTO_sesion_titer_models.md son de 33.0.0;\n"
                  "       contrasta las firmas de abajo antes de usar los "
                  "resultados.\n")
    except (ValueError, IndexError):
        pass
    for nombre in ("TiterModel", "SubstitutionModel", "TreeModel", "TiterCollection"):
        cls = getattr(tm, nombre, None)
        print(f"{nombre}: {'OK' if cls else 'NO EXISTE'}")
        if cls is None:
            continue
        for met in ("__init__", "prepare", "train", "make_training_set",
                    "predict_titer", "validate"):
            f = getattr(cls, met, None)
            if f is not None and callable(f):
                try:
                    print(f"    {met}{inspect.signature(f)}")
                except (ValueError, TypeError):
                    print(f"    {met}(?)")
    print("\n--- make_training_set, fuente real ---")
    try:
        print(inspect.getsource(tm.TiterModel.make_training_set))
    except OSError:
        print("(no se pudo leer la fuente)")
    print("--- SubstitutionModel.prepare, fuente real ---")
    try:
        print(inspect.getsource(tm.SubstitutionModel.prepare))
    except (OSError, AttributeError):
        print("(no se pudo leer la fuente)")


# ------------------------------------------------------------------- partición
def split_dejando_fuera(strain):
    """Sustituto de make_training_set que aparta UN virus concreto.

    Replica el cierre del metodo original: recalcula el censo de cepas sobre el
    conjunto de entrenamiento, que es lo que prepare() usa despues.
    """
    def _make_training_set(self, *a, **kw):
        self.test_titers, self.train_titers = {}, {}
        for key, val in self.titers.titers_normalized.items():
            if key[0] == strain:
                self.test_titers[key] = val
            else:
                self.train_titers[key] = val
        self.sera, self.ref_strains, self.test_strains = \
            self.titers.strain_census(self.train_titers)
    return _make_training_set


def construir(tm, sero):
    """Crea un SubstitutionModel. La firma se comprueba antes de llamar."""
    aln, tit, tree = ALIGN.format(sero), TITERS.format(sero), TREE.format(sero)
    for p in (aln, tit, tree):
        if not os.path.exists(p):
            raise SystemExit(f"falta {p} (correr desde la raiz del repo)")
    params = inspect.signature(tm.SubstitutionModel.__init__).parameters
    kw = {}
    if "alignment" in params:
        kw["alignment"] = aln
    if "titers" in params:
        kw["titers"] = tit
    if "gene_names" in params:
        kw["gene_names"] = [GENE]
    if "tree" in params:
        kw["tree"] = tree
    if not kw:
        raise SystemExit(f"firma inesperada de SubstitutionModel.__init__: {params}")
    return tm.SubstitutionModel(**kw)


def candidatos(modelo):
    """Virus que se pueden apartar: los de prueba que NO son referencia."""
    ref = set(getattr(modelo, "ref_strains", []) or [])
    test = set(getattr(modelo, "test_strains", []) or [])
    return sorted(test - ref), sorted(ref)


# ------------------------------------------------------------------------ LOO
def correr(sero, limite=None, verbose=True):
    tm = cargar_titer_model()

    base = construir(tm, sero)
    base.prepare(training_fraction=1.0)
    fuera, refs = candidatos(base)
    if limite:
        fuera = fuera[:limite]
    print(f"\n=== DENV{sero}: {len(fuera)} virus a apartar "
          f"({len(refs)} referencias se conservan siempre) ===")

    filas, fallos = [], 0
    for n, strain in enumerate(fuera, 1):
        try:
            m = construir(tm, sero)
            m.make_training_set = types.MethodType(split_dejando_fuera(strain), m)
            m.prepare(training_fraction=1.0)
            m.train()
            # modelo nulo: solo potencia del suero (avidez del virus apartado = 0)
            for key, val in m.test_titers.items():
                virus, suero = key[0], key[1]
                obs = val[0] if isinstance(val, (list, tuple)) else val
                pred = m.predict_titer(virus, suero)
                pot = m.serum_potency.get(suero) if hasattr(m, "serum_potency") else None
                filas.append(dict(serotype=sero, strain=virus,
                                  ref_virus=suero[0] if isinstance(suero, tuple) else "",
                                  serum_id=suero[1] if isinstance(suero, tuple) else suero,
                                  observado=obs, predicho=pred, nulo=pot))
        except Exception as e:                      # noqa: BLE001
            fallos += 1
            if verbose and fallos <= 3:
                print(f"  [fallo] {strain}: {type(e).__name__}: {e}")
        if verbose and n % 10 == 0:
            print(f"  {n}/{len(fuera)}  pares acumulados: {len(filas)}")
    if fallos:
        print(f"  {fallos} virus no se pudieron ajustar (se omiten)")
    return filas


def resumir(filas):
    """RMSE del modelo frente al nulo, en log2."""
    out = []
    for sero in sorted({f["serotype"] for f in filas}):
        sub = [f for f in filas if f["serotype"] == sero
               and f["predicho"] is not None and f["observado"] is not None]
        if not sub:
            continue
        obs = np.array([float(f["observado"]) for f in sub])
        pre = np.array([float(f["predicho"]) for f in sub])
        nul = np.array([float(f["nulo"]) if f["nulo"] is not None else np.nan
                        for f in sub])
        rmse = float(np.sqrt(((pre - obs) ** 2).mean()))
        ok = ~np.isnan(nul)
        rmse_n = float(np.sqrt(((nul[ok] - obs[ok]) ** 2).mean())) if ok.any() else np.nan
        out.append(dict(serotype=sero, n_pares=len(sub),
                        n_virus=len({f["strain"] for f in sub}),
                        rmse_loo=rmse, rmse_nulo=rmse_n,
                        ganancia=rmse_n - rmse if rmse_n == rmse_n else np.nan,
                        skill=1 - rmse / rmse_n if rmse_n == rmse_n and rmse_n else np.nan))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inspect", action="store_true",
                    help="solo comprobar la API de augur y salir")
    ap.add_argument("--sero", type=int, nargs="*", default=[1, 2, 3, 4])
    ap.add_argument("--limit", type=int, default=None,
                    help="apartar solo los primeros N virus (prueba rapida)")
    ap.add_argument("--repo", default=".",
                    help="raiz del repositorio (por defecto, el directorio actual)")
    a = ap.parse_args()

    if a.repo != ".":
        os.chdir(a.repo)
    if not os.path.isdir("resources"):
        sys.exit(f"No veo 'resources/' en {os.getcwd()}. Usa --repo /ruta/al/repo")

    try:
        import augur  # noqa: F401
    except ImportError:
        sys.exit("augur no esta en este entorno de Python. Activa el mismo "
                 "entorno desde el que corre la CLI de augur.")

    inspeccionar()
    if a.inspect:
        return

    filas = []
    for s in a.sero:
        filas += correr(s, limite=a.limit)
    if not filas:
        sys.exit("\nNo se obtuvo ningun par. Revisa la salida de --inspect.")

    os.makedirs("output", exist_ok=True)
    with open(OUT_PAIRS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)

    res = resumir(filas)
    with open(OUT_SUMMARY, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0]))
        w.writeheader()
        w.writerows(res)

    print("\n" + "=" * 78)
    print("LEAVE-ONE-VIRUS-OUT: predice los titulos de un virus nunca medido?")
    print("=" * 78)
    print(f"{'':7s} {'virus':>6s} {'pares':>7s} {'RMSE LOO':>9s} "
          f"{'RMSE nulo':>10s} {'ganancia':>9s} {'skill':>7s}")
    for r in res:
        print(f"DENV{r['serotype']}  {r['n_virus']:6d} {r['n_pares']:7d} "
              f"{r['rmse_loo']:9.3f} {r['rmse_nulo']:10.3f} "
              f"{r['ganancia']:9.3f} {r['skill']:7.3f}")
    print("\nTodo en unidades log2.")
    print("  'nulo'  = solo potencia del suero (avidez del virus apartado = 0).")
    print("  'skill' = 1 - RMSE_LOO/RMSE_nulo. <= 0 significa que el modelo no")
    print("            aporta nada sobre predecir el promedio de ese suero.")
    print("  Referencia: el ruido del ensayo ronda 1 log2. Un RMSE cercano a 1")
    print("  no se puede mejorar aunque el modelo fuera perfecto.")
    print(f"\npares  -> {OUT_PAIRS}")
    print(f"resumen-> {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
