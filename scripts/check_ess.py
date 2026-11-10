import sys, math, glob, os


def ess(x):
    """ESS estilo Tracer: suma de autocorrelaciones hasta que el par se hace negativo."""
    n = len(x)
    if n < 10:
        return 0.0
    m = sum(x) / n
    d = [v - m for v in x]
    var = sum(v * v for v in d) / n
    if var <= 0:
        return 0.0
    maxlag = min(n - 1, 2000)
    gamma = []
    for lag in range(maxlag):
        s = sum(d[i] * d[i + lag] for i in range(n - lag)) / n
        gamma.append(s / var)
    total = gamma[0]
    k = 1
    while k + 1 < len(gamma):
        par = gamma[k] + gamma[k + 1]
        if par <= 0:
            break
        total += 2 * par
        k += 2
    return n / total if total > 0 else 0.0


for f in sorted(glob.glob(sys.argv[1])):
    filas, hdr = [], None
    for line in open(f):
        if line.startswith("Sample"):
            hdr = line.rstrip("\n").split("\t")
            continue
        if hdr and line.strip() and not line.startswith("#"):
            p = line.rstrip("\n").split("\t")
            if len(p) == len(hdr):
                filas.append(p)
    if not hdr or len(filas) < 20:
        print(f"{os.path.basename(os.path.dirname(f)):<28} muestras insuficientes ({len(filas)})")
        continue
    # Los *_combined.log ya salen de LogCombiner con -b 10: recortar aqui otro
    # 10 % tiraba muestras buenas y bajaba el ESS entre un 4 % y un 12 %. Los
    # .log crudos de runs/ si llevan la cadena completa y hay que recortarlos.
    # Se puede forzar con un segundo argumento: la fraccion de burn-in (0 a 1).
    nombre = os.path.basename(f)
    if len(sys.argv) > 2:
        frac = float(sys.argv[2])
        motivo = "forzado por argumento"
    elif "_combined" in nombre:
        frac = 0.0
        motivo = "combinado, LogCombiner ya descarto el burn-in"
    else:
        frac = 0.10
        motivo = "cadena cruda"
    burn = int(len(filas) * frac)
    filas = filas[burn:]
    tag = nombre.replace(".log", "")
    print(f"\n=== {tag}  ({len(filas)} muestras, burn-in {frac:.0%} "
          f"[{motivo}]; estado {filas[-1][0]})")
    # BEAST llama a la altura del arbol "Tree.height", no "TreeHeight".
    for col in ["posterior", "treeLikelihood", "Tree.height", "origin", "delta", "ucldMean",
                "Re.1", "Re.2", "sProp.2"]:
        if col not in hdr:
            continue
        i = hdr.index(col)
        try:
            v = [float(r[i]) for r in filas]
        except ValueError:
            continue
        e = ess(v)
        marca = "OK" if e >= 200 else ("bajo" if e >= 100 else "MUY BAJO")
        print(f"   {col:<12} ESS={e:7.1f}  {marca}   media={sum(v)/len(v):.4g}")
