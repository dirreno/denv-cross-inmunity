<#
.SYNOPSIS
    Combina las semillas de cada analisis BDSKY con LogCombiner y resume el
    conjunto de arboles con TreeAnnotator.

.DESCRIPTION
    Toma solo s1 y s2 de cada analisis. Los directorios s3 que quedaron de un
    lote abortado contienen cadenas de menos del 5 % y NO deben entrar al
    combinado: el comando de ejemplo que imprime run_bdsky.ps1 los incluye
    literalmente, asi que aqui se enumeran las semillas de forma explicita.

    El burn-in se aplica UNA sola vez, en LogCombiner. TreeAnnotator se invoca
    con -burnin 0 porque su valor por defecto es 10 y volveria a recortar.

    Los arboles se submuestrean a -resample para el resumen: con logEvery=10000
    y 100 M de estados cada cadena aporta 10.001 arboles, y 20.000 deja ~4.500
    por semilla, mas que suficiente para un MCC. Las trazas NO se submuestrean,
    para no perder ESS.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File scripts\combinar_bdsky.ps1
    powershell -ExecutionPolicy Bypass -File scripts\combinar_bdsky.ps1 -Solo DENV2
#>
param(
    [string]   $BeastDir       = "D:\EQUIPO\Documents HDD\Tesis\programs\BEAST.v2.7.7.Windows\BEAST",
    [string[]] $Seeds          = @("1", "2"),
    [int]      $Burnin         = 10,
    [int]      $ResampleTrees  = 20000,
    [string]   $Xmx            = "12g",
    [string[]] $Solo           = @(),
    [switch]   $SoloTrazas,
    [switch]   $SoloArboles
)

$ErrorActionPreference = "Stop"

$raiz     = Split-Path -Parent $PSScriptRoot
$xmlDir   = Join-Path $raiz "resources\BEAST\lineage"
$runsDir  = Join-Path $xmlDir "runs"
$outDir   = Join-Path $xmlDir "combined"
$java     = Join-Path $BeastDir "jre\bin\java.exe"
$launcher = Join-Path $BeastDir "lib\launcher.jar"

foreach ($p in @($java, $launcher, $runsDir)) {
    if (-not (Test-Path $p)) { throw "No existe: $p" }
}
if (-not (Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

# Los seis analisis de la tesis. Se excluyen a proposito:
#   BDSKY_DENV3_3III_C1            diseno viejo, mezcla dos episodios, 40 M de 100 M
#   BDSKY_DENV3_3III_C2_relojviejo replica archivada del prior de reloj anterior
$analisis = @(
    "BDSKY_DENV2_2II_F1112",
    "BDSKY_DENV2_2III_D2",
    "BDSKY_DENV3_3III_C1_2023",
    "BDSKY_DENV3_3III_C2",
    "BDSKY_DENV3_3III_C1_hist",
    "BDSKY_DENV3_3III_C1_hist_delta"
)

if ($Solo.Count -gt 0) {
    $analisis = @($analisis | Where-Object {
        $tag = $_; $hit = $false
        foreach ($p in $Solo) { if ($tag -like "*$p*") { $hit = $true } }
        $hit
    })
}

function Ejecutar($clase, $argumentos, $etiqueta, $esperado) {
    $todos = @("-Xmx$Xmx", "-Dfile.encoding=UTF-8", "-cp", $launcher, $clase) + $argumentos
    if ($esperado -and (Test-Path $esperado)) { Remove-Item $esperado -Force }
    $t0 = Get-Date
    # OJO: NO redirigir stderr con 2>&1. LogCombiner y TreeAnnotator escriben su
    # banner y su progreso en stderr, y en PowerShell 5.1 redirigir el stderr de
    # un nativo envuelve cada linea en un ErrorRecord; con ErrorActionPreference
    # en Stop eso aborta el script aunque la herramienta termine en 0.
    $previo = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try {
        $salida = & $java @todos
        $codigo = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $previo
    }
    $dt = ((Get-Date) - $t0).TotalSeconds
    # TreeAnnotator puede terminar en 0 sin escribir nada si no parsea los
    # argumentos, asi que el codigo de salida no basta: se exige el archivo.
    $falta = $esperado -and -not (Test-Path $esperado)
    if ($codigo -ne 0 -or $falta) {
        $motivo = if ($falta) { "no escribio $(Split-Path -Leaf $esperado)" } else { "codigo $codigo" }
        Write-Host ("    FALLO ($motivo, {0:N0} s)" -f $dt) -ForegroundColor Red
        $salida | Select-Object -Last 20 | ForEach-Object { Write-Host "      $_" }
        return $false
    }
    Write-Host ("    ok {0,6:N0} s  {1}" -f $dt, $etiqueta) -ForegroundColor Green
    return $true
}

Write-Host "Combinando $($analisis.Count) analisis, semillas $($Seeds -join ',')"
Write-Host "  burn-in $Burnin %   resample de arboles $ResampleTrees   Xmx $Xmx"
Write-Host "  salida: $outDir`n"

$resumen = @()

foreach ($tag in $analisis) {
    Write-Host "=== $tag" -ForegroundColor Cyan

    $logs  = @()
    $trees = @()
    foreach ($s in $Seeds) {
        $d = Join-Path $runsDir ("{0}\s{1}" -f $tag, $s)
        $l = Join-Path $d "$tag.log"
        $t = Join-Path $d "$tag.trees"
        if (-not (Test-Path $l)) { Write-Host "    falta $l" -ForegroundColor Yellow; continue }
        # solo entran cadenas completas
        $ultimo = (Get-Content $l -Tail 1) -split "`t" | Select-Object -First 1
        if ([double]$ultimo -lt 99000000) {
            Write-Host "    s$s incompleta (estado $ultimo), se omite" -ForegroundColor Yellow
            continue
        }
        $logs  += $l
        $trees += $t
    }
    if ($logs.Count -lt 2) {
        Write-Host "    menos de dos semillas completas, se salta" -ForegroundColor Yellow
        continue
    }

    $outLog  = Join-Path $outDir "$tag`_combined.log"
    $outTree = Join-Path $outDir "$tag`_combined.trees"
    $outMcc  = Join-Path $outDir "$tag`_MCC.tree"
    $okL = $true; $okT = $true; $okM = $true

    if (-not $SoloArboles) {
        $a = @()
        foreach ($l in $logs) { $a += @("-log", $l) }
        $a += @("-b", $Burnin, "-renumber", "-o", $outLog)
        $okL = Ejecutar "beast.pkgmgmt.launcher.LogCombinerLauncher" $a "traza -> $(Split-Path -Leaf $outLog)" $outLog
    }

    if (-not $SoloTrazas) {
        $a = @()
        foreach ($t in $trees) { $a += @("-log", $t) }
        $a += @("-b", $Burnin, "-resample", $ResampleTrees, "-renumber", "-o", $outTree)
        $okT = Ejecutar "beast.pkgmgmt.launcher.LogCombinerLauncher" $a "arboles -> $(Split-Path -Leaf $outTree)" $outTree

        if ($okT) {
            # TreeAnnotator recibe los dos nombres como argumentos posicionales y
            # no los resuelve bien cuando la ruta absoluta lleva espacios
            # ("Documents HDD"): termina en 0 sin escribir nada. Se ejecuta con el
            # directorio de salida como cwd y solo los nombres de archivo.
            # -burnin 0 porque LogCombiner ya recorto; el defecto seria 10.
            Push-Location $outDir
            try {
                $a = @("-burnin", 0, "-height", "CA",
                       (Split-Path -Leaf $outTree), (Split-Path -Leaf $outMcc))
                $okM = Ejecutar "beast.pkgmgmt.launcher.TreeAnnotatorLauncher" $a "MCC -> $(Split-Path -Leaf $outMcc)" $outMcc
            } finally {
                Pop-Location
            }
        }
    }

    $resumen += [pscustomobject]@{
        Analisis = $tag
        Traza    = if ($okL -and (Test-Path $outLog))  { "{0:N0} filas" -f ((Get-Content $outLog | Measure-Object -Line).Lines) } else { "-" }
        Arboles  = if ($okT -and (Test-Path $outTree)) { "{0:N0} MB" -f ((Get-Item $outTree).Length/1MB) } else { "-" }
        MCC      = if ($okM -and (Test-Path $outMcc))  { "ok" } else { "-" }
    }
}

Write-Host "`n=== resumen"
$resumen | Format-Table -AutoSize
Write-Host "Verificar el ESS de los combinados con:"
Write-Host "  python scripts\check_ess.py `"resources/BEAST/lineage/combined/*_combined.log`""
