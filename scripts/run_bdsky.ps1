<#
    Lanza los cuatro analisis BDSKY con varias cadenas independientes en paralelo.

    Cada combinacion analisis x semilla corre en su propio directorio, asi que los
    .log y .trees no se pisan. Al terminar, combinar las cadenas de un mismo
    analisis con LogCombiner (ver el final de este archivo).

    Uso (Windows PowerShell 5.1):
        powershell -ExecutionPolicy Bypass -File scripts\run_bdsky.ps1
        powershell -ExecutionPolicy Bypass -File scripts\run_bdsky.ps1 -Seeds 1,2,3,4
        powershell -ExecutionPolicy Bypass -File scripts\run_bdsky.ps1 -Threads 2 -MaxParalelo 6
#>
param(
    # Con -File, PowerShell pasa "1,2" como UNA cadena en vez de un array, asi
    # que estos llegan como [string[]] y se normalizan mas abajo.
    [string[]] $Seeds    = @("1", "2"),
    [int]   $Threads     = 1,
    [int]   $MaxParalelo = 8,
    [string]$XmxPorProceso = "2g",
    # cpu = BEAGLE forzado a CPU con SSE. Es el predeterminado porque en esta
    #   maquina el modo automatico elige OpenCL sobre la Radeon integrada y va
    #   ~14x mas lento (61 min/M frente a 4,4 medidos en 3III_C2).
    # gpu = fuerza la RTX 4070: solo gana ~10% en el analisis mas pesado y nada
    #   en los demas, y una sola tarjeta no sirve para 8 cadenas a la vez.
    # java = sin BEAGLE. auto = dejar elegir a BEAST (NO usar aqui).
    [ValidateSet("cpu", "gpu", "java", "auto")]
    [string]$Beagle      = "cpu",
    # Patrones para correr solo algunos analisis. CUIDADO: el script pasa
    # -overwrite, asi que lanzarlo sin filtro mientras hay cadenas en marcha
    # borra sus .log. Ej: -Solo DENV3_3III_C1_hist,DENV3_3III_C1_2023
    [string[]]$Solo      = @(),
    # Permite relanzar un analisis que ya esta corriendo, borrando su .log.
    [switch]$Pisar,
    [string]$BeastDir    = "D:\EQUIPO\Documents HDD\Tesis\programs\BEAST.v2.7.7.Windows\BEAST"
)

$ErrorActionPreference = "Stop"

# Normaliza los argumentos de lista: con -File llegan como una sola cadena.
function Lista($v) { @($v | ForEach-Object { $_ -split "," } | Where-Object { $_ -ne "" }) }
$Seeds = @(Lista $Seeds | ForEach-Object { [int]$_ })
$Solo  = Lista $Solo

$raiz     = Split-Path -Parent $PSScriptRoot
$xmlDir   = Join-Path $raiz "resources\BEAST\lineage"
$runsDir  = Join-Path $xmlDir "runs"
$java     = Join-Path $BeastDir "jre\bin\java.exe"
$launcher = Join-Path $BeastDir "lib\launcher.jar"

foreach ($p in @($java, $launcher, $xmlDir)) {
    if (-not (Test-Path $p)) { throw "No existe: $p" }
}
New-Item -ItemType Directory -Force -Path $runsDir | Out-Null

# Los de DENV2 son 10.723 nt y tardan mucho mas: van primero para que no
# queden solos al final ocupando un nucleo mientras el resto ya termino.
$analisis = @(
    "BDSKY_DENV2_2II_F1112",
    "BDSKY_DENV2_2III_D2",
    "BDSKY_DENV3_3III_C1_2023",
    "BDSKY_DENV3_3III_C2",
    "BDSKY_DENV3_3III_C1_hist",
    "BDSKY_DENV3_3III_C1_hist_delta"
)

if ($Solo.Count -gt 0) {
    $todos = $analisis
    $analisis = @($todos | Where-Object {
        $tag = $_
        $hit = $false
        foreach ($p in $Solo) { if ($tag -like "*$p*") { $hit = $true } }
        $hit
    })
    if ($analisis.Count -eq 0) {
        throw ("El filtro -Solo ({0}) no coincide con ninguno de: {1}" -f
               ($Solo -join ", "), ($todos -join ", "))
    }
}

# El script pasa -overwrite, asi que relanzar un analisis que ya esta en marcha
# le borra el .log. Se comprueba la linea de comando de cada java vivo y se
# aborta solo si hay colision real; -Pisar la ignora a proposito.
$enMarcha = @(Get-CimInstance Win32_Process -Filter "Name='java.exe'" |
    ForEach-Object {
        $m = [regex]::Match($_.CommandLine, 'BDSKY_[A-Za-z0-9_]+\.xml')
        if ($m.Success) { $m.Value -replace '\.xml$', '' }
    })
$choque = @($analisis | Where-Object { $enMarcha -contains $_ })
if ($choque.Count -gt 0 -and -not $Pisar) {
    throw (("Ya estan corriendo: {0}. Relanzarlos con -overwrite borraria sus .log. " +
            "Usa -Solo para excluirlos, o -Pisar si de verdad quieres reiniciarlos.") -f ($choque -join ", "))
}

$cola = foreach ($a in $analisis) { foreach ($s in $Seeds) { [pscustomobject]@{ Tag = $a; Seed = $s } } }

$flagsBeagle = switch ($Beagle) {
    "cpu"  { @("-beagle_CPU", "-beagle_SSE") }
    "gpu"  { @("-beagle_GPU") }
    "java" { @("-java") }
    default { @() }
}

Write-Host "Lanzando $($cola.Count) cadenas ($($analisis.Count) analisis x $($Seeds.Count) semillas)"
Write-Host "  $Threads hilo(s) por cadena, maximo $MaxParalelo simultaneas, -Xmx$XmxPorProceso cada una"
Write-Host "  BEAGLE: $Beagle$(if ($flagsBeagle) { ' -> ' + ($flagsBeagle -join ' ') })"
Write-Host ""

$procs = @()
foreach ($job in $cola) {
    while (@($procs | Where-Object { -not $_.HasExited }).Count -ge $MaxParalelo) {
        Start-Sleep -Seconds 20
    }

    $dir = Join-Path $runsDir ("{0}\s{1}" -f $job.Tag, $job.Seed)
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    Copy-Item (Join-Path $xmlDir ($job.Tag + ".xml")) $dir -Force

    # Start-Process no entrecomilla los elementos de -ArgumentList, asi que una
    # ruta con espacios hay que entrecomillarla aqui o java la parte en dos.
    $argumentos = @(
        "-Xss256m", "-Xmx$XmxPorProceso", "-Dfile.encoding=UTF-8",
        "-cp", ('"{0}"' -f $launcher), "beast.pkgmgmt.launcher.BeastLauncher"
    )
    $argumentos += $flagsBeagle
    $argumentos += @("-threads", $Threads, "-overwrite", "-seed", $job.Seed,
                     ($job.Tag + ".xml"))

    $p = Start-Process -FilePath $java -ArgumentList $argumentos `
            -WorkingDirectory $dir -NoNewWindow -PassThru `
            -RedirectStandardOutput (Join-Path $dir "beast.out") `
            -RedirectStandardError  (Join-Path $dir "beast.err")

    $procs += $p
    Write-Host ("  [{0,5}] {1} semilla {2}" -f $p.Id, $job.Tag, $job.Seed)
}

Write-Host ""
Write-Host "Todas lanzadas. Esperando..."
$procs | Wait-Process

Write-Host ""
Write-Host "Resultados:"
foreach ($job in $cola) {
    $dir = Join-Path $runsDir ("{0}\s{1}" -f $job.Tag, $job.Seed)
    $out = Join-Path $dir "beast.out"
    $ok  = (Test-Path $out) -and (Select-String -Path $out -Pattern "^Done!" -Quiet)
    $estado = if ($ok) { "OK" } else { "FALLO" }
    Write-Host ("  {0,-6} {1} s{2}" -f $estado, $job.Tag, $job.Seed)
}

Write-Host ""
Write-Host "Para combinar las cadenas de un analisis (descartando 10% de burn-in):"
Write-Host "  logcombiner -log runs\<TAG>\s1\<TAG>.log -log runs\<TAG>\s2\<TAG>.log -log runs\<TAG>\s3\<TAG>.log -b 10 -o <TAG>_combined.log"
Write-Host "  logcombiner -log runs\<TAG>\s1\<TAG>.trees -log runs\<TAG>\s2\<TAG>.trees -log runs\<TAG>\s3\<TAG>.trees -b 10 -o <TAG>_combined.trees"
Write-Host "  treeannotator -burnin 0 -heights median <TAG>_combined.trees <TAG>_MCC.tree"
