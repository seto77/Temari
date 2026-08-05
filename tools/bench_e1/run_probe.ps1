#Requires -Version 7
<#
E1 sha256 不一致の切り分けマトリクス (probe_repro.jl を逐次実行)。
対象 2 ジョブ × {-t 4 ×3 プロセス, -t 4 gcthreads=1, -t 1 ×2, -t 2, -t 32}
× 各プロセス内 2 reps (HIGH) + QUICK 連打 (-t 4 / -t 32 × 10 reps)。
#>
[CmdletBinding()]
param(
    [string]$Julia = "julia",
    [string]$JuliaChannel = "+1.11",
    [string]$OutDir = ""
)
Set-StrictMode -Version 3
$ErrorActionPreference = "Stop"

$BenchDir = $PSScriptRoot
$Probe = Join-Path $BenchDir "probe_repro.jl"
if (-not $OutDir) {
    $OutDir = Join-Path $BenchDir ("results\probe_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
}
New-Item -ItemType Directory -Force $OutDir | Out-Null
$JB = @(); if ($JuliaChannel) { $JB += $JuliaChannel }

$jobs = @(
    @{ label = "Z006_K_E0030";  z = 6;  tag = "K";  e0 = "30.0" },
    @{ label = "Z038_L3_E0060"; z = 38; tag = "L3"; e0 = "60.0" }
)
# name, threads, gcthreads(null=既定), mode, reps
$matrix = @(
    @{ n = "t4_a";   t = 4;  g = $null; m = "high";  r = 2 },
    @{ n = "t4_b";   t = 4;  g = $null; m = "high";  r = 2 },
    @{ n = "t4_c";   t = 4;  g = $null; m = "high";  r = 2 },
    @{ n = "t4_gc1"; t = 4;  g = 1;     m = "high";  r = 2 },
    @{ n = "t1_a";   t = 1;  g = $null; m = "high";  r = 2 },
    @{ n = "t1_b";   t = 1;  g = $null; m = "high";  r = 1 },
    @{ n = "t2_a";   t = 2;  g = $null; m = "high";  r = 2 },
    @{ n = "t32_a";  t = 32; g = $null; m = "high";  r = 2 },
    @{ n = "q_t4";   t = 4;  g = $null; m = "quick"; r = 10 },
    @{ n = "q_t32";  t = 32; g = $null; m = "quick"; r = 10 }
)

$sw = [System.Diagnostics.Stopwatch]::StartNew()
foreach ($j in $jobs) {
    foreach ($mrow in $matrix) {
        $out = Join-Path $OutDir "$($j.label)__$($mrow.m)_$($mrow.n).json"
        $argv = @()
        $argv += @("-t", "$($mrow.t)")
        if ($null -ne $mrow.g) { $argv += "--gcthreads=$($mrow.g)" }
        $argv += @($Probe, "$($j.z)", $j.tag, $j.e0, $mrow.m, "$($mrow.r)", $out)
        Write-Host ("[{0:n0}s] {1} {2}" -f $sw.Elapsed.TotalSeconds, $j.label, $mrow.n)
        & $Julia @JB @argv
        if ($LASTEXITCODE -ne 0) { Write-Host "  ⚠ FAILED (exit=$LASTEXITCODE)" -ForegroundColor Red }
    }
}
Write-Host ("done in {0:n0}s -> {1}" -f $sw.Elapsed.TotalSeconds, $OutDir)
