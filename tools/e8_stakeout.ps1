# e8_stakeout.ps1 — E8 待ち伏せドライバ (260806Cl 追加)
#
# 16P×2T 相当の負荷を張りつつ、全ワーカで同一行 (既定 Z=6 K E0=275) を反復計算。
# パス間で F.hex を突合し、不一致を検出したら STOP を置いて全ワーカを止め、
# 当該ペアのサイドカーを自動突合して「どの ε ノードから違うか / 縮約後から
# 違うか」を印字して終了する。上限 = MaxPasses (ワーカごと) と MaxHours。
#
#   powershell -File tools\e8_stakeout.ps1            # 既定 16P×2T, 2h
#   powershell -File tools\e8_stakeout.ps1 -Workers 8 -MaxHours 0.5
#
# ⚠ 実行はフリート A/B 完了後、コーディネータの許可を得てから。
param(
    [int]$Workers = 16,
    [int]$ThreadsPer = 2,
    [int]$Z = 6,
    [string]$Channel = "K",
    [double]$E0 = 275.0,
    [int]$MaxPasses = 200,
    [double]$MaxHours = 2.0,
    [string]$JuliaChannel = "+1.11",
    [string]$OutRoot = "$env:TEMP\e8_stakeout"
)
$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$workerJl = Join-Path $repo "tools\e8_worker.jl"
$root = Join-Path $OutRoot ("e8_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Force $root | Out-Null
"run root: $root"
"workers : $Workers x t$ThreadsPer   target: Z=$Z $Channel E0=$E0   limits: $MaxPasses passes/worker, $MaxHours h"

# ---- ワーカ起動 ----
$procs = @()
for ($w = 1; $w -le $Workers; $w++) {
    $wd = Join-Path $root ("w{0:D2}" -f $w)
    $p = Start-Process julia -PassThru -WindowStyle Hidden `
        -RedirectStandardOutput (Join-Path $root ("w{0:D2}.out.log" -f $w)) `
        -RedirectStandardError  (Join-Path $root ("w{0:D2}.err.log" -f $w)) `
        -ArgumentList @($JuliaChannel, "-t", "$ThreadsPer", "--gcthreads=1",
                        $workerJl, $wd, "$Z", $Channel, "$E0", "$MaxPasses")
    $procs += $p
}
"launched $($procs.Count) workers (pids: $($procs.Id -join ','))"

function Stop-AllWorkers([object[]]$procs, [string]$root) {
    New-Item -ItemType File -Force (Join-Path $root "STOP") | Out-Null
    "STOP placed; waiting up to 10 min for current passes to finish..."
    $t0 = Get-Date
    while ((($procs | Where-Object { -not $_.HasExited }).Count -gt 0) -and
           ((Get-Date) - $t0).TotalMinutes -lt 10) { Start-Sleep -Seconds 10 }
    foreach ($p in $procs) {
        if (-not $p.HasExited) { try { Stop-Process -Id $p.Id -Force } catch {} }
    }
}

function Get-UlpDiff([string]$hexA, [string]$hexB) {
    $a = [System.Convert]::ToUInt64($hexA, 16)
    $b = [System.Convert]::ToUInt64($hexB, 16)
    if ($a -ge $b) { return $a - $b } else { return $b - $a }
}

function Compare-Sidecars([string]$refPass, [string]$badPass) {
    $sr = Get-ChildItem $refPass -Filter "e8_pid*.json" | Sort-Object Name
    $sb = Get-ChildItem $badPass -Filter "e8_pid*.json" | Sort-Object Name
    if ($sr.Count -eq 0 -or $sb.Count -eq 0) {
        "!! サイドカーが見つからない (ref=$($sr.Count), bad=$($sb.Count)) — 計装ビルドで実行したか確認"
        return
    }
    if ($sr.Count -ne $sb.Count) { "!! サイドカー数不一致: $($sr.Count) vs $($sb.Count)" }
    $n = [Math]::Min($sr.Count, $sb.Count)
    for ($k = 0; $k -lt $n; $k++) {
        $jr = Get-Content $sr[$k].FullName -Raw | ConvertFrom-Json
        $jb = Get-Content $sb[$k].FullName -Raw | ConvertFrom-Json
        "--- sidecar[$k]: ne=$($jr.ne) nK=$($jr.nK) jt=$($jr.julia_threads)/$($jb.julia_threads) bt=$($jr.blas_threads)/$($jb.blas_threads)"
        if ($jr.eps_sha -ne $jb.eps_sha) { "  ★eps (求積ノード) 自体が相違 → 上流 (gl01/eigvals=LAPACK stev) 起因" }
        if ($jr.we_sha -ne $jb.we_sha)   { "  ★we (求積重み) 自体が相違 → 上流 (gl01/eigvals) 起因" }
        $diffSlices = @()
        for ($i = 0; $i -lt $jr.slice_sha.Count; $i++) {
            if ($jr.slice_sha[$i] -ne $jb.slice_sha[$i]) { $diffSlices += ($i + 1) }
        }
        $nDiff = 0
        $Nr = $jr.N_hex -split ","
        $Nb = $jb.N_hex -split ","
        $nK = [Math]::Min($Nr.Count, $Nb.Count)
        for ($i = 0; $i -lt $nK; $i++) { if ($Nr[$i] -ne $Nb[$i]) { $nDiff++ } }
        if ($diffSlices.Count -gt 0) {
            "  判定 (a): ε スライス相違 → ノード内部起因 (eps_worker 内のワークスペース/BLAS 疑い)"
            "    相違 ε ノード (1-based, $($diffSlices.Count)/$($jr.ne)): $($diffSlices -join ',')"
            "    縮約後 N の相違: $nDiff/$nK 点"
        }
        elseif ($nDiff -gt 0) {
            "  判定 (b): 全 ε スライス一致・縮約後 N のみ相違 ($nDiff/$nK 点)"
            "    → 縮約順序起因 (N = dNde' * we の BLAS gemv 'T')。index-order 固定の自前ループ化で修正可能"
        }
        else {
            "  判定 (c): スライスも N も一致 → 相違は N→F 後段 (F=N./N[1] は決定的なはず)"
            "    → メモリ破壊系 (既知の Windows GC 問題の同族) を疑う。要エスカレーション"
        }
    }
}

# ---- 監視ループ ----
$deadline = (Get-Date).AddHours($MaxHours)
$refPass = $null
$refF = $null
$seen = @{}
$nCompared = 0
$mismatch = $false
while ($true) {
    Start-Sleep -Seconds 15
    $doneMarks = Get-ChildItem $root -Recurse -Filter DONE -ErrorAction SilentlyContinue
    foreach ($dm in $doneMarks) {
        $pass = $dm.DirectoryName
        if ($seen.ContainsKey($pass)) { continue }
        $seen[$pass] = $true
        $f = Get-Content (Join-Path $pass "F.hex")
        if ($null -eq $refF) { $refF = $f; $refPass = $pass; continue }
        $nCompared++
        $diffIdx = @()
        for ($i = 0; $i -lt [Math]::Min($refF.Count, $f.Count); $i++) {
            if ($refF[$i] -ne $f[$i]) { $diffIdx += $i }
        }
        if ($refF.Count -ne $f.Count) { $diffIdx += -1 }
        if ($diffIdx.Count -gt 0) {
            $mismatch = $true
            ""
            "================ MISMATCH DETECTED ================"
            "ref : $refPass"
            "bad : $pass"
            "F 相違点 (0-based): $($diffIdx.Count)/$($refF.Count)"
            foreach ($i in ($diffIdx | Select-Object -First 8)) {
                if ($i -ge 0) {
                    $ulp = Get-UlpDiff $refF[$i] $f[$i]
                    "  s#$i  ref=$($refF[$i])  bad=$($f[$i])  |dULP|=$ulp"
                }
            }
            $n0r = Get-Content (Join-Path $refPass "N0.hex")
            $n0b = Get-Content (Join-Path $pass "N0.hex")
            "N0: $(if ($n0r -eq $n0b) { 'IDENTICAL' } else { "DIFFER ($n0r vs $n0b)" })"
            ""
            Compare-Sidecars $refPass $pass
            break
        }
    }
    if ($mismatch) { break }
    $alive = ($procs | Where-Object { -not $_.HasExited }).Count
    "$(Get-Date -Format HH:mm:ss)  passes compared: $nCompared  workers alive: $alive"
    if ($alive -eq 0) { "all workers done — no mismatch in $nCompared comparisons"; break }
    if ((Get-Date) -gt $deadline) { "time limit ($MaxHours h) reached — no mismatch in $nCompared comparisons"; break }
}
Stop-AllWorkers $procs $root
"run root kept for inspection: $root"
exit $(if ($mismatch) { 2 } else { 0 })
