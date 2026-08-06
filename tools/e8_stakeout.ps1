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
    # "Z:ch:E0" のカンマ区切り。ワーカへラウンドロビン配分 (16 ワーカ 2 標的 → 8+8)
    [string]$Targets = "6:K:275,38:L3:40",
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
$targetList = @()
foreach ($t in ($Targets -split ",")) {
    $f = $t -split ":"
    $targetList += ,@([int]$f[0], $f[1].ToUpper(), [double]$f[2])
}
"run root: $root"
"workers : $Workers x t$ThreadsPer (BELOW_NORMAL)   targets: $Targets   limits: $MaxPasses passes/worker, $MaxHours h"

# ---- ワーカ起動 (ラウンドロビンで標的を配分、優先度 BELOW_NORMAL) ----
$procs = @()
$wdKey = @{}
for ($w = 1; $w -le $Workers; $w++) {
    $tgt = $targetList[($w - 1) % $targetList.Count]
    $key = "Z$($tgt[0])_$($tgt[1])_E0$($tgt[2])"
    $wd = Join-Path $root ("w{0:D2}_{1}" -f $w, $key)
    $wdKey[$wd] = $key
    $p = Start-Process julia -PassThru -WindowStyle Hidden `
        -RedirectStandardOutput (Join-Path $root ("w{0:D2}.out.log" -f $w)) `
        -RedirectStandardError  (Join-Path $root ("w{0:D2}.err.log" -f $w)) `
        -ArgumentList @($JuliaChannel, "-t", "$ThreadsPer", "--gcthreads=1",
                        $workerJl, $wd, "$($tgt[0])", $tgt[1], "$($tgt[2])", "$MaxPasses")
    try { $p.PriorityClass = [System.Diagnostics.ProcessPriorityClass]::BelowNormal } catch {}
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

# 相違ノードの分布 (ULP = ビットパターンの UInt64 距離、rel = |Δ|/max(|a|,|b|)) を
# 印字する。相殺増幅 (小さい値のノードだけ rel が跳ねる) の検証用。
function Show-DiffStats([string[]]$hexA, [string[]]$hexB, [string]$label) {
    $n = [Math]::Min($hexA.Count, $hexB.Count)
    $rows = @()
    $maxUlp = [uint64]0; $iMaxUlp = -1
    $maxRel = 0.0; $iMaxRel = -1
    for ($i = 0; $i -lt $n; $i++) {
        if ($hexA[$i] -eq $hexB[$i]) { continue }
        $ua = [System.Convert]::ToUInt64($hexA[$i], 16)
        $ub = [System.Convert]::ToUInt64($hexB[$i], 16)
        $ulp = if ($ua -ge $ub) { $ua - $ub } else { $ub - $ua }
        $va = [System.BitConverter]::UInt64BitsToDouble($ua)
        $vb = [System.BitConverter]::UInt64BitsToDouble($ub)
        $den = [Math]::Max([Math]::Abs($va), [Math]::Abs($vb))
        $rel = if ($den -gt 0) { [Math]::Abs($va - $vb) / $den } else { 0.0 }
        if ($ulp -gt $maxUlp) { $maxUlp = $ulp; $iMaxUlp = $i }
        if ($rel -gt $maxRel) { $maxRel = $rel; $iMaxRel = $i }
        $rows += ("    {0}#{1}  ref={2} bad={3}  |dULP|={4}  rel={5:E2}  (ref={6:E6})" -f
                  $label, $i, $hexA[$i], $hexB[$i], $ulp, $rel, $va)
    }
    if ($rows.Count -eq 0) { "  $label : 相違なし"; return }
    "  $label 相違 $($rows.Count)/$n 点   maxULP=$maxUlp @#$iMaxUlp   maxRel=$('{0:E2}' -f $maxRel) @#$iMaxRel"
    $rows | Select-Object -First 60
    if ($rows.Count -gt 60) { "    ... (残り $($rows.Count - 60) 点は省略)" }
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
        $alignSame = ($jr.dNde_ptr_mod64 -eq $jb.dNde_ptr_mod64) -and
                     ($jr.dNde_ptr_mod4096 -eq $jb.dNde_ptr_mod4096) -and
                     ($jr.we_ptr_mod64 -eq $jb.we_ptr_mod64) -and
                     ($jr.N_ptr_mod64 -eq $jb.N_ptr_mod64)
        "    align ref/bad: dNde mod64 $($jr.dNde_ptr_mod64)/$($jb.dNde_ptr_mod64), mod4096 $($jr.dNde_ptr_mod4096)/$($jb.dNde_ptr_mod4096); we mod64 $($jr.we_ptr_mod64)/$($jb.we_ptr_mod64); N mod64 $($jr.N_ptr_mod64)/$($jb.N_ptr_mod64)  → $(if ($alignSame) { '同一' } else { '★相違' })"
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
        if ($nDiff -gt 0) { Show-DiffStats $Nr $Nb "N" }
        if ($diffSlices.Count -gt 0) {
            "  判定 (a): ε スライス相違 → ノード内部起因 (eps_worker 内のワークスペース/BLAS 疑い)"
            "    相違 ε ノード (1-based, $($diffSlices.Count)/$($jr.ne)): $($diffSlices -join ',')"
            "    縮約後 N の相違: $nDiff/$nK 点"
        }
        elseif ($nDiff -gt 0) {
            "  判定 (b): 全 ε スライス一致・縮約後 N のみ相違 ($nDiff/$nK 点)"
            "    → 縮約の文脈依存丸め (N = dNde' * we の BLAS gemv 'T')。index-order 固定の自前ループ化で修正可能"
            if (-not $alignSame) {
                "    細分 (b1): ポインタ整列も相違 → 整列依存 peeling 仮説の強い証拠 (GC 配置揺れ)"
            }
            else {
                "    細分 (b2): 整列は同一なのに N 相違 → BLAS スレッド分割 (部分和結合位置) / その他"
            }
        }
        else {
            "  判定 (c): スライスも N も一致 → 相違は N→F 後段 (F=N./N[1] は決定的なはず)"
            "    → メモリ破壊系 (既知の Windows GC 問題の同族) を疑う。要エスカレーション"
        }
    }
}

# ---- 監視ループ (標的ごとに参照パスを持つ) ----
$deadline = (Get-Date).AddHours($MaxHours)
$refPass = @{}
$refF = @{}
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
        $wd = Split-Path $pass -Parent
        $key = $wdKey[$wd]
        if ($null -eq $key) { $key = "?" }
        $f = Get-Content (Join-Path $pass "F.hex")
        if (-not $refF.ContainsKey($key)) { $refF[$key] = $f; $refPass[$key] = $pass; continue }
        $nCompared++
        $rF = $refF[$key]
        $diffIdx = @()
        for ($i = 0; $i -lt [Math]::Min($rF.Count, $f.Count); $i++) {
            if ($rF[$i] -ne $f[$i]) { $diffIdx += $i }
        }
        if ($rF.Count -ne $f.Count) { $diffIdx += -1 }
        if ($diffIdx.Count -gt 0) {
            $mismatch = $true
            ""
            "================ MISMATCH DETECTED ================"
            "target: $key"
            "ref : $($refPass[$key])"
            "bad : $pass"
            Show-DiffStats $rF $f "F"
            $n0r = Get-Content (Join-Path $refPass[$key] "N0.hex")
            $n0b = Get-Content (Join-Path $pass "N0.hex")
            "N0: $(if ($n0r -eq $n0b) { 'IDENTICAL' } else { "DIFFER ($n0r vs $n0b)" })"
            ""
            Compare-Sidecars $refPass[$key] $pass
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
