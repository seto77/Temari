#Requires -Version 7
<#
.SYNOPSIS
    コード版フリート A/B — handout 配備判断用ベンチマーク。

.DESCRIPTION
    版 A (baseline): 指定コミット (既定 d688a35 = Phase 1 まで、P2-1 なし =
      配備中 handout のカーネル相当) をスクラッチパッドに clone/checkout。
      gen_production.jl・SCF キャッシュ・ハーネスは main の版をコピー
      (= handout の実運用ペアリングと同じ組合せ。S_GRID/e0_grid の同一性も保証)。
    版 B: main 現状 (読み取りのみ)。
    構成は本番フリート形態 8P×4T のみ。同一ジョブ表 (E1 と同じ 6ch×Rows) を
    A,B,A,B と交互に PassesEach パスずつ流し、rows/min 比と版間 F sha256
    (48 ジョブ全一致 = 配備ゲート) を測る。不一致時は F_hex から ULP 解剖し、
    両版の単発再実行 (probe_repro.jl) で「E8 型の一過性フリップ」か
    「実カーネル差」かを自動切り分けする。

    ⚠ 実行中はマシンの全 32 論理コアを ~40 分飽和させます (ワーカは
      BELOW_NORMAL)。既存 CPU 負荷が -MaxCpuPct を超えていれば開始を拒否。

.EXAMPLE
    pwsh -File tools\bench_e1\run_ab.ps1 -Smoke     # 配管検証 (数分)
    pwsh -File tools\bench_e1\run_ab.ps1            # 本計測 (~40 分)
#>
[CmdletBinding()]
param(
    [string]$Julia = "julia",
    [string]$JuliaChannel = "+1.11",   # 本番ピン 1.11.9
    [int]$Rows = 8,
    [int]$PassesEach = 2,              # 版ごとのパス数 (A,B 交互に実行)
    [string]$BaselineRef = "d688a35",  # 版 A のコミット
    [string]$CloneDir = "",            # 版 A の clone 先 (既定: scratchpad 相当の TEMP)
    [switch]$SkipPrep,                 # clone 準備済みならスキップ
    [switch]$Smoke,                    # 縮小配管検証 (1 ジョブ、2P×2T、各 1 パス)
    [switch]$Force,
    [int]$MaxCpuPct = 20,
    [int]$ConfigTimeoutMin = 180,
    [int]$StartupTimeoutMin = 20,
    [int]$CooldownSec = 45,
    [int]$StallKillMin = 10            # 260806Cl: ログ/出力停滞がこの分数を超えたら
                                       # 当該パスの全ワーカを kill (GC ウェッジ対策。
                                       # verA_p2 で実測: プロセス生存のまま計算停止)
)

Set-StrictMode -Version 3
$ErrorActionPreference = "Stop"

$BenchDir = $PSScriptRoot
$RepoRoot = (Resolve-Path (Join-Path $BenchDir "..\..")).Path
if (-not $CloneDir) { $CloneDir = Join-Path $env:TEMP "temari_ab_baseline" }
$JB = @(); if ($JuliaChannel) { $JB += $JuliaChannel }

function Write-Banner([string]$msg, [string]$color = "Cyan") {
    Write-Host ("=" * 72) -ForegroundColor $color
    Write-Host $msg -ForegroundColor $color
    Write-Host ("=" * 72) -ForegroundColor $color
}
function Format-CliArg([string]$a) {
    if ($a -match '[\s"]') { '"' + ($a -replace '"', '\"') + '"' } else { $a }
}
function Read-JsonFile([string]$p) { Get-Content $p -Raw | ConvertFrom-Json }
function HexToDouble([string]$h) {
    [BitConverter]::Int64BitsToDouble([BitConverter]::ToInt64(
        [BitConverter]::GetBytes([Convert]::ToUInt64($h, 16)), 0))
}
function HexToInt64([string]$h) {
    [BitConverter]::ToInt64([BitConverter]::GetBytes([Convert]::ToUInt64($h, 16)), 0)
}

# --------------------------------------------------------------------
# 0. 前提検査
# --------------------------------------------------------------------
$verLine = (& $Julia @JB --version) 2>&1 | Select-Object -First 1
if ($LASTEXITCODE -ne 0) { throw "julia の起動に失敗: $Julia $JB" }
$juliaVer = ($verLine -replace "julia version\s*", "").Trim()
Write-Host "Julia: $verLine"
if ($juliaVer -ne "1.11.9") {
    Write-Host "⚠ 本番ピン 1.11.9 と不一致 ($juliaVer)。A/B の配備判断材料としては無効" -ForegroundColor Yellow
}

$samples = foreach ($i in 1..3) {
    (Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage -Average).Average
    Start-Sleep -Seconds 2
}
$cpuLoad = [math]::Round(($samples | Measure-Object -Average).Average, 1)
Write-Host "既存 CPU 負荷: $cpuLoad % (3 サンプル平均)"
if ($cpuLoad -gt $MaxCpuPct -and -not $Force) {
    Write-Banner "中止: 既存 CPU 負荷 $cpuLoad % > $MaxCpuPct %。作者がマシンを使用中の可能性。`nアイドルになってから再実行してください (-Force で強行可)。" "Red"
    exit 2
}
Write-Banner ("⚠ コード版フリート A/B を開始します: 約 40 分、全 32 論理コアを飽和させます。`n" +
              "  ワーカは BELOW_NORMAL 優先度。計測中は他の作業をしないでください。") "Yellow"
Start-Sleep -Seconds 3

# --------------------------------------------------------------------
# 1. 版 A の準備 (clone + checkout + 配備相当ペアリング)
# --------------------------------------------------------------------
if (-not $SkipPrep) {
    if (-not (Test-Path (Join-Path $CloneDir ".git"))) {
        Write-Host "clone: $RepoRoot -> $CloneDir"
        git clone -q $RepoRoot $CloneDir
        if ($LASTEXITCODE -ne 0) { throw "git clone 失敗" }
    }
    git -C $CloneDir checkout -q --detach -f $BaselineRef
    if ($LASTEXITCODE -ne 0) { throw "git checkout $BaselineRef 失敗" }
    # 配備相当ペアリング: 運転系 (gen_production) とハーネスは main の版、
    # 計測対象カーネル (ionization.jl) だけが baseline に留まる。
    # これで S_GRID / e0_grid / HIGH_SETTINGS の定義が両版で同一になる。
    Copy-Item (Join-Path $RepoRoot "src\gen_production.jl") (Join-Path $CloneDir "src\") -Force
    New-Item -ItemType Directory -Force (Join-Path $CloneDir "tools\bench_e1") | Out-Null
    foreach ($f in @("bench_worker.jl", "probe_repro.jl")) {
        Copy-Item (Join-Path $BenchDir $f) (Join-Path $CloneDir "tools\bench_e1\") -Force
    }
    # SCF キャッシュを同一バイトで共有 (等温性 + 入力同一性)
    Copy-Item (Join-Path $RepoRoot "src\atom_cache_jl111_*.jls") (Join-Path $CloneDir "src\") -Force
}
$refA = (git -C $CloneDir rev-parse --short HEAD)
$refB = (git -C $RepoRoot rev-parse --short HEAD)
$ionizationClean = -not (git -C $CloneDir status --porcelain -- src/ionization.jl)
Write-Host "版 A: $CloneDir @ $refA (ionization.jl 無改変: $ionizationClean)"
Write-Host "版 B: $RepoRoot @ $refB"
if (-not $ionizationClean) { throw "版 A の src/ionization.jl が改変されている — baseline として無効" }

$WorkerA = Join-Path $CloneDir "tools\bench_e1\bench_worker.jl"
$WorkerB = Join-Path $BenchDir "bench_worker.jl"

# --------------------------------------------------------------------
# 2. 実行ディレクトリ + plan (両版で独立生成 → バイト一致を検証) + warmcache
# --------------------------------------------------------------------
$stamp  = Get-Date -Format "yyyyMMdd_HHmmss"
$RunDir = Join-Path $BenchDir ("results\" + $(if ($Smoke) { "absmoke_" } else { "ab_" }) + $stamp)
$DirA = Join-Path $RunDir "verA"; $DirB = Join-Path $RunDir "verB"
New-Item -ItemType Directory -Force $DirA, $DirB | Out-Null
Write-Host "`n実行ディレクトリ: $RunDir"

$planArgs = @("plan"); $planArgsB = @($DirB) + @("--rows", "$Rows"); $planArgsA = @($DirA) + @("--rows", "$Rows")
if ($Smoke) { $planArgsB += "--smoke"; $planArgsA += "--smoke" }
& $Julia @JB -t 1 $WorkerB @planArgs @planArgsB
if ($LASTEXITCODE -ne 0) { throw "plan (B) 失敗" }
& $Julia @JB -t 1 $WorkerA @planArgs @planArgsA
if ($LASTEXITCODE -ne 0) { throw "plan (A) 失敗" }
$pA = Get-Content (Join-Path $DirA "plan.json") -Raw
$pB = Get-Content (Join-Path $DirB "plan.json") -Raw
if ($pA -ne $pB) { throw "plan.json が版間で不一致 — S_GRID/e0_grid の定義が揃っていない。A/B は無効" }
Write-Host "plan: 両版で独立生成しバイト一致を確認 (同一ジョブ表)"
$nJobs = [int](ConvertFrom-Json $pB).n_jobs

Write-Host "`n--- SCF キャッシュのプリウォーム (両版、計測外) ---"
& $Julia @JB -t 8 --gcthreads=1 $WorkerA warmcache $DirA
if ($LASTEXITCODE -ne 0) { throw "warmcache (A) 失敗" }
& $Julia @JB -t 8 --gcthreads=1 $WorkerB warmcache $DirB
if ($LASTEXITCODE -ne 0) { throw "warmcache (B) 失敗" }

# --------------------------------------------------------------------
# 3. フリートパス実行 (run_e1.ps1 の実証済みロジックの A/B 版)
# --------------------------------------------------------------------
$fleetProcs   = $Smoke ? 2 : 8
$fleetThreads = $Smoke ? 2 : 4
if ($Smoke) { $PassesEach = 1 }

function Invoke-FleetPass([string]$ver, [int]$pass, [string]$suffix = "") {
    $worker = $ver -eq "A" ? $WorkerA : $WorkerB
    $verDir = $ver -eq "A" ? $DirA : $DirB
    $cfgName = "ver${ver}_${fleetProcs}p${fleetThreads}t_p$pass$suffix"
    $cfgDir = Join-Path $verDir "configs\$cfgName"
    New-Item -ItemType Directory -Force (Join-Path $cfgDir "logs") | Out-Null
    Write-Host "`n=== $cfgName (worker: $worker) ==="

    $procs = @()
    $spawn0 = Get-Date
    for ($w = 0; $w -lt $fleetProcs; $w++) {
        $argv = @()
        if ($JuliaChannel) { $argv += $JuliaChannel }
        $argv += @("-t", "$fleetThreads", $worker, "work", $cfgDir, "w$w")
        $p = Start-Process -FilePath $Julia `
            -ArgumentList ($argv | ForEach-Object { Format-CliArg $_ }) `
            -WorkingDirectory (Split-Path $worker) -WindowStyle Hidden -PassThru `
            -RedirectStandardOutput (Join-Path $cfgDir "logs\w$w.out.log") `
            -RedirectStandardError  (Join-Path $cfgDir "logs\w$w.err.log")
        $procs += $p
    }
    $deadline = (Get-Date).AddMinutes($StartupTimeoutMin)
    while ($true) {
        $ready = @(Get-ChildItem (Join-Path $cfgDir "ready") -Filter *.flag `
                   -ErrorAction SilentlyContinue).Count
        if ($ready -ge $fleetProcs) { break }
        $dead = @($procs | Where-Object { $_.HasExited })
        if ($dead.Count -gt 0) {
            $procs | Where-Object { -not $_.HasExited } |
                ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
            throw "${cfgName}: ワーカがバリア前に死亡 (exit=$($dead[0].ExitCode))。logs\ 参照"
        }
        if ((Get-Date) -gt $deadline) {
            $procs | ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
            throw "${cfgName}: 起動タイムアウト"
        }
        Start-Sleep -Milliseconds 500
    }
    $startupS = ((Get-Date) - $spawn0).TotalSeconds
    Write-Host ("  全 {0} ワーカ ready ({1:n1} s)。計測開始" -f $fleetProcs, $startupS)

    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    New-Item -ItemType File (Join-Path $cfgDir "go.flag") | Out-Null
    $tEnd = (Get-Date).AddMinutes($ConfigTimeoutMin)
    $timedOut = $false
    $wedged = $false          # 260806Cl: 停滞 watchdog による強制終了
    $wedgedAtExit = $false    # 260806Cl: 全ジョブ完了後にプロセスだけ残った型
    $lastMsg = Get-Date
    $lastActivity = Get-Date  # 260806Cl: out/done/logs の最新更新時刻
    while (@($procs | Where-Object { -not $_.HasExited }).Count -gt 0) {
        # 260806Cl: 停滞 watchdog。verA_p2 実測のウェッジ (プロセス生存・CPU 凍結・
        # ログ 45 分停滞) を N 分で検出して kill する。最重ジョブの無音区間は
        # 実測 ~3 分 (版 A) なので N=10 は十分な余裕
        $act = Get-ChildItem (Join-Path $cfgDir "out"), (Join-Path $cfgDir "done"), `
                             (Join-Path $cfgDir "logs") -File -ErrorAction SilentlyContinue |
               Measure-Object -Property LastWriteTime -Maximum
        if ($act.Count -gt 0 -and $act.Maximum -gt $lastActivity) { $lastActivity = $act.Maximum }
        $stalledMin = ((Get-Date) - $lastActivity).TotalMinutes
        if ($stalledMin -gt $StallKillMin) {
            $outN  = @(Get-ChildItem (Join-Path $cfgDir "out")  -Filter *.json -EA SilentlyContinue).Count
            $doneN = @(Get-ChildItem (Join-Path $cfgDir "done") -Filter *.json -EA SilentlyContinue).Count
            $alive = @($procs | Where-Object { -not $_.HasExited })
            Write-Host ("  ⚠ STALL-WEDGE: 出力停滞 {0:n1} 分 (jobs {1}/{2}, done {3}/{4}, 生存 {5} プロセス) — kill" `
                -f $stalledMin, $outN, $nJobs, $doneN, $fleetProcs, $alive.Count) -ForegroundColor Red
            $alive | ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
            if ($outN -eq $nJobs -and $doneN -eq $fleetProcs) {
                # 計算は全て完了済み・終了処理だけ固まった型 → 計測はワーカ時刻で救済可
                $wedgedAtExit = $true
                Write-Host "    (exit-wedge: 全ジョブ+全ワーカ統計は完了済み。wall はワーカ時刻から再構成)" -ForegroundColor Yellow
            } else {
                $wedged = $true
            }
            break
        }
        if ((Get-Date) -gt $tEnd) {
            $procs | Where-Object { -not $_.HasExited } |
                ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
            $timedOut = $true
            Write-Host "  ⚠ タイムアウト — 強制終了 (GC wedge の可能性)" -ForegroundColor Red
            break
        }
        if (((Get-Date) - $lastMsg).TotalSeconds -ge 20) {
            $doneN = @(Get-ChildItem (Join-Path $cfgDir "out") -Filter *.json `
                       -ErrorAction SilentlyContinue).Count
            Write-Host ("  ... {0}/{1} jobs  {2:n0} s" -f $doneN, $nJobs, $sw.Elapsed.TotalSeconds)
            $lastMsg = Get-Date
        }
        Start-Sleep -Seconds 2
    }
    $sw.Stop()
    $jobResults = @(Get-ChildItem (Join-Path $cfgDir "out") -Filter *.json `
                    -ErrorAction SilentlyContinue |
                    ForEach-Object { Read-JsonFile $_.FullName })
    $ok = (-not $timedOut) -and (-not $wedged) -and ($jobResults.Count -eq $nJobs)
    $wallS = [math]::Round($sw.Elapsed.TotalSeconds, 3)
    if ($wedgedAtExit) {
        # 260806Cl: ワーカ自身が記録した epoch (t_go〜t_exit) から真の span を再構成
        $ws = @(Get-ChildItem (Join-Path $cfgDir "done") -Filter *.json |
                ForEach-Object { Read-JsonFile $_.FullName })
        $wallS = [math]::Round((($ws | Measure-Object -Property t_exit -Maximum).Maximum -
                                ($ws | Measure-Object -Property t_go -Minimum).Minimum), 3)
        Write-Host ("    exit-wedge 救済 wall = {0:n1} s" -f $wallS) -ForegroundColor Yellow
    }
    $result = [ordered]@{
        name = $cfgName; version = $ver; pass = $pass
        engine_commit = ($ver -eq "A" ? $refA : $refB)
        procs = $fleetProcs; threads = $fleetThreads
        ok = $ok; timed_out = $timedOut
        wedged = $wedged; wedged_at_exit = $wedgedAtExit
        wall_s = $wallS; startup_s = [math]::Round($startupS, 3)
        rows_per_min = ($wallS -gt 0) ? [math]::Round(60.0 * $jobResults.Count / $wallS, 3) : 0
        busy_s_total = [math]::Round(($jobResults | Measure-Object -Property wall_s -Sum).Sum, 3)
        n_jobs_done = $jobResults.Count; n_jobs_expected = $nJobs
        jobs = $jobResults
    }
    $result | ConvertTo-Json -Depth 8 | Set-Content (Join-Path $RunDir "ab_$cfgName.json")
    Write-Host ("  {0}: wall={1:n1} s  rows/min={2}  jobs={3}/{4}  {5}" -f $cfgName, $wallS,
        $result.rows_per_min, $jobResults.Count, $nJobs, ($ok ? "OK" : "FAILED")) `
        -ForegroundColor ($ok ? "Green" : "Red")
    return [pscustomobject]$result
}

$sequence = @()
for ($p = 1; $p -le $PassesEach; $p++) { $sequence += @("A", "B") | ForEach-Object { @{ ver = $_; pass = $p } } }
$runs = @()
$cool = $Smoke ? 3 : $CooldownSec
for ($si = 0; $si -lt $sequence.Count; $si++) {
    $step = $sequence[$si]
    $r = Invoke-FleetPass $step.ver $step.pass
    $runs += $r
    if (-not $r.ok) {
        # 260806Cl: ウェッジ/失敗パスは 1 回だけ自動再走。2 連続で中断 (要人間判断)
        Write-Host "  ⚠ パス失敗 (wedged=$($r.wedged)) — cooldown 後に 1 回だけ自動再走" -ForegroundColor Yellow
        Start-Sleep -Seconds $cool
        $r2 = Invoke-FleetPass $step.ver $step.pass "_retry"
        $runs += $r2
        if (-not $r2.ok) {
            throw "同一パス ver$($step.ver) p$($step.pass) が 2 回連続で失敗/ウェッジ — 中断。results と logs を保全済み"
        }
    }
    if ($si -lt $sequence.Count - 1) {
        Write-Host "  (cooldown $cool s)"
        Start-Sleep -Seconds $cool
    }
}

# --------------------------------------------------------------------
# 4. 判定: スループット比 + 版間ビット同一性 (+ 不一致の自動切り分け)
# --------------------------------------------------------------------
$runsA = @($runs | Where-Object { $_.version -eq "A" -and $_.ok })
$runsB = @($runs | Where-Object { $_.version -eq "B" -and $_.ok })
# 260806Cl: (version, pass) ごとに ok な実行が 1 本あれば成立 (再走で回復した場合を含む)
$allOk = $true
foreach ($step in $sequence) {
    if (@($runs | Where-Object { $_.version -eq $step.ver -and $_.pass -eq $step.pass -and $_.ok }).Count -eq 0) {
        $allOk = $false
    }
}

# sha 集約: version -> jobid -> [sha per pass]
$shaMap = @{ A = @{}; B = @{} }
$hexMap = @{ A = @{}; B = @{} }      # jobid -> F_hex (代表 = 最初のパス)
foreach ($r in $runs) {
    foreach ($j in $r.jobs) {
        $v = $r.version
        if (-not $shaMap[$v].ContainsKey($j.id)) { $shaMap[$v][$j.id] = @(); $hexMap[$v][$j.id] = $j.F_hex }
        $shaMap[$v][$j.id] += $j.F_sha256
    }
}
$selfBad = @{}
foreach ($v in @("A", "B")) {
    $selfBad[$v] = @($shaMap[$v].Keys | Where-Object {
        @($shaMap[$v][$_] | Sort-Object -Unique).Count -gt 1 })
}
$crossBad = @($shaMap.A.Keys | Where-Object {
    $shaMap.B.ContainsKey($_) -and
    (@($shaMap.A[$_] + $shaMap.B[$_] | Sort-Object -Unique).Count -gt 1) })

# 不一致の ULP 解剖 + 単発再実行での切り分け (E8 型フリップ vs 実カーネル差)
$plan = ConvertFrom-Json $pB
$forensics = @()
foreach ($id in $crossBad) {
    $job = $plan.jobs | Where-Object { $_.id -eq $id } | Select-Object -First 1
    $fa = $hexMap.A[$id]; $fb = $hexMap.B[$id]
    $nd = 0; $maxUlp = [long]0; $maxRel = 0.0
    for ($i = 0; $i -lt $fa.Count; $i++) {
        if ($fa[$i] -eq $fb[$i]) { continue }
        $nd++
        $ulp = [math]::Abs((HexToInt64 $fa[$i]) - (HexToInt64 $fb[$i]))
        $da = HexToDouble $fa[$i]; $db = HexToDouble $fb[$i]
        $rel = [math]::Abs(($db - $da) / [math]::Max([math]::Abs($da), 1e-300))
        if ($ulp -gt $maxUlp) { $maxUlp = $ulp }
        if ($rel -gt $maxRel) { $maxRel = $rel }
    }
    Write-Host "`n[forensics] ${id}: 版間不一致 nodes=$nd maxULP=$maxUlp maxRel=$([math]::Round($maxRel,15)) — 単発再実行で切り分け" -ForegroundColor Yellow
    $reruns = @{}
    foreach ($v in @("A", "B")) {
        $probe = $v -eq "A" ? (Join-Path $CloneDir "tools\bench_e1\probe_repro.jl") : (Join-Path $BenchDir "probe_repro.jl")
        $outp = Join-Path $RunDir "recheck_${id}_ver$v.json"
        & $Julia @JB -t 4 $probe "$($job.z)" $job.tag "$($job.e0_keV)" "high" "2" $outp | Out-Null
        $reruns[$v] = ($LASTEXITCODE -eq 0) ? @((Read-JsonFile $outp).reps | ForEach-Object sha256) : @("RERUN_FAILED")
    }
    $rerunIdentical = (@($reruns.A + $reruns.B | Sort-Object -Unique).Count -eq 1)
    $forensics += [pscustomobject]@{
        id = $id; nodes_diff = $nd; max_ulp = $maxUlp; max_rel = $maxRel
        rerun_identical = $rerunIdentical
        classification = $rerunIdentical ?
            "transient-flip (E8 型: フリート負荷時のみ、単発では両版一致)" :
            "REAL-KERNEL-DIFF (単発でも両版不一致 — 配備ブロッカー)"
    }
}

# ゲート判定
$realDiffs = @($forensics | Where-Object { -not $_.rerun_identical })
$gate = if (-not $allOk) { "NO-GO (実行失敗パスあり)" }
        elseif ($crossBad.Count -eq 0) { "GO ($($shaMap.A.Keys.Count) ジョブ全パス版間ビット同一)" }
        elseif ($realDiffs.Count -eq 0) { "GO 条件付き (不一致 $($crossBad.Count) 件は全て一過性フリップと分類、単発では両版一致)" }
        else { "NO-GO (実カーネル差 $($realDiffs.Count) 件)" }

$wallA = ($runsA | Measure-Object -Property wall_s -Average).Average
$wallB = ($runsB | Measure-Object -Property wall_s -Average).Average
$ratio = ($wallB -gt 0) ? [math]::Round($wallA / $wallB, 3) : 0

$md = [System.Collections.Generic.List[string]]::new()
$md.Add("# コード版フリート A/B ($stamp)$($Smoke ? ' — SMOKE (配管検証のみ)' : '')")
$md.Add("")
$md.Add("- 版 A (baseline) = $refA (Phase 1 まで、P2-1 なし = handout カーネル相当) + main の gen_production/ハーネス")
$md.Add("- 版 B = main @ $refB")
$md.Add("- 構成: ${fleetProcs}P×${fleetThreads}T、ジョブ $nJobs、パス A,B 交互 ×$PassesEach、Julia $juliaVer、開始時負荷 $cpuLoad %")
$md.Add("")
$md.Add("| pass | version | wall [s] | rows/min | jobs | status |")
$md.Add("|---|---|---:|---:|---|---|")
foreach ($r in $runs) {
    $st = $r.ok ? ($r.wedged_at_exit ? "OK (exit-wedge 救済)" : "OK") :
          ($r.wedged ? "WEDGED" : ($r.timed_out ? "TIMEOUT" : "FAILED"))
    $md.Add("| $($r.name) | $($r.version)@$($r.engine_commit) | $($r.wall_s) | $($r.rows_per_min) | $($r.n_jobs_done)/$($r.n_jobs_expected) | $st |")
}
$md.Add("")
if ($runsA.Count -and $runsB.Count) {
    $md.Add("**スループット比 B/A = $ratio** (wall 平均 A=$([math]::Round($wallA,1)) s / B=$([math]::Round($wallB,1)) s)")
}
$md.Add("")
$md.Add("版内自己一致: A 不一致 $($selfBad.A.Count) 件, B 不一致 $($selfBad.B.Count) 件$($selfBad.A.Count + $selfBad.B.Count -gt 0 ? " — $(@($selfBad.A + $selfBad.B) -join ', ')" : '')")
$md.Add("版間ビット同一性: 不一致 $($crossBad.Count) / $($shaMap.A.Keys.Count) ジョブ")
foreach ($f in $forensics) {
    $md.Add("- $($f.id): nodes=$($f.nodes_diff) maxULP=$($f.max_ulp) maxRel=$($f.max_rel) → $($f.classification)")
}
$md.Add("")
$md.Add("## 判定: $gate")
$md -join "`n" | Set-Content (Join-Path $RunDir "summary_ab.md")

[ordered]@{
    stamp = $stamp; verA = "$refA"; verB = "$refB"; julia = $juliaVer
    n_jobs = $nJobs; passes_each = $PassesEach
    wall_avg_A = $wallA; wall_avg_B = $wallB; throughput_ratio_BA = $ratio
    self_mismatch_A = $selfBad.A; self_mismatch_B = $selfBad.B
    cross_mismatch = $crossBad; forensics = $forensics; gate = $gate
} | ConvertTo-Json -Depth 5 | Set-Content (Join-Path $RunDir "ab_verdict.json")

Write-Host ""
Write-Banner "A/B 完了: $RunDir"
Get-Content (Join-Path $RunDir "summary_ab.md") | Write-Host
exit ($allOk ? 0 : 1)
