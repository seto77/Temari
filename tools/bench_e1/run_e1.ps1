#Requires -Version 7
<#
.SYNOPSIS
    実験 E1: プロセス並列 vs スレッド並列の公平ベンチマークのドライバ。

.DESCRIPTION
    docs/notes/speedup_audit_2026-08-05.md P4-1 の再実験。同一のジョブ表 (plan.json、
    本番 e0_grid 由来の (Z, tag, E0) 行) を 4 つのプロセス×スレッド構成で流し、
    バリア同期後の壁時計時間を比較する。詳細は README.md。

    ⚠ このスクリプトは計測中マシンの全 32 論理コアを飽和させます。
      有効な数字を得るにはマシンが完全にアイドルであること (起動時に検査)。

.EXAMPLE
    pwsh -File tools\bench_e1\run_e1.ps1 -Smoke        # 配管の煙試験 (数分・数字は無意味)
    pwsh -File tools\bench_e1\run_e1.ps1               # 本計測 (既定 -Rows 8、目安 0.5-1.5 h)
    pwsh -File tools\bench_e1\run_e1.ps1 -Rows 2       # 1/4 スケールの下見
#>
[CmdletBinding()]
param(
    [string]$Julia = "julia",
    [string]$JuliaChannel = "+1.11",   # juliaup チャネル (本番ピン 1.11.9)。"" で PATH の julia
    [int]$Rows = 8,                    # チャネルあたりの E0 行数 (6 チャネル × Rows = ジョブ数)
    [switch]$Smoke,                    # 煙試験モード: 極小ジョブ + 縮小 2 構成で配管のみ検証
    [switch]$Force,                    # CPU アイドル検査をスキップ (非推奨)
    [int]$MaxCpuPct = 20,              # これを超える既存負荷があれば中止
    [int]$ConfigTimeoutMin = 180,      # 1 構成の上限 (GC wedge 対策の watchdog)
    [int]$StartupTimeoutMin = 20,      # 全ワーカが ready になるまでの上限
    [int]$CooldownSec = 45,            # 構成間の休止 (熱の持ち越し軽減)
    [switch]$Reverse,                  # 構成を逆順 (c4→c1) で実行 (熱・順序バイアスの対照)
    [int]$Passes = 1,                  # 各構成の実行回数。2 以上で sha チェックは
                                       # 「構成内パス間の自己一致」に切り替わる
                                       # (エンジンに実行順依存がある場合の運用モード)
    [string]$OutRoot = ""
)

Set-StrictMode -Version 3
$ErrorActionPreference = "Stop"

$BenchDir = $PSScriptRoot
$RepoRoot = (Resolve-Path (Join-Path $BenchDir "..\..")).Path
$SrcDir   = Join-Path $RepoRoot "src"
$Worker   = Join-Path $BenchDir "bench_worker.jl"
if (-not $OutRoot) { $OutRoot = Join-Path $BenchDir "results" }

$JuliaBase = @()
if ($JuliaChannel) { $JuliaBase += $JuliaChannel }

function Write-Banner([string]$msg, [string]$color = "Cyan") {
    Write-Host ("=" * 72) -ForegroundColor $color
    Write-Host $msg -ForegroundColor $color
    Write-Host ("=" * 72) -ForegroundColor $color
}

# --------------------------------------------------------------------
# 0. 前提検査: Julia バージョン / CPU アイドル
# --------------------------------------------------------------------
$verLine = (& $Julia @JuliaBase --version) 2>&1 | Select-Object -First 1
if ($LASTEXITCODE -ne 0) { throw "julia の起動に失敗: $Julia $JuliaBase" }
$juliaVer = ($verLine -replace "julia version\s*", "").Trim()
Write-Host "Julia: $verLine  (`"$Julia $($JuliaBase -join ' ')`")"
if ($juliaVer -ne "1.11.9") {
    Write-Host "⚠ 本番ピンは 1.11.9 (audit の掟 / atom_cache_jl111_*)。現在 $juliaVer。" -ForegroundColor Yellow
    Write-Host "  ・SCF キャッシュはバージョン別 (jl<maj><min>) なので warmcache が別途作り直します" -ForegroundColor Yellow
    Write-Host "  ・1.12/Windows には既知の並列 GC segfault があります (gen_production.jl 冒頭)" -ForegroundColor Yellow
    Write-Host "  ・この計測値は本番構成の判断材料として弱くなります (README 参照)" -ForegroundColor Yellow
}

# CPU アイドル検査 (Win32_Processor.LoadPercentage はロケール非依存)
$samples = foreach ($i in 1..3) {
    (Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage -Average).Average
    Start-Sleep -Seconds 2
}
$cpuLoad = [math]::Round(($samples | Measure-Object -Average).Average, 1)
Write-Host "既存 CPU 負荷: $cpuLoad % (3 サンプル平均)"
if ($cpuLoad -gt $MaxCpuPct) {
    if ($Force) {
        Write-Host "⚠ 負荷 $cpuLoad % > $MaxCpuPct % だが -Force 指定のため続行 (数字の妥当性は保証されない)" -ForegroundColor Yellow
    } else {
        Write-Banner "中止: 既存 CPU 負荷 $cpuLoad % > $MaxCpuPct %。`nベンチはアイドルなマシンでのみ有効です。他の作業を止めて再実行してください (-Force で強行可)。" "Red"
        exit 2
    }
}

Write-Banner ("⚠ 警告: この実行はマシンの全 32 論理コアを長時間飽和させます。`n" +
              "  ワーカは BELOW_NORMAL 優先度ですが、計測を汚すので実行中は他の作業をしないこと。`n" +
              "  中断は Ctrl+C (起動済み julia が残った場合は Stop-Process -Name julia)。") "Yellow"
Start-Sleep -Seconds 3

# --------------------------------------------------------------------
# 1. 実行ディレクトリ + plan + warmcache
# --------------------------------------------------------------------
$stamp  = Get-Date -Format "yyyyMMdd_HHmmss"
$RunDir = Join-Path $OutRoot ($(if ($Smoke) { "smoke_" } else { "e1_" }) + $stamp)
New-Item -ItemType Directory -Force $RunDir | Out-Null
Write-Host "`n実行ディレクトリ: $RunDir"

$planArgs = @("plan", $RunDir, "--rows", "$Rows")
if ($Smoke) { $planArgs += "--smoke" }
& $Julia @JuliaBase -t 1 $Worker @planArgs
if ($LASTEXITCODE -ne 0) { throw "plan 失敗 (exit=$LASTEXITCODE)" }
$plan  = Get-Content (Join-Path $RunDir "plan.json") -Raw | ConvertFrom-Json
$nJobs = [int]$plan.n_jobs

Write-Host "`n--- SCF キャッシュのプリウォーム (計測外。全構成が同じ温度で走るための前提) ---"
$twc = [System.Diagnostics.Stopwatch]::StartNew()
& $Julia @JuliaBase -t 8 --gcthreads=1 $Worker warmcache $RunDir
if ($LASTEXITCODE -ne 0) { throw "warmcache 失敗 (exit=$LASTEXITCODE)" }
$twc.Stop()
Write-Host ("warmcache: {0:n1} s" -f $twc.Elapsed.TotalSeconds)

# --------------------------------------------------------------------
# 2. 構成の定義 (E1 の 4 構成 / smoke は縮小 2 構成)
# --------------------------------------------------------------------
if ($Smoke) {
    $Configs = @(
        [pscustomobject]@{ name = "smokeA_2p2t";     procs = 2; threads = 2;  gcthreads = $null },
        [pscustomobject]@{ name = "smokeB_1p4t_gc1"; procs = 1; threads = 4;  gcthreads = 1 }
    )
} else {
    $Configs = @(
        [pscustomobject]@{ name = "c1_8p4t";     procs = 8;  threads = 4;  gcthreads = $null },
        [pscustomobject]@{ name = "c2_1p32t";    procs = 1;  threads = 32; gcthreads = $null },
        [pscustomobject]@{ name = "c3_16p2t";    procs = 16; threads = 2;  gcthreads = $null },
        [pscustomobject]@{ name = "c4_8p4t_gc1"; procs = 8;  threads = 4;  gcthreads = 1 }
    )
}
if ($Reverse) { [array]::Reverse($Configs); Write-Host "構成は逆順で実行 ($(($Configs | ForEach-Object name) -join ' -> '))" }

# 実行メタデータ
$cpu = Get-CimInstance Win32_Processor | Select-Object -First 1
$commit = (& git -C $RepoRoot rev-parse --short HEAD 2>$null)
[ordered]@{
    experiment = "E1"; stamp = $stamp; smoke = [bool]$Smoke
    rows = $Rows; n_jobs = $nJobs
    julia = "$Julia $($JuliaBase -join ' ')"; julia_version = $juliaVer
    repo_commit = "$commit"; cpu = $cpu.Name.Trim()
    logical_processors = $cpu.NumberOfLogicalProcessors
    cpu_load_at_start_pct = $cpuLoad
    configs = $Configs | ForEach-Object { "$($_.name)" }
} | ConvertTo-Json -Depth 4 | Set-Content (Join-Path $RunDir "meta.json")

# --------------------------------------------------------------------
# 3. 構成の実行
# --------------------------------------------------------------------
function Format-CliArg([string]$a) {
    if ($a -match '[\s"]') { '"' + ($a -replace '"', '\"') + '"' } else { $a }
}

function Invoke-BenchConfig($cfg, [string]$cfgName) {
    $cfgDir = Join-Path $RunDir "configs\$cfgName"
    New-Item -ItemType Directory -Force (Join-Path $cfgDir "logs") | Out-Null
    Write-Host "`n=== config ${cfgName}: $($cfg.procs) proc x $($cfg.threads) threads" `
        ($null -ne $cfg.gcthreads ? "--gcthreads=$($cfg.gcthreads)" : "(gcthreads 既定)") "==="

    $procs = @()
    $spawn0 = Get-Date
    for ($w = 0; $w -lt $cfg.procs; $w++) {
        $argv = @()
        if ($JuliaChannel) { $argv += $JuliaChannel }
        $argv += @("-t", "$($cfg.threads)")
        if ($null -ne $cfg.gcthreads) { $argv += "--gcthreads=$($cfg.gcthreads)" }
        $argv += @($Worker, "work", $cfgDir, "w$w")
        $p = Start-Process -FilePath $Julia `
            -ArgumentList ($argv | ForEach-Object { Format-CliArg $_ }) `
            -WorkingDirectory $SrcDir -WindowStyle Hidden -PassThru `
            -RedirectStandardOutput (Join-Path $cfgDir "logs\w$w.out.log") `
            -RedirectStandardError  (Join-Path $cfgDir "logs\w$w.err.log")
        $procs += $p
    }

    # バリア: 全ワーカが include + JIT ウォームアップを終えるのを待つ
    $deadline = (Get-Date).AddMinutes($StartupTimeoutMin)
    while ($true) {
        $ready = @(Get-ChildItem (Join-Path $cfgDir "ready") -Filter *.flag `
                   -ErrorAction SilentlyContinue).Count
        if ($ready -ge $cfg.procs) { break }
        $dead = @($procs | Where-Object { $_.HasExited })
        if ($dead.Count -gt 0) {
            $procs | Where-Object { -not $_.HasExited } |
                ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
            throw "config ${cfgName}: ワーカがバリア前に死亡 (exit=$($dead[0].ExitCode))。logs\ を参照"
        }
        if ((Get-Date) -gt $deadline) {
            $procs | ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
            throw "config ${cfgName}: 起動タイムアウト ($StartupTimeoutMin min)"
        }
        Start-Sleep -Milliseconds 500
    }
    $startupS = ((Get-Date) - $spawn0).TotalSeconds
    Write-Host ("  全 {0} ワーカ ready ({1:n1} s)。計測開始" -f $cfg.procs, $startupS)

    # go → 全プロセス終了 までが計測区間
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    New-Item -ItemType File (Join-Path $cfgDir "go.flag") | Out-Null
    $tEnd = (Get-Date).AddMinutes($ConfigTimeoutMin)
    $timedOut = $false
    $lastMsg = Get-Date
    while (@($procs | Where-Object { -not $_.HasExited }).Count -gt 0) {
        if ((Get-Date) -gt $tEnd) {
            $procs | Where-Object { -not $_.HasExited } |
                ForEach-Object { Stop-Process -Id $_.Id -Force -ErrorAction SilentlyContinue }
            $timedOut = $true
            Write-Host "  ⚠ タイムアウト ($ConfigTimeoutMin min) — GC wedge の可能性 (IMPORT.md)。強制終了" -ForegroundColor Red
            break
        }
        if (((Get-Date) - $lastMsg).TotalSeconds -ge 15) {
            $doneN = @(Get-ChildItem (Join-Path $cfgDir "out") -Filter *.json `
                       -ErrorAction SilentlyContinue).Count
            Write-Host ("  ... {0}/{1} jobs  {2:n0} s" -f $doneN, $nJobs, $sw.Elapsed.TotalSeconds)
            $lastMsg = Get-Date
        }
        Start-Sleep -Seconds 2
    }
    $sw.Stop()
    $wallS = $sw.Elapsed.TotalSeconds

    # 結果の収集
    $jobResults = @(Get-ChildItem (Join-Path $cfgDir "out") -Filter *.json `
                    -ErrorAction SilentlyContinue |
                    ForEach-Object { Get-Content $_.FullName -Raw | ConvertFrom-Json })
    $workerStats = @(Get-ChildItem (Join-Path $cfgDir "done") -Filter *.json `
                     -ErrorAction SilentlyContinue |
                     ForEach-Object { Get-Content $_.FullName -Raw | ConvertFrom-Json })
    $workSpanS = $null
    if ($jobResults.Count -gt 0) {
        $t0 = ($jobResults | Measure-Object -Property t_start -Minimum).Minimum
        $t1 = ($jobResults | Measure-Object -Property t_end -Maximum).Maximum
        $workSpanS = [math]::Round($t1 - $t0, 3)
    }
    $busyS = [math]::Round(($jobResults | Measure-Object -Property wall_s -Sum).Sum, 3)
    $ok = (-not $timedOut) -and ($jobResults.Count -eq $nJobs)

    $result = [ordered]@{
        config = $cfgName; base_config = $cfg.name
        procs = $cfg.procs; threads = $cfg.threads
        gcthreads = $cfg.gcthreads
        ok = $ok; timed_out = $timedOut
        wall_s = [math]::Round($wallS, 3)          # go.flag → 最後のワーカ終了 (主指標)
        startup_s = [math]::Round($startupS, 3)    # spawn → 全員 ready (参考)
        work_span_s = $workSpanS                   # ワーカ時刻での最初の job 開始 → 最後の job 終了
        busy_s_total = $busyS                      # Σ 行時間 (レーン内直列換算の総仕事)
        n_jobs_done = $jobResults.Count; n_jobs_expected = $nJobs
        jobs = $jobResults; workers = $workerStats
    }
    $result | ConvertTo-Json -Depth 8 | Set-Content (Join-Path $RunDir "config_$cfgName.json")
    Write-Host ("  {0}: wall={1:n1} s  jobs={2}/{3}  {4}" -f $cfgName, $wallS,
        $jobResults.Count, $nJobs, ($ok ? "OK" : "FAILED")) `
        -ForegroundColor ($ok ? "Green" : "Red")
    return $ok
}

$allOk = $true
$totalRuns = $Configs.Count * $Passes
$runIdx = 0
for ($ci = 0; $ci -lt $Configs.Count; $ci++) {
    for ($pass = 1; $pass -le $Passes; $pass++) {
        $runIdx++
        $cfgName = ($Passes -gt 1) ? "$($Configs[$ci].name)_p$pass" : $Configs[$ci].name
        $ok = Invoke-BenchConfig $Configs[$ci] $cfgName
        if (-not $ok) { $allOk = $false }
        if ($runIdx -lt $totalRuns) {
            $cool = $Smoke ? 3 : $CooldownSec
            Write-Host "  (cooldown $cool s)"
            Start-Sleep -Seconds $cool
        }
    }
}

# --------------------------------------------------------------------
# 4. サマリ (summary.md): 表 + ビット同一性クロスチェック
# --------------------------------------------------------------------
$cfgFiles = @(Get-ChildItem (Join-Path $RunDir "config_*.json") | Sort-Object Name)
$results  = @($cfgFiles | ForEach-Object { Get-Content $_.FullName -Raw | ConvertFrom-Json })
$baseline = $results | Where-Object { $_.base_config -in @("c1_8p4t", "smokeA_2p2t") } |
            Sort-Object config | Select-Object -First 1

$md = [System.Collections.Generic.List[string]]::new()
$md.Add("# E1 benchmark summary ($stamp)$($Smoke ? ' — SMOKE (数字は無意味、配管検証のみ)' : '')")
$md.Add("")
$md.Add("- CPU: $($cpu.Name.Trim()) ($($cpu.NumberOfLogicalProcessors) logical)")
$md.Add("- Julia: $juliaVer (`$Julia $($JuliaBase -join ' ')`)  / 本番ピン 1.11.9 $($juliaVer -eq '1.11.9' ? 'OK' : '**不一致 ⚠**')")
$md.Add("- repo commit: $commit / jobs: $nJobs (rows/ch=$Rows) / 開始時 CPU 負荷: $cpuLoad %")
$md.Add("- 計測区間 = バリア (全ワーカ JIT ウォームアップ済み) → 最終ワーカ終了。起動・JIT・SCF ウォームは計測外")
$md.Add("")
$md.Add("| config | P x T | gcthreads | wall [s] | rows/min | speedup vs $($baseline ? $baseline.config : '-') | jobs | Σrow [s] | status |")
$md.Add("|---|---|---|---:|---:|---:|---|---:|---|")
foreach ($r in $results) {
    $rpm = ($r.wall_s -gt 0) ? [math]::Round(60.0 * $r.n_jobs_done / $r.wall_s, 2) : 0
    $sp  = ($baseline -and $r.wall_s -gt 0) ? [math]::Round($baseline.wall_s / $r.wall_s, 3) : "-"
    $gct = ($null -ne $r.gcthreads) ? "$($r.gcthreads)" : "default"
    $md.Add("| $($r.config) | $($r.procs) x $($r.threads) | $gct | $([math]::Round($r.wall_s,1)) | $rpm | $sp | $($r.n_jobs_done)/$($r.n_jobs_expected) | $([math]::Round($r.busy_s_total,1)) | $($r.ok ? 'OK' : 'FAILED') |")
}
$md.Add("")

# ビット同一性: 同じジョブの F の sha256 比較。
# Passes==1: 構成間クロスチェック (不一致 = exit 3)。
# Passes>=2: 主チェックは「構成内のパス間自己一致」に切り替え (エンジンに実行順
#            依存の非決定性がある場合の運用モード)。クロスは参考情報として残す。
$hashByJob = @{}
foreach ($r in $results) {
    foreach ($j in $r.jobs) {
        if (-not $hashByJob.ContainsKey($j.id)) { $hashByJob[$j.id] = @{} }
        $hashByJob[$j.id][$r.config] = $j.F_sha256
    }
}
$mismatch = @($hashByJob.Keys | Where-Object {
    @($hashByJob[$_].Values | Sort-Object -Unique).Count -gt 1 })
if ($mismatch.Count -eq 0) {
    $md.Add("F の sha256 クロスチェック: 全 $($hashByJob.Count) ジョブが全構成 (全パス) でビット同一 — OK")
} else {
    $tagx = ($Passes -ge 2) ? " (参考情報。Passes>=2 では自己一致が主チェック)" : ""
    $md.Add("**⚠ F の sha256 構成間不一致: $($mismatch.Count) ジョブ**$tagx — $($mismatch -join ', ')")
}
if ($Passes -ge 2) {
    $md.Add("")
    $md.Add("構成内自己一致 (パス間):")
    foreach ($grp in ($results | Group-Object base_config)) {
        $byJob = @{}
        foreach ($r in $grp.Group) {
            foreach ($j in $r.jobs) {
                if (-not $byJob.ContainsKey($j.id)) { $byJob[$j.id] = @() }
                $byJob[$j.id] += $j.F_sha256
            }
        }
        $bad = @($byJob.Keys | Where-Object {
            @($byJob[$_] | Sort-Object -Unique).Count -gt 1 })
        if ($bad.Count -eq 0) {
            $md.Add("- $($grp.Name): $($byJob.Count) ジョブ自己一致 OK")
        } else {
            $md.Add("- **$($grp.Name): $($bad.Count) ジョブが自己不一致** — $($bad -join ', ')")
        }
    }
}
$md.Add("")
$md.Add("注意: 構成は直列実行 (順序 = 表の順$($Reverse ? '、-Reverse 指定で逆順' : ''))。熱・電力状態の持ち越しは cooldown $CooldownSec s で軽減しているが完全ではない。結論を出す前に順序を変えた再実行を 1 回推奨。")
$md -join "`n" | Set-Content (Join-Path $RunDir "summary.md")

Write-Host ""
Write-Banner "完了: $RunDir`n  summary.md / config_*.json を参照"
Get-Content (Join-Path $RunDir "summary.md") | Write-Host
if ($mismatch.Count -gt 0 -and $Passes -lt 2) {
    Write-Host "⚠ ビット同一性の不一致あり" -ForegroundColor Red; exit 3
}
exit ($allOk ? 0 : 1)
