<#
jobq bootstrap (tools/jobq/PROTOCOL.md section 10) - registers this PC as a jobq worker host.

  powershell -NoProfile -ExecutionPolicy Bypass -File \\10.31.108.5\jobq\setup\bootstrap.ps1
      [-Slots N] [-Threads T] [-Remove] [-DryRun] [-Root R] [-Spool S] [-Local L] [-User DOMAIN\name]

Normally started by double-clicking register.cmd / unregister.cmd at the share root; those elevate
themselves and hand this script -Root and -User.

-Root  the share ROOT a person browses: register.cmd, unregister.cmd, README.txt, setup\, code\, spool\.
-Spool where machines write (default <Root>\spool). Split on purpose: the root has to stay readable
       to a person looking for the file to double-click, so nothing machine-written may land there.
-User  the account the workers must run as, as seen *before* UAC. juliaup and the Credential Manager
       are per-user, so if UAC elevated a different administrator this script stops rather than
       register tasks whose profile has no julia and no saved share credentials. With -Remove it is
       only a warning: unregistering is machine-wide administrator work and must stay possible on a
       PC whose worker account no longer exists.
-Local the local working root (default C:\jobq); a non-default value is passed to the workers as
       JOBQ_LOCAL in the task action.

Run elevated while logged on as the worker account. Asks for that account's password once (Task
Scheduler "run whether logged on or not"). Re-run = update (idempotent; WORKER_ID is kept): worker
tasks for slots < Slots are re-registered in place (a running worker keeps running, with the
worker.conf it read at start); slots >= Slots are stopped (whole process tree) and unregistered.
-Remove = stop (process tree) + unregister all jobq tasks, retired_utc in the host record.
-DryRun = print every action, change nothing, no prompt.
Windows PowerShell 5.1 and pwsh 7. ASCII only, LF line endings.
#>
[CmdletBinding()]
param(
  [int]$Slots = 0,
  [int]$Threads = 0,
  [switch]$Remove,
  [switch]$DryRun,
  [string]$Root = '\\10.31.108.5\jobq',
  [string]$Spool = '',
  [string]$Local = 'C:\jobq',
  [string]$User = ''
)
Set-StrictMode -Version 2
$ErrorActionPreference = 'Stop'
# A trailing '\' (tab completion adds one) would turn the closing quote of -Root "...\" in the nastest command line into \".
$Root = $Root.TrimEnd('\', '/')
$Local = $Local.TrimEnd('\', '/')
$Spool = $Spool.TrimEnd('\', '/')
if (-not $Root -or -not $Local) { throw '-Root and -Local must not be empty' }
if (-not $Spool) { $Spool = Join-Path $Root 'spool' }        # PROTOCOL section 1.2: SPOOL defaults to ROOT/spool
$Utf8NoBom = New-Object System.Text.UTF8Encoding($false)
$WorkerTaskGlob = 'jobq-worker-s*'
$NasTestTask = 'jobq-nastest'

# ---------------------------------------------------------------- helpers
function Say([string]$m) { Write-Host "[bootstrap] $m" }
function Step([string]$what, [scriptblock]$do) {      # every mutating action goes through here
  if ($DryRun) { Write-Host "[dry-run] would: $what" } else { Say $what; & $do }
}
function UtcNow { (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ') }
function Prop($obj, [string]$name, $default) {       # PSCustomObject property with default
  if ($null -ne $obj -and $obj.PSObject.Properties[$name]) { return $obj.$name } else { return $default }
}
function Write-LfFile([string]$path, [string]$text) { # UTF-8 without BOM, LF only
  [IO.File]::WriteAllText($path, ($text -replace "`r`n", "`n"), $Utf8NoBom)
}
function Write-JsonAtomic([string]$path, $obj) {      # tmp + rename (PROTOCOL section 1)
  $dir = Split-Path -Parent $path
  if (-not (Test-Path -LiteralPath $dir)) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }
  $tmp = Join-Path $dir ('.' + (Split-Path -Leaf $path) + '.tmp')
  Write-LfFile $tmp (($obj | ConvertTo-Json -Depth 8) + "`n")
  Move-Item -LiteralPath $tmp -Destination $path -Force
}
function ConvertTo-MsysPath([string]$p) {             # C:\jobq -> /c/jobq, \\host\share -> //host/share
  $p = $p.TrimEnd('\', '/')
  if ($p -match '^([A-Za-z]):(.*)$') { return '/' + $Matches[1].ToLower() + ($Matches[2] -replace '\\', '/') }
  return ($p -replace '\\', '/')
}
function Invoke-Native([string]$exe, [string[]]$argv) { # stderr merged; exit code in $script:NativeExit
  $old = $ErrorActionPreference; $ErrorActionPreference = 'Continue'
  try { $out = (& $exe @argv 2>&1 | ForEach-Object { "$_" }) -join "`n"; $script:NativeExit = $LASTEXITCODE }
  finally { $ErrorActionPreference = $old }
  return $out
}
function Invoke-NativeStreaming([string]$exe, [string[]]$argv) {  # same, but prints as it goes
  # A winget install on a fresh PC takes minutes. Buffering its output until it returns makes the
  # cmd window register.cmd opened look hung, and an operator who kills it gets a half-installed PC.
  $old = $ErrorActionPreference; $ErrorActionPreference = 'Continue'
  try { & $exe @argv 2>&1 | ForEach-Object { Write-Host "    $_" }; $script:NativeExit = $LASTEXITCODE }
  finally { $ErrorActionPreference = $old }
}
function Update-SessionPath {                          # pick up PATH entries added by winget
  $env:Path = [Environment]::GetEnvironmentVariable('Path', 'Machine') + ';' + [Environment]::GetEnvironmentVariable('Path', 'User')
}
function Find-GitBash {                                # default path -> registry -> git.exe on PATH
  $cands = @('C:\Program Files\Git\bin\bash.exe')
  foreach ($k in 'HKLM:\SOFTWARE\GitForWindows', 'HKCU:\SOFTWARE\GitForWindows') {
    $ip = Get-ItemProperty -Path $k -Name InstallPath -ErrorAction SilentlyContinue
    if ($ip) { $cands += (Join-Path $ip.InstallPath 'bin\bash.exe') }
  }
  $git = Get-Command git.exe -ErrorAction SilentlyContinue
  if ($git) { $cands += (Join-Path (Split-Path -Parent (Split-Path -Parent $git.Source)) 'bin\bash.exe') }
  foreach ($c in $cands) { if (Test-Path -LiteralPath $c) { return $c } }
  return $null
}
function Install-IfMissing([string]$id, [scriptblock]$present) {
  if (& $present) { Say "$id already installed"; return }
  Step "winget install --id $id -e --accept-source-agreements --accept-package-agreements" {
    Invoke-NativeStreaming 'winget' @('install', '--id', $id, '-e', '--accept-source-agreements', '--accept-package-agreements')
    if ($script:NativeExit -ne 0) { throw "winget install $id failed (exit $script:NativeExit)" }
    Update-SessionPath
  }
}
function Install-JuliaChannel([string]$ver) {          # juliaup add only when the channel is absent
  $have = $false
  if (Get-Command juliaup.exe -ErrorAction SilentlyContinue) {
    $have = (Invoke-Native 'juliaup' @('status')) -match ('(?m)^\s*\*?\s*' + [regex]::Escape($ver) + '\s')
  }
  if ($have) { Say "julia $ver already installed (juliaup status)"; return }
  Step "juliaup add $ver" {
    Write-Host (Invoke-Native 'juliaup' @('add', $ver))
    if ($script:NativeExit -ne 0) { throw "juliaup add $ver failed (exit $script:NativeExit)" }
  }
}
function Get-ExistingWorkerId {
  if (-not (Test-Path -LiteralPath $confPath)) { return $null }
  foreach ($line in [IO.File]::ReadAllLines($confPath)) {
    if ($line -match '^WORKER_ID=([a-z0-9][a-z0-9-]{0,40})\s*$') { return $Matches[1] }
  }
  return $null
}
function New-WorkerId {                                # <hostname lower, [^a-z0-9] -> '-'>-<8 hex>
  $h = ($hostName.ToLower() -replace '[^a-z0-9]', '-').Trim('-')
  if ($h.Length -gt 32) { $h = $h.Substring(0, 32) }
  if ($h -notmatch '^[a-z0-9]') { $h = 'host' }
  return $h + '-' + [guid]::NewGuid().ToString('N').Substring(0, 8)
}
function Get-WorkerSlot([string]$taskName) {           # jobq-worker-s<k> -> k, anything else -> -1
  if ($taskName -match '^jobq-worker-s(\d+)$') { return [int]$Matches[1] } else { return -1 }
}
function ConvertTo-ShWord([string]$s) {               # one bash word: plain chars as-is, otherwise single-quoted
  if ($s -match '^[A-Za-z0-9_./:=-]+$') { return $s } else { return "'" + ($s -replace "'", "'\''") + "'" }
}
function Get-MsysProcTable {                           # ps -W = "PID PPID PGID WINPID TTY UID STIME COMMAND"
  # MSYS ids survive exec, Win32 ParentProcessId does not, so this is the only view of the whole worker tree.
  $exe = Find-GitBash
  if (-not $exe) { return $null }
  $txt = Invoke-Native $exe @('-c', 'ps -W')
  if ($script:NativeExit -ne 0 -or -not $txt) { return $null }
  $rows = @()
  foreach ($line in ($txt -split "`n")) {
    if ($line -match '^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s') {
      $rows += [pscustomobject]@{ MPid = [int]$Matches[1]; MPpid = [int]$Matches[2]; MPgid = [int]$Matches[3]; Winpid = [uint32]$Matches[4] }
    }
  }
  if ($rows.Count -eq 0) { return $null }
  return $rows
}
function Stop-WorkerTree($task) {
  # Task Scheduler's Stop ends only the process it started (the bin\bash.exe launcher); worker.sh and julia keep
  # running, the slot's status file keeps ticking and the reaper never reaps. Killing the rest is not a Win32 tree
  # walk: bash -lc "<worker.sh> <k>" EXECs the script, so (measured on this host) the surviving worker carries
  # "<worker.sh> <k>" *without* -lc and its ParentProcessId points at the exec'd-away launcher, and a child started
  # with "cmd &" (julia) is Win32-orphaned too. Only the MSYS ids survive exec. So: match the command payload ->
  # ps -W -> the whole MSYS process group -> taskkill /T /F per WINPID -> wait <= 30 s, re-killing anything that
  # still carries the payload (PROTOCOL sections 5.3 / 11; same semantics as kill_tree in tools/lane_watchdog.sh).
  $actArgs = $null
  try { $actArgs = [string]$task.Actions[0].Arguments } catch { $actArgs = $null }
  if (-not $actArgs) { Say "  $($task.TaskName): no action arguments - process tree not searched"; return }
  $payload = $actArgs.Trim()
  if ($payload -match '^-[A-Za-z]*c\s+"(.*)"$') { $payload = $Matches[1] }                              # -lc "<cmd>" -> <cmd>
  $payload = $payload -replace '^(?:[A-Za-z_][A-Za-z0-9_]*=(?:''[^'']*''|"[^"]*"|\S*)\s+)+', ''         # VAR=... prefixes become environment on exec
  if (-not $payload) { Say "  $($task.TaskName): no command in [$actArgs] - process tree not searched"; return }
  $re = '(^|[\s"])' + [regex]::Escape($payload) + '"?\s*$'   # slot-specific: "...worker.sh 1" cannot match "...worker.sh 10"
  $procs = @(Get-CimInstance Win32_Process | Select-Object ProcessId, ParentProcessId, Name, CommandLine, CreationDate)
  $seeds = @($procs | Where-Object { $_.Name -eq 'bash.exe' -and $_.CommandLine -and $_.CommandLine -match $re })   # the launcher and every exec/fork of it (PROTOCOL section 10 step 4)
  if ($seeds.Count -eq 0) { Say "  $($task.TaskName): no process carries [$payload]"; return }
  $tree = @{}                                          # winpid -> CreationDate (everything we intend to kill)
  $byPid = @{}; foreach ($p in $procs) { $byPid[[uint32]$p.ProcessId] = $p }
  $add = { param($id) if (-not $tree.ContainsKey($id) -and $id -ne [uint32]$PID -and $byPid.ContainsKey($id)) { $tree[$id] = $byPid[$id].CreationDate } }
  foreach ($s in $seeds) { & $add ([uint32]$s.ProcessId) }
  $rows = Get-MsysProcTable
  if ($rows) {                                         # whole MSYS process group of every seed (survives exec)
    $byWin = @{}; foreach ($r in $rows) { $byWin[$r.Winpid] = $r }
    $pgids = @{}
    foreach ($id in @($tree.Keys)) { if ($byWin.ContainsKey($id) -and $byWin[$id].MPgid -gt 0) { $pgids[$byWin[$id].MPgid] = $true } }
    foreach ($r in $rows) { if ($pgids.ContainsKey($r.MPgid)) { & $add $r.Winpid } }
  } else { Say "  $($task.TaskName): ps -W unavailable - falling back to the Win32 parent chain (exec'd children may survive)" }
  $pending = @($tree.Keys)                             # plus Win32 descendants (native children of julia etc.)
  while ($pending.Count -gt 0) {
    $id = $pending[0]; $pending = @($pending | Select-Object -Skip 1)
    if (-not $byPid.ContainsKey($id)) { continue }
    foreach ($c in @($procs | Where-Object { $_.ParentProcessId -eq $id -and $_.CreationDate -ge $byPid[$id].CreationDate })) {
      $cid = [uint32]$c.ProcessId
      if (-not $tree.ContainsKey($cid)) { & $add $cid; $pending += $cid }
    }
  }
  Say "  $($task.TaskName): taskkill /T /F on $($tree.Count) process(es) carrying or grouped with [$payload]: $((@($tree.Keys) | ForEach-Object { "$($byPid[$_].Name):$_" }) -join ' ')"
  foreach ($id in @($tree.Keys)) {
    Invoke-Native 'taskkill.exe' @('/PID', "$id", '/T', '/F') | Out-Null
    if ($script:NativeExit -ne 0 -and $script:NativeExit -ne 128) { Say "  taskkill $id exit $script:NativeExit" }   # 128 = already gone
  }
  $alive = @()
  for ($i = 0; $i -lt 30; $i++) {
    $now = @(Get-CimInstance Win32_Process | Select-Object ProcessId, Name, CommandLine, CreationDate)
    $alive = @($now | Where-Object { $tree.ContainsKey([uint32]$_.ProcessId) -and $tree[[uint32]$_.ProcessId] -eq $_.CreationDate })
    $fresh = @($now | Where-Object { $_.Name -eq 'bash.exe' -and $_.CommandLine -and $_.CommandLine -match $re -and -not $tree.ContainsKey([uint32]$_.ProcessId) })
    if ($alive.Count -eq 0 -and $fresh.Count -eq 0) { break }
    foreach ($f in $fresh) {                           # forked after the snapshot
      $tree[[uint32]$f.ProcessId] = $f.CreationDate
      Invoke-Native 'taskkill.exe' @('/PID', "$($f.ProcessId)", '/T', '/F') | Out-Null
    }
    Start-Sleep -Seconds 1
  }
  if ($alive.Count -gt 0) { Say "  WARNING: $($task.TaskName): still alive after 30 s: $(($alive | ForEach-Object { "$($_.Name):$($_.ProcessId)" }) -join ' ')" }
}
function Remove-JobqTasks([string[]]$names) {          # stop -> kill the worker's process tree -> unregister. The slot's claim is left to the reaper.
  # The reaper notices within claim_timeout because the slot's status file stops ticking (PROTOCOL section 7).
  foreach ($n in $names) {
    foreach ($t in @(Get-ScheduledTask -TaskName $n -ErrorAction SilentlyContinue)) {
      $isWorker = (Get-WorkerSlot $t.TaskName) -ge 0
      $verb = if ($isWorker) { 'stop task + kill its process tree + unregister' } else { 'stop + unregister task' }
      Step "$verb $($t.TaskName)" {
        if ($t.State -eq 'Running') { Stop-ScheduledTask -TaskName $t.TaskName -ErrorAction SilentlyContinue }
        if ($isWorker) { Stop-WorkerTree $t }
        Unregister-ScheduledTask -TaskName $t.TaskName -Confirm:$false
      }
    }
  }
}
function Register-JobqTask([string]$name, $action, $trigger, $settings, [string]$desc) {
  # -User/-Password => logon type Password (runs whether the user is logged on or not)
  $p = @{ TaskName = $name; Action = $action; Settings = $settings; User = $userName; Password = $plainPw; Description = $desc; Force = $true }
  if ($trigger) { $p.Trigger = $trigger }
  Register-ScheduledTask @p | Out-Null
}

# ---------------------------------------------------------------- identity / paths
$me = [Security.Principal.WindowsIdentity]::GetCurrent()
$isAdmin = (New-Object Security.Principal.WindowsPrincipal($me)).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
if (-not $isAdmin) {
  if ($DryRun) { Say 'WARNING: not elevated (tolerated for -DryRun)' } else { throw 'Run from an elevated (administrator) PowerShell.' }
}
$userName = $me.Name
$hostName = $env:COMPUTERNAME
# -User carries the account that was logged on before UAC (register.cmd passes %USERDOMAIN%\%USERNAME%).
# UAC may elevate as a different administrator, and juliaup, the julia channels and the Credential
# Manager entries for the share all live in the *user* profile: registering the worker tasks for an
# account that has none of them produces a host that fails every job, so stop here instead.
if ($User) {
  $want = $User.Trim()
  if ($want -and -not $want.Equals($userName, [StringComparison]::OrdinalIgnoreCase)) {
    # -Remove only stops tasks, kills their process trees and unregisters them - all machine-wide
    # administrator work that touches nothing in the worker's profile. Blocking it would leave a PC
    # whose worker account is gone with no way to clean itself up by double-click, so warn and go on.
    if ($Remove) {
      Say "WARNING: unregister.cmd was elevated as $userName but the worker account is $want - removing anyway (removal is machine-wide)"
    } else {
      # Plain Write-Host, not throw: this is read by a person in the cmd window register.cmd leaves
      # open at `pause`, and a PowerShell error record buries the instruction in a stack trace.
      Write-Host ''
      Write-Host 'FAIL: wrong account.'
      Write-Host "  register.cmd was elevated as $userName but the worker account is $want."
      Write-Host "  Log on as $want, or make $want a local administrator, then run register.cmd again."
      Write-Host '  (juliaup and the Credential Manager are per-user: tasks registered for the wrong'
      Write-Host '   profile would find no julia and no saved share credentials.)'
      Write-Host ''
      exit 1
    }
  } else { Say "user check ok: running as $userName" }
}
$setupSrc = Join-Path $Root 'setup'
$setupDst = Join-Path $Local 'setup'
$confPath = Join-Path $Local 'worker.conf'
$hostsDir = Join-Path $Spool 'hosts'                   # PROTOCOL section 1.2: the ledger is machine-written -> SPOOL
$spoolMsys = ConvertTo-MsysPath $Spool
$localMsys = ConvertTo-MsysPath $Local                 # /c/jobq = worker.sh's built-in default
Say "root=$Root spool=$Spool local=$Local user=$userName host=$hostName dry-run=$DryRun"

# ---------------------------------------------------------------- -Remove
if ($Remove) {
  Remove-JobqTasks @($WorkerTaskGlob, $NasTestTask)
  $wid = Get-ExistingWorkerId
  if (-not $wid) { Say "no WORKER_ID in $confPath - tasks removed, no host record to retire"; exit 0 }
  $p = Join-Path $hostsDir "$wid.json"
  # Everything from here on is bookkeeping on the share, and the machine-local work is already done:
  # the tasks are stopped, their process trees killed and the tasks unregistered. Neither an
  # unreachable share nor a corrupt ledger can undo that - and both are *likely* on the PC an operator
  # is unregistering (it lost access to the NAS, or it is being decommissioned). $ErrorActionPreference
  # is 'Stop', so letting either throw ends the script here: the raw error record and unregister.cmd's
  # "[FAIL] ... exit code 1" would report a removal that in fact succeeded (PROTOCOL section 10.1: that
  # cmd window is the whole PASS/FAIL report), and the stale slot status files below would be skipped.
  # Measured 2026-08-21: Write-JsonAtomic's New-Item on \\localhost\<nonexistent share> and
  # ConvertFrom-Json on a truncated record both terminate. So: warn, keep going, exit 0 - the same
  # pattern the NAS-test failure path uses further down.
  $rec = [ordered]@{}
  try {
    if (Test-Path -LiteralPath $p) {
      foreach ($pp in (Get-Content -LiteralPath $p -Raw | ConvertFrom-Json).PSObject.Properties) { $rec[$pp.Name] = $pp.Value }
    }
  } catch {
    Say "WARNING: cannot read $p ($($_.Exception.Message)) - writing a minimal record instead"
    $rec = [ordered]@{}
  }
  if ($rec.Count -eq 0) { $rec['worker_id'] = $wid; $rec['hostname'] = $hostName }
  $rec['retired_utc'] = UtcNow
  Step "write retired_utc into $p (tmp + rename)" {
    try { Write-JsonAtomic $p $rec }
    catch { Say "WARNING: could not record retired_utc in $p ($($_.Exception.Message)) - this PC IS unregistered (tasks stopped and removed); only the ledger on the share is left unchanged" }
  }
  # The slots are gone, so their status files can only mislead 'queuectl hosts' and the reaper's
  # liveness check (a missing status file is read as silence, which is exactly the truth now).
  # -Filter is a filesystem wildcard, so re-test the name against the exact <worker_id>-s<slot> form
  # before deleting: "abc-s*" would also match a *different* worker whose id happens to start "abc-s".
  $slotRe = '^' + [regex]::Escape($wid) + '-s\d+\.status\.json$'
  foreach ($s in @(Get-ChildItem -LiteralPath $hostsDir -Filter "$wid-s*.status.json" -ErrorAction SilentlyContinue |
                   Where-Object { $_.Name -match $slotRe })) {
    Step "delete stale slot status $($s.Name)" { Remove-Item -LiteralPath $s.FullName -Force -ErrorAction SilentlyContinue }
  }
  Say "retired worker_id=$wid ($Local is left in place; a claim still in running/ is reaped once its status tick stops)"
  exit 0
}

# ---------------------------------------------------------------- PIN.json
$pinPath = Join-Path $setupSrc 'PIN.json'
if (-not (Test-Path -LiteralPath $pinPath)) { throw "PIN.json not found: $pinPath (is $Root reachable?)" }
$pin = Get-Content -LiteralPath $pinPath -Raw -Encoding UTF8 | ConvertFrom-Json
# nastest.ps1 is one of the files deploy_setup.sh puts in ROOT/setup (PROTOCOL section 1.1); it is run
# from LOCAL/setup like every other program (NAS programs are never executed).
$nasTestSrc = Join-Path $setupSrc 'nastest.ps1'
if (-not (Test-Path -LiteralPath $nasTestSrc)) { throw "nastest.ps1 not found: $nasTestSrc (run tools/jobq/deploy_setup.sh first)" }
$juliaVer = [string](Prop $pin 'julia_version' '1.11.9')
$slotFraction = [double](Prop $pin 'slot_fraction' 0.75)
if ($Threads -le 0) { $Threads = [int](Prop $pin 'threads_default' 2) }

# ---------------------------------------------------------------- hardware -> slots
$cpus = @(Get-CimInstance Win32_Processor)
$coresPhys = [int]($cpus | Measure-Object -Property NumberOfCores -Sum).Sum
$coresLog = [int]($cpus | Measure-Object -Property NumberOfLogicalProcessors -Sum).Sum
$cpuName = ($cpus[0].Name -replace '\s+', ' ').Trim()
$ramGb = [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 1)
if ($Slots -le 0) { $Slots = [int][math]::Max(1, [math]::Floor($coresPhys * $slotFraction / $Threads)) }
Say "cpu='$cpuName' cores=$coresPhys/$coresLog ram=${ramGb}GB -> slots=$Slots threads=$Threads (slot_fraction=$slotFraction) julia=$juliaVer"

# ---------------------------------------------------------------- step 1: Git, juliaup, julia channel
Install-IfMissing 'Git.Git' { [bool](Find-GitBash) }
Install-IfMissing 'Julialang.Juliaup' { [bool](Get-Command juliaup.exe -ErrorAction SilentlyContinue) }
Install-JuliaChannel $juliaVer
$bash = Find-GitBash
if (-not $bash) {
  if ($DryRun) { $bash = 'C:\Program Files\Git\bin\bash.exe'; Say "bash.exe not found; dry-run assumes $bash" }
  else { throw 'bash.exe not found after Git install (expected C:\Program Files\Git\bin\bash.exe)' }
}

# ---------------------------------------------------------------- step 2: LOCAL, setup copy, worker.conf
$workerId = Get-ExistingWorkerId
if ($workerId) { Say "keeping WORKER_ID=$workerId from $confPath" } else { $workerId = New-WorkerId; Say "new WORKER_ID=$workerId" }
if ($workerId -notmatch '^[a-z0-9][a-z0-9-]{0,40}$') { throw "invalid worker_id: $workerId" }
# STATUS_INTERVAL was LEASE_INTERVAL before the leases/ directory was deleted (PROTOCOL section 9);
# read the old PIN key as a fallback so an un-updated PIN.json still yields the intended value.
$conf = @(
  "JOBQ_ROOT=$(ConvertTo-MsysPath $Root)",
  "JOBQ_SPOOL=$spoolMsys",
  "JOBQ_LOCAL=$localMsys",
  "WORKER_ID=$workerId",
  "SLOTS=$Slots",
  "THREADS=$Threads",
  "STALL_SECONDS=$(Prop $pin 'stall_seconds' 7200)",
  "MAX_ATTEMPTS=$(Prop $pin 'max_attempts' 5)",
  "HEARTBEAT_INTERVAL=$(Prop $pin 'heartbeat_interval' 180)",
  "RETRY_BACKOFF=$(Prop $pin 'retry_backoff' 30)",
  "DEGRADED_SLEEP=$(Prop $pin 'degraded_sleep' 600)"
) -join "`n"
# code/ (extracted content-addressed trees, PROTOCOL section 1.3) replaced repos/ when the workers
# stopped using git. hosts/ is created here so the very first status write has somewhere to land.
Step "create $Local\{setup,logs,state,work,code} and $hostsDir; copy $setupSrc\* -> $setupDst" {
  foreach ($d in 'setup', 'logs', 'state', 'work', 'code') { New-Item -ItemType Directory -Path (Join-Path $Local $d) -Force | Out-Null }
  New-Item -ItemType Directory -Path $hostsDir -Force | Out-Null
  Copy-Item -Path (Join-Path $setupSrc '*') -Destination $setupDst -Recurse -Force
  if (-not (Test-Path -LiteralPath (Join-Path $setupDst 'nastest.ps1'))) { throw "copy did not produce $setupDst\nastest.ps1" }
}
$confChanged = (Test-Path -LiteralPath $confPath) -and (([IO.File]::ReadAllText($confPath) -replace "`r`n", "`n") -ne ($conf + "`n"))
Step "write $confPath (LF)" { Write-LfFile $confPath ($conf + "`n") }
if ($DryRun) { Write-Host ($conf -replace '(?m)^', '    ') }

# ---------------------------------------------------------------- step 3: NAS test task (runs under the saved credentials)
$nasTestScript = Join-Path $setupDst 'nastest.ps1'      # the distributed copy (never the one on the share)
$nasTestLog = Join-Path $Local 'logs\nastest.log'
$plainPw = $null
if ($DryRun) { Say "[dry-run] would prompt once for the password of $userName (Read-Host -AsSecureString)" }
else {
  $sec = Read-Host -AsSecureString -Prompt "Password for $userName (saved by Task Scheduler so tasks run whether logged on or not)"
  if ($sec.Length -eq 0) { throw 'Empty password: a task that runs whether logged on or not needs one.' }
  $bstr = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($sec)
  try { $plainPw = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($bstr) } finally { [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($bstr) }
}

$psExe = Join-Path $env:SystemRoot 'System32\WindowsPowerShell\v1.0\powershell.exe'
# -Spool as well as -Root: the test has to touch what the workers touch. A share whose root is
# readable but whose spool is not (ACL, read-only, quota) would otherwise pass and fail every job.
$nasArgs = "-NoProfile -ExecutionPolicy Bypass -File `"$nasTestScript`" -Root `"$Root`" -Spool `"$Spool`" -Log `"$nasTestLog`""
$nasTest = [ordered]@{ passed = $false; checked_utc = $null; whoami = $null; userprofile = $null; last_task_result = $null; detail = $null }
Step "register + start task $NasTestTask = `"$psExe`" $nasArgs ; wait <= 120 s for $nasTestLog" {
  if (Test-Path -LiteralPath $nasTestLog) { Remove-Item -LiteralPath $nasTestLog -Force }
  Remove-JobqTasks @($NasTestTask)
  $act = New-ScheduledTaskAction -Execute $psExe -Argument $nasArgs
  $set = New-ScheduledTaskSettingsSet -ExecutionTimeLimit (New-TimeSpan -Minutes 5) -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries -MultipleInstances IgnoreNew
  Register-JobqTask $NasTestTask $act $null $set 'jobq NAS access test (bootstrap.ps1)'
  Start-ScheduledTask -TaskName $NasTestTask
  $deadline = (Get-Date).AddSeconds(120)
  $text = $null
  while ((Get-Date) -lt $deadline) {
    Start-Sleep -Seconds 2
    if (Test-Path -LiteralPath $nasTestLog) {
      $text = [IO.File]::ReadAllText($nasTestLog)
      if ($text -match '(?m)^RESULT ') { break }
    }
  }
  $info = Get-ScheduledTaskInfo -TaskName $NasTestTask
  $nasTest['last_task_result'] = $info.LastTaskResult
  $nasTest['checked_utc'] = UtcNow
  if ($text) {
    if ($text -match '(?m)^whoami=(.*)$') { $nasTest['whoami'] = $Matches[1].Trim() }
    if ($text -match '(?m)^userprofile=(.*)$') { $nasTest['userprofile'] = $Matches[1].Trim() }
    $nasTest['passed'] = [bool]($text -match '(?m)^RESULT PASS')
    $nasTest['detail'] = ($text -replace "`r", '').Trim()
  } else {
    $nasTest['detail'] = "no log within 120 s (task state $((Get-ScheduledTask -TaskName $NasTestTask).State), LastTaskResult $($info.LastTaskResult))"
  }
}

# ---------------------------------------------------------------- host record (section 10 step 6)
$hostRecPath = Join-Path $hostsDir "$workerId.json"
function Build-HostRecord {
  $existing = $null
  if (Test-Path -LiteralPath $hostRecPath) { try { $existing = Get-Content -LiteralPath $hostRecPath -Raw | ConvertFrom-Json } catch { $existing = $null } }
  $now = UtcNow
  return [ordered]@{
    schema = 1; worker_id = $workerId; hostname = $hostName; cpu = $cpuName
    cores_physical = $coresPhys; cores_logical = $coresLog; ram_gb = $ramGb
    slots = $Slots; threads = $Threads; julia_version = $juliaVer
    registered_utc = (Prop $existing 'registered_utc' $now); updated_utc = $now
    bootstrap_user = $userName; root = $Root; spool = $Spool; local = $Local; bash = $bash
    nas_test = $nasTest
    # No 'gates' field: the per-host join gate was abolished on 2026-08-21 (PROTOCOL section 6.5).
    # Any CPU may join; agreement between machines is measured afterwards with tools/agreement_check.py.
    # A 'gates' object left over from a ledger written before that date is deliberately NOT carried
    # forward - nothing reads it, and keeping it would suggest a join gate still exists.
    # 'cpu' stays, but as provenance for the mixed-provenance summary (section 6.5.4), not as a verdict.
  }
}
if (-not $DryRun -and -not $nasTest['passed']) {
  Say "NAS TEST FAILED - workers are NOT registered. Log: $nasTestLog"
  Write-Host $nasTest['detail']
  try { Write-JsonAtomic $hostRecPath (Build-HostRecord); Say "failure recorded in $hostRecPath" }
  catch { Say "could not write $hostRecPath : $($_.Exception.Message)" }
  $nasHost = if ($Root -match '^\\\\([^\\]+)') { $Matches[1] } else { $Root }
  Say "Hint: if the NAS needs other credentials, run as this user: cmdkey /add:$nasHost /user:... /pass:...  then re-run."
  exit 1
}

# ---------------------------------------------------------------- step 4: worker tasks
# Re-run = update in place: Register-ScheduledTask -Force rewrites the definition while a running instance keeps its
# processes, and Start-ScheduledTask is a no-op on a running task (IgnoreNew). Stopping first would orphan worker.sh
# and julia (the scheduler ends only the bin\bash.exe launcher) and the fresh worker would RECOVER the slot's claim
# while the orphan is still writing to the same work dir. Only slots >= Slots are stopped (whole tree) + unregistered.
$workerSh = ConvertTo-MsysPath (Join-Path $Local 'setup\worker.sh')
$stale = @(Get-ScheduledTask -TaskName $WorkerTaskGlob -ErrorAction SilentlyContinue |
  Where-Object { $s = Get-WorkerSlot $_.TaskName; $s -lt 0 -or $s -ge $Slots } | ForEach-Object { $_.TaskName })
if ($stale.Count -gt 0) { Remove-JobqTasks $stale }
if ($confChanged -and @(Get-ScheduledTask -TaskName $WorkerTaskGlob -ErrorAction SilentlyContinue | Where-Object { $_.State -eq 'Running' }).Count -gt 0) {
  Say 'NOTE: worker.conf changed - running workers keep the values they read at start (to restart them now: -Remove, then run again)'
}
for ($k = 0; $k -lt $Slots; $k++) {
  $name = "jobq-worker-s$k"
  $delay = 'PT' + (60 + 60 * $k) + 'S'
  $cmd = (ConvertTo-ShWord $workerSh) + " $k"
  if ($localMsys -ne '/c/jobq') { $cmd = 'JOBQ_LOCAL=' + (ConvertTo-ShWord $localMsys) + ' ' + $cmd }   # worker.sh finds LOCAL/worker.conf via env JOBQ_LOCAL (default /c/jobq)
  $argStr = "-lc `"$cmd`""
  Step "register $name = `"$bash`" $argStr (AtStartup delay $delay; PT0S limit; batteries ok; not idle-only; IgnoreNew; StartWhenAvailable; restart 999 x PT1M; in place if it exists) + Start-ScheduledTask (no-op if running)" {
    $act = New-ScheduledTaskAction -Execute $bash -Argument $argStr
    $trg = New-ScheduledTaskTrigger -AtStartup
    $trg.Delay = $delay
    # -Priority 7 (BelowNormal) is the scheduler default; it is stated here as a decision, not left
    # to chance: an interactive user of the PC must always win against an 8-slot fleet of julia.
    $set = New-ScheduledTaskSettingsSet -ExecutionTimeLimit ([TimeSpan]::Zero) -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
      -RunOnlyIfIdle:$false -MultipleInstances IgnoreNew -StartWhenAvailable -RestartCount 999 -RestartInterval (New-TimeSpan -Minutes 1) `
      -Priority 7
    Register-JobqTask $name $act $trg $set "jobq worker slot $k (bootstrap.ps1)"
    Start-ScheduledTask -TaskName $name
    Say "  $name state=$((Get-ScheduledTask -TaskName $name).State)"
  }
}

$plainPw = $null                                       # the decoded password dies with the last Register-JobqTask above

# ---------------------------------------------------------------- step 5: power
Step 'powercfg /change standby-timeout-ac 0 ; powercfg /hibernate off' {
  Write-Host (Invoke-Native 'powercfg' @('/change', 'standby-timeout-ac', '0'))
  Write-Host (Invoke-Native 'powercfg' @('/hibernate', 'off'))
}

# ---------------------------------------------------------------- step 6: host record
Step "write host record $hostRecPath (tmp + rename)" { Write-JsonAtomic $hostRecPath (Build-HostRecord) }
if ($DryRun) { Write-Host (((Build-HostRecord) | ConvertTo-Json -Depth 8) -replace '(?m)^', '    ') }

Say "done: worker_id=$workerId slots=$Slots threads=$Threads julia=$juliaVer bash=$bash"
if ($DryRun) { Say 'dry-run: nothing was changed' }
exit 0
