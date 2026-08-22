<#
jobq NAS access test (tools/jobq/PROTOCOL.md section 10.2 step 3).

  powershell -NoProfile -ExecutionPolicy Bypass -File <LOCAL>\setup\nastest.ps1 -Root <ROOT> [-Spool <SPOOL>] -Log <logfile>

Registered and started by bootstrap.ps1 as the scheduled task "jobq-nastest" so that it runs
under the *worker account's* saved credentials, not under the elevated shell that ran bootstrap:
that is the only way to find out whether the account the workers will run as can actually reach
the share (SMB credentials, drive letters and the Credential Manager are all per-user).

What it proves, in the order a worker needs it:
  whoami / USERPROFILE   which account and profile the scheduler really used
  Test-Path ROOT         the share root is reachable at all
  create / rename / read / delete in SPOOL\hosts
                         the four operations the whole protocol is built on, exercised where the
                         machine-written files live (SPOOL), not at the human-facing root - the
                         root may be readable while the spool is not (ACL, quota, read-only share).

Writes a plain log for a human, and exits 0 (PASS) / 1 (FAIL). bootstrap.ps1 parses this log:
lines "whoami=", "userprofile=", and the marker "RESULT PASS" / "RESULT FAIL" - keep them.
Windows PowerShell 5.1 and pwsh 7. ASCII only, LF line endings.
#>
param(
  [string]$Root,
  [string]$Spool = '',
  [string]$Log
)
Set-StrictMode -Version 2
$ErrorActionPreference = 'Stop'
# Tab completion adds a trailing '\', which turns the closing quote of -Root "...\" into an escaped
# quote when the task's command line is parsed; trim here too so a hand-typed run behaves.
$Root = "$Root".TrimEnd('\', '/')
$Spool = "$Spool".TrimEnd('\', '/')
if (-not $Spool) { $Spool = Join-Path $Root 'spool' }

$L = New-Object System.Collections.Generic.List[string]
$rc = 1
try {
  $L.Add('started_utc=' + (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ'))
  $L.Add('whoami=' + [Security.Principal.WindowsIdentity]::GetCurrent().Name)
  $L.Add('userprofile=' + $env:USERPROFILE)
  $L.Add("root=$Root")
  $L.Add("spool=$Spool")
  $ok = Test-Path -LiteralPath $Root
  $L.Add("test_path_root=$ok")
  if (-not $ok) { throw "root not reachable: $Root" }
  # The spool is what every worker writes to; the root only has to be readable.
  $hosts = Join-Path $Spool 'hosts'
  if (-not (Test-Path -LiteralPath $hosts)) {
    New-Item -ItemType Directory -Path $hosts -Force -ErrorAction SilentlyContinue | Out-Null
  }
  $hostsOk = Test-Path -LiteralPath $hosts -PathType Container
  $L.Add("test_path_spool_hosts=$hostsOk ($hosts)")
  # Say which half is broken: a share whose root answers Test-Path while the spool does not is the
  # interesting failure (ACL, read-only export, a file in the way), and the raw WriteAllText
  # exception below would not name it.
  if (-not $hostsOk) { throw "spool not usable (root is reachable, $hosts is not): check permissions on $Spool" }
  $tag = '.nastest-' + "$env:COMPUTERNAME".ToLower() + '-' + [guid]::NewGuid().ToString('N').Substring(0, 8)
  $f1 = Join-Path $hosts ($tag + '.tmp'); $f2 = Join-Path $hosts ($tag + '.ok')
  $payload = "jobq nastest $tag"
  [IO.File]::WriteAllText($f1, $payload); $L.Add("create=ok ($f1)")
  Move-Item -LiteralPath $f1 -Destination $f2; $L.Add('rename=ok')
  if ([IO.File]::ReadAllText($f2) -ne $payload) { throw 'readback mismatch' }
  $L.Add('read=ok')
  Remove-Item -LiteralPath $f2; $L.Add('delete=ok')
  $rc = 0
  $L.Add('RESULT PASS')
} catch {
  $L.Add('error=' + $_.Exception.Message)
  $L.Add('RESULT FAIL')
}
$L.Add("exit_code=$rc")
try {
  if ($Log) {
    $d = Split-Path -Parent $Log
    if ($d -and -not (Test-Path -LiteralPath $d)) { New-Item -ItemType Directory -Path $d -Force | Out-Null }
    [IO.File]::WriteAllText($Log, (($L -join "`n") + "`n"), (New-Object System.Text.UTF8Encoding($false)))
  }
} catch { }
Write-Host ($L -join "`n")
exit $rc
