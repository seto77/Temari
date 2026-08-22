<#
Disable Windows PROCESS_POWER_THROTTLING_EXECUTION_SPEED for one process.

Used only for the D317-10 2x2 experiment.  The default fleet configuration does not call this
script.  -SelfTest applies HighQoS to this PowerShell process and immediately returns control to
the system, proving the API and access rights without leaving a persistent setting.
#>
[CmdletBinding(DefaultParameterSetName = 'Process')]
param(
  [Parameter(Mandatory = $true, ParameterSetName = 'Process')]
  [ValidateRange(1, 2147483647)]
  [int]$ProcessId,
  [Parameter(Mandatory = $true, ParameterSetName = 'Tree')]
  [ValidateRange(1, 2147483647)]
  [int]$TreeRootProcessId,
  [Parameter(Mandatory = $true, ParameterSetName = 'SelfTest')]
  [switch]$SelfTest
)
Set-StrictMode -Version 2
$ErrorActionPreference = 'Stop'

if (-not ('Jobq.NativePowerPolicy' -as [type])) {
  Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;

namespace Jobq {
  [StructLayout(LayoutKind.Sequential)]
  public struct ProcessPowerThrottlingState {
    public UInt32 Version;
    public UInt32 ControlMask;
    public UInt32 StateMask;
  }

  public static class NativePowerPolicy {
    public const UInt32 PROCESS_SET_INFORMATION = 0x0200;
    public const UInt32 PROCESS_QUERY_LIMITED_INFORMATION = 0x1000;
    public const UInt32 PROCESS_POWER_THROTTLING_CURRENT_VERSION = 1;
    public const UInt32 PROCESS_POWER_THROTTLING_EXECUTION_SPEED = 0x1;
    public const Int32 ProcessPowerThrottling = 4;

    [DllImport("kernel32.dll", SetLastError = true)]
    public static extern IntPtr OpenProcess(UInt32 desiredAccess, bool inheritHandle, UInt32 processId);

    [DllImport("kernel32.dll", SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool SetProcessInformation(
      IntPtr process,
      Int32 informationClass,
      ref ProcessPowerThrottlingState information,
      UInt32 informationSize);

    [DllImport("kernel32.dll", SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool GetProcessInformation(
      IntPtr process,
      Int32 informationClass,
      ref ProcessPowerThrottlingState information,
      UInt32 informationSize);

    [DllImport("kernel32.dll", SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool CloseHandle(IntPtr handle);
  }
}
'@
}

function Set-ExecutionSpeedPolicy([int]$Id, [bool]$ManagedByJobq) {
  $h = [Jobq.NativePowerPolicy]::OpenProcess(
    ([Jobq.NativePowerPolicy]::PROCESS_SET_INFORMATION -bor
     [Jobq.NativePowerPolicy]::PROCESS_QUERY_LIMITED_INFORMATION), $false, [uint32]$Id)
  if ($h -eq [IntPtr]::Zero) {
    throw [ComponentModel.Win32Exception]::new(
      [Runtime.InteropServices.Marshal]::GetLastWin32Error(), "OpenProcess($Id) failed")
  }
  try {
    $s = New-Object Jobq.ProcessPowerThrottlingState
    $s.Version = [Jobq.NativePowerPolicy]::PROCESS_POWER_THROTTLING_CURRENT_VERSION
    if ($ManagedByJobq) {
      # HighQoS: take control of EXECUTION_SPEED and turn that throttling mechanism off.
      $s.ControlMask = [Jobq.NativePowerPolicy]::PROCESS_POWER_THROTTLING_EXECUTION_SPEED
      $s.StateMask = 0
    } else {
      # Restore system-managed behavior.
      $s.ControlMask = 0
      $s.StateMask = 0
    }
    $size = [uint32][Runtime.InteropServices.Marshal]::SizeOf($s)
    if (-not [Jobq.NativePowerPolicy]::SetProcessInformation(
        $h, [Jobq.NativePowerPolicy]::ProcessPowerThrottling, [ref]$s, $size)) {
      throw [ComponentModel.Win32Exception]::new(
        [Runtime.InteropServices.Marshal]::GetLastWin32Error(), "SetProcessInformation($Id) failed")
    }
    $q = New-Object Jobq.ProcessPowerThrottlingState
    $q.Version = [Jobq.NativePowerPolicy]::PROCESS_POWER_THROTTLING_CURRENT_VERSION
    if (-not [Jobq.NativePowerPolicy]::GetProcessInformation(
        $h, [Jobq.NativePowerPolicy]::ProcessPowerThrottling, [ref]$q, $size)) {
      throw [ComponentModel.Win32Exception]::new(
        [Runtime.InteropServices.Marshal]::GetLastWin32Error(), "GetProcessInformation($Id) failed")
    }
    if ($ManagedByJobq -and
        (($q.ControlMask -band [Jobq.NativePowerPolicy]::PROCESS_POWER_THROTTLING_EXECUTION_SPEED) -eq 0 -or
         ($q.StateMask -band [Jobq.NativePowerPolicy]::PROCESS_POWER_THROTTLING_EXECUTION_SPEED) -ne 0)) {
      throw ("HighQoS readback failed for pid {0}: ControlMask=0x{1:X} StateMask=0x{2:X}" -f
        $Id, $q.ControlMask, $q.StateMask)
    }
    return ("control=0x{0:X} state=0x{1:X}" -f $q.ControlMask, $q.StateMask)
  } finally {
    [void][Jobq.NativePowerPolicy]::CloseHandle($h)
  }
}

if ($SelfTest) {
  $verified = Set-ExecutionSpeedPolicy $PID $true
  [void](Set-ExecutionSpeedPolicy $PID $false)
  Write-Output "disable_ecoqos selftest PASS pid=$PID readback=$verified"
} elseif ($PSCmdlet.ParameterSetName -eq 'Tree') {
  # Git Bash's $! is an MSYS process. /proc/<pid>/winpid is the bash wrapper, whose descendants are
  # julialauncher.exe and then the real julia.exe. Apply HighQoS to both, and do not report success
  # until the real compute process has been seen.
  $deadline = (Get-Date).AddSeconds(10)
  $applied = @{}
  $readback = @{}
  $sawJulia = $false
  do {
    $rows = @(Get-CimInstance Win32_Process | Select-Object ProcessId, ParentProcessId, Name)
    $ids = @{ ([uint32]$TreeRootProcessId) = $true }
    $changed = $true
    while ($changed) {
      $changed = $false
      foreach ($row in $rows) {
        $id = [uint32]$row.ProcessId
        if (-not $ids.ContainsKey($id) -and $ids.ContainsKey([uint32]$row.ParentProcessId)) {
          $ids[$id] = $true; $changed = $true
        }
      }
    }
    foreach ($row in $rows) {
      $id = [uint32]$row.ProcessId
      if (-not $ids.ContainsKey($id) -or $row.Name -notin @('julialauncher.exe', 'julia.exe') -or $applied.ContainsKey($id)) { continue }
      try {
        $readback[$id] = Set-ExecutionSpeedPolicy ([int]$id) $true
        $applied[$id] = $row.Name
        if ($row.Name -ieq 'julia.exe') { $sawJulia = $true }
      } catch {
        # A launcher may exit between the process snapshot and OpenProcess. The real julia must still
        # be observed and changed; otherwise the caller gets a nonzero result.
        if ($row.Name -ieq 'julia.exe') { throw }
      }
    }
    if (-not $sawJulia) { Start-Sleep -Milliseconds 50 }
  } while (-not $sawJulia -and (Get-Date) -lt $deadline)
  if (-not $sawJulia) { throw "real julia.exe was not found below Win32 pid $TreeRootProcessId within 10 s" }
  Write-Output ("HighQoS applied tree_root={0} targets={1}" -f $TreeRootProcessId,
    (($applied.GetEnumerator() | Sort-Object Name | ForEach-Object {
      "$($_.Value):$($_.Key)[$($readback[$_.Key])]" }) -join ','))
} else {
  $verified = Set-ExecutionSpeedPolicy $ProcessId $true
  Write-Output "HighQoS applied pid=$ProcessId readback=$verified"
}
