@echo off
setlocal EnableExtensions
rem ===========================================================================
rem  jobq register.cmd - double-click this to register THIS PC as a jobq worker.
rem  Lives at the share ROOT (\\10.31.108.5\jobq). PROTOCOL.md section 10.1.
rem
rem  ASCII only on purpose: this window uses the console code page (932 here),
rem  so non-ASCII text would come out as mojibake. The Japanese explanation is
rem  in README.txt next to this file.
rem
rem  This file MUST stay CRLF (cmd.exe needs it) - see .gitattributes.
rem  It never relies on the current directory: a cmd window cannot have a UNC
rem  current directory, so every path below is built from %~dp0.
rem
rem  usage: register.cmd [-Slots N] [-Threads T] [-DryRun]
rem         (extra arguments are passed straight through to bootstrap.ps1;
rem          do not use arguments that need quoting)
rem  internal: --orig-user DOMAIN\name   is added by the self-elevation step
rem ===========================================================================

set "JOBQ_SELF=%~f0"
set "JOBQ_HERE=%~dp0"
if "%JOBQ_HERE:~-1%"=="\" set "JOBQ_HERE=%JOBQ_HERE:~0,-1%"
set "JOBQ_BOOTSTRAP=%JOBQ_HERE%\setup\bootstrap.ps1"
set "JOBQ_TITLE=register"

rem ---- arguments -----------------------------------------------------------
rem  The worker account is the account that is logged on NOW, because juliaup
rem  and the Credential Manager are per-user. UAC may elevate as a DIFFERENT
rem  administrator, so the original name is captured here and passed through.
set "JOBQ_USER=%USERDOMAIN%\%USERNAME%"
set "JOBQ_RELAUNCHED=0"
set "JOBQ_EXTRA="
if /i "%~1"=="--orig-user" (
  set "JOBQ_USER=%~2"
  set "JOBQ_RELAUNCHED=1"
  shift
  shift
)
:collect
if "%~1"=="" goto :collected
set JOBQ_EXTRA=%JOBQ_EXTRA% %1
shift
goto :collect
:collected

rem ---- are we elevated? ----------------------------------------------------
rem  net session is the documented probe, but it also fails when the Server
rem  service is stopped; fltmc is a second opinion that needs no service.
set "JOBQ_ELEVATED=0"
net session >nul 2>&1
if not errorlevel 1 set "JOBQ_ELEVATED=1"
if "%JOBQ_ELEVATED%"=="0" (
  fltmc >nul 2>&1
  if not errorlevel 1 set "JOBQ_ELEVATED=1"
)
if "%JOBQ_ELEVATED%"=="1" goto :run
if "%JOBQ_RELAUNCHED%"=="1" (
  echo [WARN] still not elevated after a UAC round-trip - continuing anyway,
  echo        bootstrap.ps1 will say what is missing.
  goto :run
)

rem ---- re-launch this file elevated ----------------------------------------
rem  cmd.exe /c ""<self>" --orig-user "<user>" <extra>"
rem  The outer pair of quotes is what cmd /c strips; the \" are for the argv
rem  parser between cmd and powershell (a bare " would end the -Command string
rem  and the quoting inside would be lost).
echo Requesting administrator rights - please confirm the UAC prompt.
powershell -NoProfile -Command "Start-Process -FilePath cmd.exe -Verb RunAs -ArgumentList '/c \"\"%JOBQ_SELF%\" --orig-user \"%JOBQ_USER%\"%JOBQ_EXTRA%\"'"
if errorlevel 1 (
  echo.
  echo [FAIL] could not start the elevated window ^(UAC cancelled, or no
  echo        administrator account available^). Nothing was changed.
  echo.
  pause
  exit /b 1
)
exit /b 0

rem ---- do the work ---------------------------------------------------------
:run
if not exist "%JOBQ_BOOTSTRAP%" (
  echo.
  echo [FAIL] not found: %JOBQ_BOOTSTRAP%
  echo        Is the share reachable from THIS account? Open %JOBQ_HERE%
  echo        in Explorer and try again.
  echo.
  pause
  exit /b 1
)
echo ==== jobq %JOBQ_TITLE% ====
echo   share root  : %JOBQ_HERE%
echo   worker user : %JOBQ_USER%
echo   bootstrap   : %JOBQ_BOOTSTRAP%
echo   extra args  :%JOBQ_EXTRA%
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%JOBQ_BOOTSTRAP%" -Root "%JOBQ_HERE%" -User "%JOBQ_USER%"%JOBQ_EXTRA%
set "JOBQ_RC=%ERRORLEVEL%"
echo.
if "%JOBQ_RC%"=="0" (
  echo [PASS] jobq %JOBQ_TITLE% finished.
) else (
  echo [FAIL] jobq %JOBQ_TITLE% failed with exit code %JOBQ_RC%.
  echo        Read the lines above. Logs: C:\jobq\logs\
)
echo.
pause
exit /b %JOBQ_RC%
