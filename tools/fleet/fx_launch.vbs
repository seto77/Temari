' Launch launch_fx.sh (dataset-factors v1 generation + endpoint stage 2) with no console window.
' Same mechanism as ep_launch.vbs (validated 2026-08-12..14). Keep this file ASCII-only.
Set sh = CreateObject("WScript.Shell")
sh.Run """C:\Program Files\Git\bin\bash.exe"" -lc ""bash /c/tmp/temari_factors_2026-08-16/launch_fx.sh""", 0, False
