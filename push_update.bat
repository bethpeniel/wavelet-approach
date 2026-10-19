@echo off
rem Commits every change in this folder and pushes it to GitHub
rem (https://github.com/bethpeniel/wavelet-approach).
rem Usage: double-click, or run  push_update.bat "what changed"
cd /d "%~dp0"
set MSG=%~1
if "%MSG%"=="" set MSG=Update wavelet app (%date% %time%)
git add -A
git diff --cached --quiet && (echo Nothing to commit - GitHub is already up to date. & goto :end)
git commit -m "%MSG%"
git push
:end
pause
