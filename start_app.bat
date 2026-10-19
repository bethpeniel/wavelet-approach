@echo off
rem Starts the Wavelet Approach Shiny server on port 4321 (if not already running)
rem and opens open_app.html in the default browser.
set RSCRIPT="C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
set APPDIR=%~dp0

netstat -ano | findstr /R /C:":4321 .*LISTENING" >nul
if %errorlevel%==0 (
  echo Server already running on port 4321.
) else (
  echo Starting the wavelet app server on port 4321 ...
  start "Wavelet app server" /min %RSCRIPT% -e "shiny::runApp('%APPDIR:\=/%', port = 4321, launch.browser = FALSE)"
  timeout /t 8 /nobreak >nul
)
start "" "%APPDIR%open_app.html"
