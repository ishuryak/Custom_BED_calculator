@echo off
rem Starts the calculator in the default browser.
rem Uses RSCRIPT if set, else Rscript on the PATH, else the newest R under Program Files.
cd /d "%~dp0"
set "RS=%RSCRIPT%"
if not defined RS for /f "delims=" %%i in ('where Rscript 2^>nul') do if not defined RS set "RS=%%i"
if not defined RS for /f "delims=" %%d in ('dir /b /ad /o-n "%ProgramFiles%\R\R-*" 2^>nul') do if not defined RS if exist "%ProgramFiles%\R\%%d\bin\Rscript.exe" set "RS=%ProgramFiles%\R\%%d\bin\Rscript.exe"
if not defined RS (
  echo Rscript was not found. Install R, or set RSCRIPT to the full path of Rscript.exe.
  pause
  exit /b 1
)
echo Using "%RS%"
"%RS%" -e "shiny::runApp('.', launch.browser = TRUE)"
pause
