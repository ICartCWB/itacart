@ECHO OFF
REM Minimal build script for the itacart documentation on Windows.
REM The build imports the package from ..\src, so no installation is needed.

pushd %~dp0

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=source
set BUILDDIR=_build

if "%1" == "" goto help
if "%1" == "help" goto help
if "%1" == "strict" goto strict
if "%1" == "live" goto live
if "%1" == "coverage" goto coverage
if "%1" == "clean" goto clean

%SPHINXBUILD% -b %1 %SOURCEDIR% %BUILDDIR%\%1
goto end

:help
echo.html      build the HTML documentation into %BUILDDIR%\html
echo.strict    build with warnings as errors, as CI does
echo.live      rebuild and reload on save (needs sphinx-autobuild)
echo.coverage  list public names that carry no docstring
echo.clean     remove the build output and the generated pages
goto end

:strict
%SPHINXBUILD% -b html -W --keep-going %SOURCEDIR% %BUILDDIR%\html
goto end

:live
sphinx-autobuild %SOURCEDIR% %BUILDDIR%\html --watch ..\src
goto end

:coverage
%SPHINXBUILD% -b coverage %SOURCEDIR% %BUILDDIR%\coverage
type %BUILDDIR%\coverage\python.txt
goto end

:clean
rmdir /s /q %BUILDDIR% 2>nul
rmdir /s /q %SOURCEDIR%\api\generated 2>nul
rmdir /s /q %SOURCEDIR%\_generated 2>nul
goto end

:end
popd
