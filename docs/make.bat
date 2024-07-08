@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

set SPHINXBUILD=sphinx-build
set SOURCEDIR=source
set BUILDDIR=build
set PACKAGE_DIR=..\src

if not "%SPHINXBUILD%"=="" (
    %SPHINXBUILD% >NUL 2>NUL
)
if errorlevel 9009 (
    echo.
    echo The 'sphinx-build' command was not found. Make sure you have Sphinx
    echo installed, then set the SPHINXBUILD environment variable to point
    echo to the full path of the 'sphinx-build' executable. Alternatively, you
    echo may add the Sphinx directory to PATH.
    echo.
    echo If you don't have Sphinx installed, grab it from
    echo https://www.sphinx-doc.org/
    exit /b 1
)

REM Custom target for generating API docs
if "%1" == "apidoc" (
    sphinx-apidoc --templatedir="source/_templates/public" -f -o "%SOURCEDIR%/public_api" "%PACKAGE_DIR%"
    sphinx-apidoc --private --templatedir="source/_templates/private" -f -o "%SOURCEDIR%/private_api" "%PACKAGE_DIR%"
    goto end
)

REM Python script to customize .rst files
if "%1" == "customize-rst" (
    python rename.py
    goto end
)

if "%1" == "" goto help

REM Modify the html target to depend on the new apidoc and customize-rst targets
if "%1" == "html" (
    call :apidoc
    call :customize-rst
)

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end

:apidoc
sphinx-apidoc -o "%SOURCEDIR%/public_api" "%PACKAGE_DIR%"
sphinx-apidoc --private -o "%SOURCEDIR%/private_api" "%PACKAGE_DIR%"
goto :eof

:customize-rst
python prepend_api_designation.py
goto :eof

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%

:end
popd