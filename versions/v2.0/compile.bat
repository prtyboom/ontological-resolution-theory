@echo off
echo ========================================
echo  Compiling ORT v2.0
echo ========================================
echo.

echo [1/3] First LaTeX pass...
pdflatex -interaction=nonstopmode ORT_v2.0.tex

echo.
echo [2/3] Second LaTeX pass (for references)...
pdflatex -interaction=nonstopmode ORT_v2.0.tex

echo.
echo [3/3] Cleaning auxiliary files...
del *.aux 2>nul
del *.log 2>nul
del *.out 2>nul
del *.toc 2>nul
del *.bbl 2>nul
del *.blg 2>nul
del *.synctex.gz 2>nul
del *.fdb_latexmk 2>nul
del *.fls 2>nul

echo.
echo ========================================
echo  Done! Output: ORT_v2.0.pdf
echo ========================================
echo.

if exist ORT_v2.0.pdf (
    echo SUCCESS: PDF generated
    start ORT_v2.0.pdf
) else (
    echo ERROR: PDF not created - check for LaTeX errors
)

pause