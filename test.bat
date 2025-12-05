@echo off
setlocal enabledelayedexpansion

echo Building graph...
python src/naive/dbg_indexer.py build -i sequences_file_list.txt -k 31 -o graph.cdbg
if %errorlevel% neq 0 (
    echo [ERROR] Build failed
    pause
    exit /b
)

echo Querying graph...
python src/naive/dbg_indexer.py query -q data/query.fa -i graph.cdbg -o res_query.txt
if %errorlevel% neq 0 (
    echo [ERROR] Query failed
    pause
    exit /b
)

echo Checking results...
fc /W res_query.txt check_res_query.txt >nul

if %errorlevel%==0 (
    echo Test PASSED
) else (
    echo Test FAILED query results are different
)

echo.
echo Done.
pause
