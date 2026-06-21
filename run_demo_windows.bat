@echo off
title Smith-Waterman Local Alignment Demo
echo ==========================================================
echo           SMITH-WATERMAN LOCAL SEQUENCE ALIGNMENT
echo ==========================================================
echo.
echo [STEP 1] Running Sequential Baseline...
echo ----------------------------------------------------------
seq_align.exe 5000 5000 42
echo.
echo [STEP 2] Running Parallel MPI Version (3 Processes)...
echo ----------------------------------------------------------
"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" -n 3 mpi_align.exe 5000 5000 42 128
echo.
echo ==========================================================
echo Demo completed. Press any key to exit.
echo ==========================================================
pause > nul
