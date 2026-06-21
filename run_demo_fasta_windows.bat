@echo off
title Smith-Waterman FASTA Local Alignment Demo
echo ==========================================================
echo           SMITH-WATERMAN FASTA LOCAL ALIGNMENT DEMO
echo ==========================================================
echo.
echo Dataset: SARS-CoV-2 (COVID-19) vs MERS-CoV Genomes
echo.
echo [STEP 1] Running Sequential Baseline...
echo ----------------------------------------------------------
seq_align.exe data\sars_cov_2.fasta data\mers_cov.fasta
echo.
echo [STEP 2] Running Parallel MPI Version (3 Processes)...
echo ----------------------------------------------------------
"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" -n 3 mpi_align.exe data\sars_cov_2.fasta data\mers_cov.fasta 128
echo.
echo ==========================================================
echo Demo completed. Press any key to exit.
echo ==========================================================
pause > nul
