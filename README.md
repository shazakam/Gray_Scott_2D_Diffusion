# Gray Scott Diffusion Parallel Fortran implementation

This repository contains the code for the final MA40177 - Scientific Computing Coursework written in Fortran90.
The code is meant to be run in parallel using MPI Routines.
visualise.py contains the code to create an image from the approximated solution at time T.


There are two jobcsript files:
1. jobscript.slm
2. testing_jobscript.slm

jobscript.slm is used to execute the code written for the questions 1-10 in the assignment
testing_jobscript.slm is used to execute the testing.f90 program (on a single processor).

As with the jobscripts there are two main programs:
1. grayscott.f90
2. testing.f90

grayscott.f90 contains the optimised implementation (Q10) using the sparsegather_opt.f90 subroutine
testing.f90 contians a few simple tests for the intitial.f90, create_matrix.f90 and matmult.f90 subroutines.

To compile grayscott.f90 run "make grayscott"
To execute the grayscott executable run "sbatch jobscript.slm"

To compile testing.f90 run "make testing"
To execute the testing executable run "sbatch testing_jobscript.slm"

###################
