# Knight Exchange Problem Solver

## Project Overview
This repository contains implementations of parallel algorithms for solving the Knight Exchange Problem on a chessboard.

## Problem Definition
- Exchange positions of white and black knights on a chessboard
- Input: Chessboard size (m,n), number of knights (k), initial positions (W, B)
- Goal: Find the shortest sequence of moves to exchange knight positions

## Algorithms Implemented
1. Sequential Algorithm (Branch and Bound Depth-First Search)
2. Task Parallelism using OpenMP
3. Data Parallelism using OpenMP
4. Task Parallelism using MPI

## Performance Evaluation
- Metrics: Speed-up, Computational Cost, Parallel Efficiency
- Test instances: 7x3, 4x5, and 8x3 chessboards

## Key Findings
- OpenMP: Optimal performance with 8 threads
- MPI: Best balance achieved with 24 threads (8 per slave node)
- Diminishing returns observed with higher thread counts
