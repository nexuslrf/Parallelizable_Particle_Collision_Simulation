#!/bin/bash
src="collision_mpi"
n=16
mpicc ${src}.c -o ${src} -lm
echo "5000 20000 1 50 perf" > inputs_p.txt
echo "5000 20000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "20000 20000 1 50 perf" > inputs_p.txt
echo "20000 20000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "30000 20000 1 50 perf" > inputs_p.txt
echo "30000 20000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "40000 20000 1 50 perf" > inputs_p.txt
echo "40000 20000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "50000 20000 1 50 perf" > inputs_p.txt
echo "50000 20000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""

