#!/bin/bash
src="collision_mpi"
n=16
mpicc ${src}.c -o ${src} -lm
echo "10000 5000 1 50 perf" > inputs_p.txt
echo "10000 5000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 10000 1 50 perf" > inputs_p.txt
echo "10000 10000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 100000 1 50 perf" > inputs_p.txt
echo "10000 100000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 1000000 1 50 perf" > inputs_p.txt
echo "10000 1000000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 10000000 1 50 perf" > inputs_p.txt
echo "10000 10000000 1 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
