#!/bin/bash
src="collision_mpi"
n=16
mpicc ${src}.c -o ${src} -lm
echo "10000 20000 4 50 perf" > inputs_p.txt
echo "10000 20000 4 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 20000 64 50 perf" > inputs_p.txt
echo "10000 20000 64 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 20000 128 50 perf" > inputs_p.txt
echo "10000 20000 128 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 20000 256 50 perf" > inputs_p.txt
echo "10000 20000 256 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 20000 512 50 perf" > inputs_p.txt
echo "10000 20000 512 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
echo "10000 20000 1024 50 perf" > inputs_p.txt
echo "10000 20000 1024 50 perf"
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs_p.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""

