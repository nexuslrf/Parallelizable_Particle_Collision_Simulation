#!/bin/bash
src="collision_mpi"
mpicc ${src}.c -o ${src} -lm
n=2
echo "n=2"
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
n=3
echo "n=3"
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
n=4
echo "n=4"
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
n=8
echo "n=8"
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
n=16
echo "n=16"
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
n=32
echo "n=32"
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
n=64
echo "n=64"
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
mpirun -np ${n} ./${src} < inputs.txt > inputs_mpi.out
tail -n 2 inputs_mpi.out
echo ""
