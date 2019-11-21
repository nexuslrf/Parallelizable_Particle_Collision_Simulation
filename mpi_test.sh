#!/bin/bash
src=$1
mpicc ${src}.c -o ${src} -lm
is=(1 2 3)
ns=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
# for i in ${is[@]}; do
#     for n in ${ns[@]}; do
#         echo "n=${n}"
#         mpirun -np ${n} ./${src} < testcases_visualizer/testcase${i}.in > inputs_mpi.out
#         tail -n 2 inputs_mpi.out
#         echo ""
#     done
# done
for n in ${ns[@]}; do
    echo "n=${n}"
    mpirun -np ${n} ./${src} < inputs_20.txt > inputs_mpi.out
    tail -n 2 inputs_mpi.out
    echo ""
done