#!/bin/bash
run_program=$2
input_file=$1
num_block=$3
num_thread=$4
nvcc -arch sm_70 -o collision_cuda ${run_program}.cu
./collision_cuda ${num_block} ${num_thread} < $input_file > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < $input_file > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < $input_file > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
# nvprof  ./collision_cuda ${num_block} ${num_thread} < $input_file > outputs_${run_program}_${input_file}
# tail -n 2 outputs_${run_program}_${input_file}
# echo ""