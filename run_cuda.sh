#!/bin/bash
run_program=$2
input_file=$1
num_block=$3
num_thread=$4
nvcc -arch sm_70 -o collision_cuda ${run_program}.cu
#cp $input_file inputs.txt
./collision_cuda ${num_block} ${num_thread}
cp outputs.txt outputs_${run_program}_${input_file}
