#!/bin/bash
run_program=$2
input_file=$1
num_block=$3
num_thread=$4
nvcc -o run_prog ${run_program}.cu
cp $input_file inputs.txt
./run_prog ${num_block} ${num_thread}
cp outputs.txt outputs_${run_program}_${input_file}