#!/bin/bash
run_program=$2
input_file=$1
nvcc -o run_prog ${run_program}.c -lm
cp $input_file inputs.txt
./run_prog
cp outputs.txt outputs_${run_program}_${input_file}