#!/bin/bash
run_program=$1
input_file=$2
gcc -fopenmp -o run_prog ${run_program}.c -lm
cp $input_file inputs.txt
perf stat -- ./run_prog
cp outputs.txt outputs_${run_program}_${input_file}