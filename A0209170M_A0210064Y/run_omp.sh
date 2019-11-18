#!/bin/bash
run_program=$2
input_file=$1
gcc -fopenmp -o run_prog ${run_program}.c -lm
./run_prog < $input_file > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./run_prog < $input_file > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./run_prog < $input_file > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
# perf stat -- ./collision_cuda < $input_file > outputs_${run_program}_${input_file}
# tail -n 2 outputs_${run_program}_${input_file}
# echo ""