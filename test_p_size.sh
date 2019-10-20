#!/bin/bash
run_program=$1
num_block=$2
num_thread=$3
nvcc -arch sm_70 -o collision_cuda ${run_program}.cu
echo "5000 20000 1 20 perf" > input_p.txt
echo "5000 20000 1 20 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "10000 20000 1 20 perf" > input_p.txt
echo "10000 20000 1 20 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "20000 20000 1 20 perf" > input_p.txt
echo "20000 20000 1 20 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "30000 20000 1 20 perf" > input_p.txt
echo "30000 20000 1 20 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "40000 20000 1 20 perf" > input_p.txt
echo "40000 20000 1 20 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "50000 20000 1 20 perf" > input_p.txt
echo "50000 20000 1 20 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""