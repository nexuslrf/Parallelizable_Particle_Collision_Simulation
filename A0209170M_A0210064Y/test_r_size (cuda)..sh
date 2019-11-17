#!/bin/bash
run_program=$1
num_block=$2
num_thread=$3
nvcc -arch sm_70 -o collision_cuda ${run_program}.cu
echo "10000 20000 1 50 perf" > input_p.txt
echo "10000 20000 1 50 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "10000 20000 4 50 perf" > input_p.txt
echo "10000 20000 4 50 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "10000 20000 64 50 perf" > input_p.txt
echo "10000 20000 64 50 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "10000 20000 128 50 perf" > input_p.txt
echo "10000 20000 128 50 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "10000 20000 256 50 perf" > input_p.txt
echo "10000 20000 256 50 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
echo "10000 20000 512 50 perf" > input_p.txt
echo "10000 20000 512 50 perf"
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""
./collision_cuda ${num_block} ${num_thread} < input_p.txt > outputs_${run_program}_${input_file}
tail -n 2 outputs_${run_program}_${input_file}
echo ""