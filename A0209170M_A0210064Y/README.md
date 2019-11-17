- Changes made  to Assignment One:

  ​     For OpenMP code, our final code is *collision_openmp.c*. We changed it and now it behaves well under all three testcases. (We tried different approaches to parallelize with OpenMP, but only the best one is uploaded this time.)

  ​     For Cuda code, not sure if there was any problem.  Our code *collision_cuda_balance.cu* looks well under testcases.  And we add more analysis for the cuda part in our report this time.   

- Our OpenMPI code *collision_mpi.c* introduces a lot of variables, which is a bit hard to read. Our basic implementation thought and some important data structure are explained in our report. 

- How to reproduce results?

  ​      For OpenMPI code, to reproduce different execution time varying with process number, run *run_mpi.sh*

  ​       To reproduce different execution time with different particle numbers, run *run_mpi_p.sh*

  ​       To reproduce different execution time with different square size, run *run_mpi_l.sh*

  ​       To reproduce different execution time with different radius, run *run_mpi_r.sh*

  ​       (The above tests are run on machine **soctf-pdc-007**)

  ​      For Cuda and OpenMP codes, reproducing methods are in the report. 

  ​      <u>Please put "timer.h" within the same folder. We use this header to help us output time consumed.</u>

  