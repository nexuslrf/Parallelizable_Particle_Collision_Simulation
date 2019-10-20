#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "timer.h"

#define RAND01 (rand()%2)
#define eps 1e-10

typedef enum {
    MODE_PRINT,
    MODE_PERF
} simulation_mode_t;

typedef struct {
    double x;
    double y;
    double vx;
    double vy;
    int colli_p;
    int colli_w;
    double x_n;
    double y_n;
} particle_t;

typedef struct
{
    int pa;
    int pb;
    double time;
} Collision;

__constant__ int n, l, r, s, bnd_far, r_sq_4;
__managed__ particle_t* particles;
__managed__ int *colli_mat, *colli_queue;
__managed__ Collision *colli_time;
__managed__ int count, real_colli;
int host_n, host_l, host_r, host_s, host_bnd_far, host_r_sq_4;
Collision *colli;


/* Current implementation: simplest one 
    @TODO optimize it?
    Every thread-> one particle compare N times
*/

__global__ void find_collisions(int step, int chunk_size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // printf("[Debug:cuda_thread] step %d thread %d\n",step, i);
    idx = idx * chunk_size;
    int i;
    particle_t *P_a, *P_b;
    for(i=idx;i<idx+chunk_size&&i<n;i++)
    {
        P_a=particles+i;
        P_a->x_n = P_a->x + P_a->vx;
        P_a->y_n = P_a->y + P_a->vy;
        double lambda, lambda_1, lambda_2, wall_colli, dx1, dy1, Delta, Dx, Dy, dDpdD, DDpDD, dDmdD;
        int cnt;
        // Case 1
        lambda_1 = lambda_2 = 2;
        wall_colli = 0;
        if(P_a->x_n<r)
        {
            lambda_1 = (r - P_a->x) / P_a->vx;
            wall_colli = 1;
        }
        else if(P_a->x_n>bnd_far)
        {
            lambda_1 = (bnd_far - P_a->x) / P_a->vx;
            wall_colli = 1;
        }

        if(P_a->y_n<r)
        {
            lambda_2 = (r - P_a->y) / P_a->vy;
            wall_colli = 1;
        }
        else if(P_a->y_n>bnd_far)
        {
            lambda_2 = (bnd_far - P_a->y) / P_a->vy;
            wall_colli = 1;
        }
        ////
        // if(i==16)
        // {
        //     printf("[Debug:cuda] %lf %lf %d\n", lambda_1, lambda_2, wall_colli);
        // }
        ////
        if(wall_colli)
        {
            cnt=atomicAdd(&count, 1); // nice!
            colli_time[cnt].pb = i;
            lambda = lambda_1-lambda_2;
            if(lambda==0) // Cornor collision!
            {
                colli_time[cnt].pa = -1; // -1 to present this case.
                colli_time[cnt].time = lambda_1;
            }
            else if(lambda<0) // x wall collision!
            {
                colli_time[cnt].pa = -2; // -2 to present this case.
                colli_time[cnt].time = lambda_1;
            }
            else if(lambda>0) // y wall collision!
            {
                colli_time[cnt].pa = -3; // -3 to present this case.
                colli_time[cnt].time = lambda_2;
            }
        }
        ///////////////
        for(int j=i+1; j<n; j++)
        {
            P_b = particles+j;
            dx1 = P_b->x - P_a->x;
            dy1 = P_b->y - P_a->y;
            // early stop
            Dx = P_b->vx - P_a->vx;
            Dy = P_b->vy - P_a->vy;
            dDpdD = dx1*Dx + dy1*Dy;
            if(dDpdD>=0) // To judge the right direction
                continue;
            // Case 2: overlap at startup:
            ////////////////
            Delta = dx1*dx1 + dy1*dy1;
            if(Delta - r_sq_4<=0 && Delta!=0)
            {
                cnt=atomicAdd(&count, 1);
                colli_time[cnt].time = 0;
                colli_time[cnt].pa = i;
                colli_time[cnt].pb = j; // pa always smaller than pb
                continue; // no need to further detect.
            }
            ////////////////
            // Case 3: Normal collision case
            ////////////////
            DDpDD = Dx*Dx + Dy*Dy;
            dDmdD = dx1*Dy - dy1*Dx;
            Delta = r_sq_4*DDpDD - dDmdD*dDmdD;
            if(Delta<=0)
                continue;
            Delta = sqrt(Delta);
            lambda = (-dDpdD - Delta)/DDpDD;
            if(lambda<1)
            {
                cnt=atomicAdd(&count, 1);
                colli_time[cnt].time = lambda;
                colli_time[cnt].pa = i;
                colli_time[cnt].pb = j;
            }
            ////////////////
        }
    }
}

__host__ void find_real_collisions()
{
    for(int i=0;i<count;i++)
    {
        colli = colli_time+i;
        ////
        // if(1 && (colli->pa == 2||colli->pb==2))
        // {
        //     printf("[Debug:inconsist] %d %d %10.8f\n",colli->pa, colli->pb, colli->time);
        // }
        ////
        if(colli->pa<0){ //wall collision
            if(!colli_mat[colli->pb])
            {
                colli_mat[colli->pb] = 1;
                colli_queue[real_colli++] = i;
                particles[colli->pb].colli_w++;
            }
        }
        else if(!colli_mat[colli->pa])
        {
            if(!colli_mat[colli->pb])
            {
                colli_mat[colli->pa] = 1;
                colli_mat[colli->pb] = 1;
                colli_queue[real_colli++] = i;
                particles[colli->pa].colli_p++;
                particles[colli->pb].colli_p++;
            }
        }
    }
}

__device__ void bound_pos(particle_t *p)
{
    double tx=0,ty=0;
    if(p->x_n>bnd_far)
        tx = (p->x_n-bnd_far)/p->vx;
    else if(p->x_n<r)
        tx = (p->x_n-r)/p->vx;
    if(p->y_n>bnd_far)
        ty = (p->y_n-bnd_far)/p->vy;
    else if(p->y_n<r)
        ty = (p->y_n-r)/p->vy;

    tx =ty = tx>ty?tx:ty;
    p->x_n = p->x_n - tx*p->vx;
    p->y_n = p->y_n - ty*p->vy;
}

__global__ void update_particle(int chunk_size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    idx = idx * chunk_size;
    int i;
    particle_t *P_a, *P_b;
    Collision *Colli;
    double Dx ,Dy, Delta, dx1, dy1, dx2, dy2, DDpDD;
    for(i=idx;i<idx+chunk_size&&i<real_colli;i++)
    {
        Colli = colli_time + colli_queue[i];
        if(Colli->pa==-1) // Cornor colli;
        {
            P_a = particles + Colli->pb;
            P_a->vx = -1*P_a->vx;
            P_a->vy = -1*P_a->vy;
            P_a->x_n = P_a->x+(1-2*Colli->time)*P_a->vx;
            P_a->y_n = P_a->y+(1-2*Colli->time)*P_a->vy;
            bound_pos(P_a);
        }
        else if(Colli->pa==-2)//  X wall colli;
        {
            P_a = particles + Colli->pb;
            P_a->vx = -1*P_a->vx;
            P_a->x_n = P_a->x+(1-2*Colli->time)*P_a->vx;
            bound_pos(P_a);
        }
        else if(Colli->pa==-3)// Y wall colli;
        {
            P_a = particles + Colli->pb;
            P_a->vy = -1*P_a->vy;
            P_a->y_n = P_a->y+(1-2*Colli->time)*P_a->vy;
            bound_pos(P_a);
        }
        else // P-P colli;
        {
            P_a = particles + Colli->pa;
            P_b = particles + Colli->pb;
            P_a->x_n = P_a->x + Colli->time*P_a->vx;
            P_a->y_n = P_a->y + Colli->time*P_a->vy;
            P_b->x_n = P_b->x + Colli->time*P_b->vx;
            P_b->y_n = P_b->y + Colli->time*P_b->vy;
            Dx = P_b->x_n - P_a->x_n;
            Dy = P_b->y_n - P_a->y_n;
            Delta = 1 - Colli->time;
            /* To reduce var:
            dx1: nv1; dy1: tv1;
            dx2: nv2; dy2: tv2;
            */
            dx1 = Dx*P_a->vx + Dy*P_a->vy;
            dy1 = Dx*P_a->vy - Dy*P_a->vx;
            dx2 = Dx*P_b->vx + Dy*P_b->vy;
            dy2 = Dx*P_b->vy - Dy*P_b->vx;
            DDpDD = Dx*Dx + Dy*Dy;
            if(DDpDD!=0)
            {
                // Update velocities
                P_a->vx = (dx2*Dx-dy1*Dy)/DDpDD;
                P_a->vy = (dx2*Dy+dy1*Dx)/DDpDD;
                P_b->vx = (dx1*Dx-dy2*Dy)/DDpDD;
                P_b->vy = (dx1*Dy+dy2*Dx)/DDpDD;
            }
            // Update position
            P_a->x_n = P_a->x_n + Delta*P_a->vx;
            P_a->y_n = P_a->y_n + Delta*P_a->vy;
            bound_pos(P_a);
            P_b->x_n = P_b->x_n + Delta*P_b->vx;
            P_b->y_n = P_b->y_n + Delta*P_b->vy;
            bound_pos(P_b);
        }
    }
}

__host__ double doubleRand(double min, double max) // return [min, max] double vars
{
    return min+(max-min)*(rand() / (double)RAND_MAX);
}

__host__ void randomise_particles()
{
    /* Implement randomisation */
    for(int i=0; i<host_n; i++)
    {
        particles[i].x = doubleRand(host_r,host_bnd_far);
        particles[i].y = doubleRand(host_r,host_bnd_far);
        particles[i].vx = (1 - 2*RAND01)*doubleRand(host_l/(double)8.0/host_r,host_l/(double)4.0);
        particles[i].vy = (1 - 2*RAND01)*doubleRand(host_l/(double)8.0/host_r,host_l/(double)4.0);
    }
}

__host__ void print_particles(int step)
{
    int i;
    for (i = 0; i < host_n; i++) {
        printf("%d %d %10.8lf %10.8lf %10.8lf %10.8lf\n", step, i, particles[i].x, particles[i].y,
            particles[i].vx, particles[i].vy);
    }
}

__host__ void print_statistics(int num_step)
{
    int i;
    for (i = 0; i < host_n; i++) {
        printf("%d %d %10.8lf %10.8lf %10.8lf %10.8lf %d %d\n", num_step, i, particles[i].x,
            particles[i].y, particles[i].vx, particles[i].vy,
            particles[i].colli_p, particles[i].colli_w);
    }
}

__host__ int compare (const void * a, const void * b)
{
    Collision *colli_A = (Collision*)a;
    Collision *colli_B = (Collision*)b;
    double cmpf = colli_A->time - colli_B->time;
    if(fabs(cmpf)<eps)
    {
        int cmpt = colli_A->pa - colli_B->pa;
        if(cmpt!=0)
            return cmpt;
        else
            return colli_A->pb - colli_B->pb;
    }
    else
        return cmpf<0?-1:1;
}

__host__ void check_cuda_errors()
{
    cudaError_t rc;
    rc = cudaGetLastError();
    if (rc != cudaSuccess)
    {
        printf("Last CUDA error %s\n", cudaGetErrorString(rc));
    }
}

int main(int argc, char** argv)
{
    StartTimer();
    int i,j=0;
    double x, y, vx, vy;
    int num_blocks, num_threads, chunk_p, chunk_c, total_threads;
    int step;
    simulation_mode_t mode;
    char mode_buf[6];

    freopen("./inputs.txt","r",stdin);
    freopen("./outputs.txt","w",stdout);
    srand(0);
    if (argc != 3) {
        printf("Usage:\n%s num_blocks num_threads\n", argv[0]);
        return 1;
    }

    num_blocks = atoi(argv[1]);
    num_threads = atoi(argv[2]);
    total_threads = num_threads * num_blocks;

    scanf("%d", &host_n);
    scanf("%d", &host_l);
    scanf("%d", &host_r);
    scanf("%d", &host_s);
    scanf("%5s", mode_buf);
    host_bnd_far = host_l - host_r;
    host_r_sq_4 = host_r * host_r * 4;
    chunk_p = (host_n-1) / total_threads + 1;

    cudaMallocManaged((void**)&particles, sizeof(particle_t) * host_n);
    cudaMallocManaged((void**)&colli_mat,sizeof(int) * host_n);
    cudaMallocManaged((void**)&colli_queue,sizeof(int) * host_n);
    cudaMallocManaged((void**)&colli_time, sizeof(Collision) * host_n*(host_n+1)/2);

    while (scanf("%d %lf %lf %lf %lf", &i, &x, &y, &vx, &vy) != EOF) {
        j++;
        particles[i].x = x;
        particles[i].y = y;
        particles[i].vx = vx;
        particles[i].vy = vy;
        particles[i].colli_p = 0;
        particles[i].colli_w = 0;
    }
    if (j==0) {
        randomise_particles();
    }
    else if(j!=host_n){
        fprintf(stderr, "Not enough particle parameters!\n");
        exit(1);
    }

    mode = strcmp(mode_buf, "print") == 0 ? MODE_PRINT : MODE_PERF;

    /* Copy to GPU constant memory */
    cudaMemcpyToSymbol(n, &host_n, sizeof(n));
    cudaMemcpyToSymbol(l, &host_l, sizeof(l));
    cudaMemcpyToSymbol(r, &host_r, sizeof(r));
    cudaMemcpyToSymbol(s, &host_s, sizeof(s));
    cudaMemcpyToSymbol(bnd_far, &host_bnd_far, sizeof(bnd_far));
    cudaMemcpyToSymbol(r_sq_4, &host_r_sq_4, sizeof(r_sq_4));
    check_cuda_errors();

    for (step = 0; step < host_s; step++) {
        if (step == 0) {
            print_particles(step);
        }
        count=0; //initialize collision numbers every step
        real_colli=0;
        cudaMemset(colli_mat, 0, sizeof(int) * host_n);
        /* Call the kernel */
        find_collisions<<<num_blocks, num_threads>>>(step, chunk_p); 
        /* Barrier */
        cudaDeviceSynchronize();
        // find real collisions
        qsort(colli_time, count, sizeof(Collision), compare);
        find_real_collisions();
        /* Call the kernel */
        chunk_c =(real_colli-1)/total_threads+1; 
        update_particle<<<num_blocks, num_threads>>>(chunk_c); 
        /* Barrier */
        cudaDeviceSynchronize();
        for(i=0;i<host_n;i++)
        {
            particles[i].x = particles[i].x_n;
            particles[i].y = particles[i].y_n;
        }
        if(mode==MODE_PRINT)
            print_particles(step+1);
    }

    print_statistics(host_s);
    double exec_time=GetTimer();
    printf("Time elapsed:%lf",exec_time);

    fclose(stdin);
    fclose(stdout);
    cudaFree(particles);
    cudaFree(colli_time);
    cudaFree(colli_mat);
    cudaFree(colli_queue);

    return 0;
}
