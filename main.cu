#include <stdio.h>
#include <string.h>
#include <math.h>

#define RAND01 (rand()%2)

typedef enum {
    MODE_PRINT,
    MODE_PERF
} simulation_mode_t;

typedef struct {
    int i;
    double x;
    double y;
    double vx;
    double vy;
    int colli_p;
    int colli_w;
    // double x_n;
    // double y_n;
} particle_t;

__constant__ int l, r, s, 
                bnd_far, r_sq_4;
__managed__ particle_t* particles;
__constant__ int n;
int host_n, host_l, host_r, host_s, 
            host_bnd_far, host_r_sq_4;

/* Current implementation: simplest one 
    @TODO optimize it?
    Every thread-> one particle compare N times
*/

__global__ void simulate_step(int num_threads)
{
    int i = blockIdx.x * num_threads + threadIdx.x;
    if(i<n)
        return;
    particle_t *P_a, *P_b, *P_colli;
    P_a = particles + i;
    /* Dummy code that does not check for collision or walls */
    double x_n = P_a->x + P_a->vx;
    double y_n = P_a->y + P_a->vy;
    for(i=0; i<n; i++)
    {

    }
}

__device__ void check_pp_colli(particle_t *P_a, particle_t *P_b, double &time)
{
    double dx1, dy1, Delta, Dx, Dy, dDpdD, dDmdD, DDpDD, lambda;
    time = 1;
    dx1 = P_b->x - P_a->x;
    dy1 = P_b->y - P_a->y;
    // Case 2: overlap at startup:
    ////////////////
    Delta = dx1*dx1 + dy1*dy1;
    if(Delta - r_sq_4<=0 && Delta!=0)
    {
        time = 0;
        cnt++;
    }
    else
    {
        Dx = P_b->vx - P_a->vx;
        Dy = P_b->vy - P_a->vy;
        dDpdD = dx1*Dx + dy1*Dy;
        if(dDpdD>=0) // To judge the right direction
            continue;
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
        // printf("[Debug:lambda]: %f\n", lambda);
        if(lambda<1)
        {
            time = lambda;
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
        particles[i].i = i;
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
        printf("%d %d %d %d %d %d\n", step, i, particles[i].x, particles[i].y,
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

int main(int argc, char** argv)
{
    int i;
    double x, y, vx, vy;
    int num_blocks, num_threads;
    int step;
    simulation_mode_t mode;
    char mode_buf[6];

    if (argc != 3) {
        printf("Usage:\n%s num_blocks num_threads\n", argv[0]);
        return 1;
    }

    num_blocks = atoi(argv[1]);
    num_threads = atoi(argv[2]);

    scanf("%d", &host_n);
    scanf("%d", &host_l);
    scanf("%d", &host_r);
    scanf("%d", &host_s);
    scanf("%5s", mode_buf);
    host_bnd_far = host_l - host_r;
    host_r_sq_4 = host_r * host_r * 4;

    cudaMallocManaged(&particles, sizeof(particle_t) * host_n);

    for (i = 0; i < host_n; i++) {
        particles[i].i = -1;
        particles[i].colli_p = 0;
        particles[i].colli_w = 0;
    }

    while (scanf("%d %lf %lf %lf %lf", &i, &x, &y, &vx, &vy) != EOF) {
        particles[i].i = i;
        particles[i].x = x;
        particles[i].y = y;
        particles[i].vx = vx;
        particles[i].vy = vy;
    }

    if (particles[0].i == -1) {
        randomise_particles(n, r, l);
    }

    mode = strcmp(mode_buf, "print") == 0 ? MODE_PRINT : MODE_PERF;

    /* Copy to GPU constant memory */
    cudaMemcpyToSymbol(n, &host_n, sizeof(n));
    cudaMemcpyToSymbol(l, &host_l, sizeof(l));
    cudaMemcpyToSymbol(r, &host_r, sizeof(r));
    cudaMemcpyToSymbol(s, &host_s, sizeof(s));
    cudaMemcpyToSymbol(bnd_far, &host_bnd_far, sizeof(bnd_far));
    cudaMemcpyToSymbol(r_sq_4, &host_r_sq_4, sizeof(r_sq_4));


    for (step = 0; step < host_s; step++) {
        if (step == 0) {
            print_particles(step);
        }

        /* Call the kernel */
        simulate_step<<<num_blocks, num_threads>>>(num_threads);

        /* Barrier */
        cudaDeviceSynchronize();
    }

    print_statistics(host_s);

    return 0;
}
