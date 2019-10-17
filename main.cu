#include <stdio.h>
#include <string.h>

typedef enum {
    MODE_PRINT,
    MODE_PERF
} simulation_mode_t;

typedef struct {
    int i;
    int x;
    int y;
    int vx;
    int vy;
    int p_collisions;
    int w_collisions;
} particle_t;

__constant__ int l, r, s;

__managed__ particle_t* particles;
__constant__ int n;
int host_n;

__global__ void simulate_step(int num_threads)
{
    int i = blockIdx.x * num_threads + threadIdx.x;
    /* Dummy code that does not check for collision or walls */
    particles[i].x += particles[i].vx;
    particles[i].y += particles[i].vy;
}

__host__ void randomise_particles()
{
    /* TODO Implement randomisation */
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
        printf("%d %d %d %d %d %d %d %d\n", num_step, i, particles[i].x,
            particles[i].y, particles[i].vx, particles[i].vy,
            particles[i].p_collisions, particles[i].w_collisions);
    }
}

int main(int argc, char** argv)
{
    int i, x, y, vx, vy;
    int num_blocks, num_threads;
    int step;
    int host_l, host_r, host_s;
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

    cudaMallocManaged(&particles, sizeof(particle_t) * host_n);

    for (i = 0; i < host_n; i++) {
        particles[i].i = -1;
        particles[i].p_collisions = 0;
        particles[i].w_collisions = 0;
    }

    while (scanf("%d %d %d %d %d", &i, &x, &y, &vx, &vy) != EOF) {
        particles[i].i = i;
        particles[i].x = x;
        particles[i].y = y;
        particles[i].vx = vx;
        particles[i].vy = vy;
    }

    if (particles[0].i == -1) {
        randomise_particles();
    }

    mode = strcmp(mode_buf, "print") == 0 ? MODE_PRINT : MODE_PERF;

    /* Copy to GPU constant memory */
    cudaMemcpyToSymbol(n, &host_n, sizeof(n));
    cudaMemcpyToSymbol(l, &host_l, sizeof(l));
    cudaMemcpyToSymbol(r, &host_r, sizeof(r));
    cudaMemcpyToSymbol(s, &host_s, sizeof(s));

    for (step = 0; step < host_s; step++) {
        if (mode == MODE_PRINT || step == 0) {
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
