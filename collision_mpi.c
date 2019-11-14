#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <xmmintrin.h>
#include <string.h>
#include <math.h>
#include "timer.h"
#include <mpi.h>

#define RAND01 (rand()%2)
#define eps 1e-10

typedef enum {
    MODE_PRINT,
    MODE_PERF
} simulation_mode_t;

typedef struct
{
    double x;
    double y;
    double vx;
    double vy;   
    int colli_p;
    int colli_w;
}  Particle;

MPI_Datatype Particle_Type;
MPI_Datatype type[3] = { 
    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
    MPI_INT, MPI_INT };
int blocklen[6] = { 1, 1, 1, 1, 1, 1 };
MPI_Aint offsets[6];
offsets[0] = offsetof(Particle, x);
offsets[1] = offsetof(Particle, y);
offsets[2] = offsetof(Particle, vx);
offsets[3] = offsetof(Particle, vy);
offsets[4] = offsetof(Particle, colli_p);
offsets[5] = offsetof(Particle, colli_w);
MPI_Type_create_struct(6, blocklen, offsets, type, &Particle_Type);
MPI_Type_commit(&Particle_Type);

typedef struct
{
    int pa;
    int pb;
    double time;
} Collision;

int n, l, r, s, bnd_far, r_sq_4, num_cmp;
int num_slave, myid, chunk_size;
Particle *particles, *P_a, *P_b;

#define MASTER_ID 0

int compare (const void * a, const void * b)
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

double doubleRand(double min, double max) // return [min, max] double vars
{
    return min+(max-min)*(rand() / (double)RAND_MAX);
}

void bound_pos(Particle *p)
{
    int bnd_far = l-r;
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

/*
Master's job: 
    * init resource; 
    * scatter data;
    * gather res;
*/
void master()
{
    StartTimer();
    // srand((unsigned)time(NULL));
    srand(0);
    int i, j, k, step;
    double x, y, vx, vy;
    simulation_mode_t mode;
    char mode_buf[6];

    scanf("%d", &n);
    scanf("%d", &l);
    scanf("%d", &r);
    scanf("%d", &s);
    scanf("%5s", mode_buf);

    num_cmp = n * (n-1) / 2;
    chunk_size = (n-1) / num_slave + 1;
    
    // Use collective func.
    MPI_Bcast(&n, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    MPI_Bcast(&l, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    MPI_Bcast(&r, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    ///////////
    particles = (Particle *)malloc(n * sizeof(Particle));
    j = 0;
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
    else if(j!=n){
        fprintf(stderr, "Not enough particle parameters!\n");
        exit(1);
    }
    mode = strcmp(mode_buf, "print") == 0 ? MODE_PRINT : MODE_PERF;
    
    /*
        start sim
    */
    for(step = 0; step < s; step++){
        // Scatter particles
        for(v)
    }
}

void slave()
{
    int i, j, k, cnt=0, num_chunk=1;
    for(i=num_slave-1; i>0; i--)
    {
        for(j=0; j<i; j++)
        {
            if(cnt%num_slave == myid)
                num_chunk++;
            cnt++;
        }
    }
    
    MPI_Bcast(&n, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    MPI_Bcast(&l, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    MPI_Bcast(&r, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);

    bnd_far = l - r;
    r_sq_4 = r * r * 4;
    num_cmp = n * (n-1) / 2;
    chunk_size = (n-1) / slave + 1;




}



int main(int argc, char ** argv)
{
    int nprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (nprocs < 2)
    {
        fprintf(stderr, "#Proc >= 2 !\n")
    }
    num_slave = nprocs - 1;
    /*
    Master-Slave Pattern
    */
    if (myid == MASTER_ID)
    {
        master();
    }
    else
    {
        myid--;
        slave();
    }
    MPI_Finalize();
    return 0;
}