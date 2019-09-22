#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <xmmintrin.h>
#include <string.h>

#define RAND01 (rand()%2)

typedef struct
{
    double x;
    double y;
    double vx;
    double vy;   
    // Speculative pts
    /////////////////
    double x_n;
    double y_n;
    double vx_n;
    double vy_n;
    /////////////////
}  Particle;
int N, L, r, S;
char mode[6];

double doubleRand(double min, double max) // return [min, max] double vars
{
    return min+(max-min)*(rand() / (double)RAND_MAX);
}
double bound_pos(double p)
{
    int bnd_far = L-r;
    if(p>bnd_far)
        return bnd_far;
    else if(p<r)
        return r;
    else
        return p;
}

int main()
{
    srand(0);
    Particle *particles;
    freopen("./inputs.txt","r",stdin);
    scanf("%d %d %d %d %s",&N, &L, &r, &S, mode);
    particles = (Particle *)malloc(N * sizeof(Particle));
    int i =0,j, t, idx;
    double x, y, vx, vy;
    while(scanf("%d %lf %lf %lf %lf", &idx,&x,&y,&vx,&vy)!=EOF)
    {
        i++;
        particles[idx].x = x;
        particles[idx].y = y;
        particles[idx].vx = vx;
        particles[idx].vy = vy;
    }
    if(i==0)
    {
        for(;i<N;i++)
        {
            particles[i].x = doubleRand(r,L-r);
            particles[i].y = doubleRand(r,L-r);
            particles[i].vx = (1 - 2*RAND01)*doubleRand(L/(double)8.0/r,L/(double)4.0);
            particles[i].vy = (1 - 2*RAND01)*doubleRand(L/(double)8.0/r,L/(double)4.0);
        }
    }
    else if(i!=N)
    {
        fprintf(stderr, "Not enough particle parameters!\n");
        exit(1);
    }
    // Print initial position
    for(i=0; i<N; i++)
        printf("0 %d %10.8lf %10.8lf %10.8lf %10.8lf\n", i, 
                    particles[i].x, particles[i].y, particles[i].vx, particles[i].vy);
    // Naive Method
    /* Basic Idea:
    1. Build collision table (N+1)x(N+1), N is the wall. To check whether the collision happens.
    2. Build collision time table for processing.
    🎯
    */
    int **colli_mat;
    colli_mat = (int**)malloc((N+1)*sizeof(int*));
    for(i=0; i<N+1; i++)
    {    
        colli_mat[i]=(int*)malloc((N+1)*sizeof(int));
        memset(colli_mat[i],0,(N+1)*sizeof(int));
    }
    double *colli_time = (double*)malloc(N*(N+1)/2*sizeof(double));
    // Start sim
    for(t=0;t<S;t++)
    {
        // Step 1: speculate no collision happen, get new pos & v.
        for(i=0; i<N; i++)
        {
            particles[i].x_n = bound_pos(particles[i].x + particles[i].vx);
            particles[i].y_n = bound_pos(particles[i].y + particles[i].vy);
        }
        // Step 2: find all possible collision independently. fill colli_mat and colli_time.
        

    }


    fclose(stdin);
    free(particles);
    return 0;
}