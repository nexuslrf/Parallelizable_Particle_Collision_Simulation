#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <xmmintrin.h>
#include <string.h>

typedef struct
{
    double x;
    double y;
    double vx;
    double vy;   
}  Particle;

double doubleRand(double min, double max) // return [min, max] double vars
{
    return min+(max-min)*(rand() / (double)RAND_MAX);
}

int main()
{
    srand(0);
    int N, L, r, S;
    char mode[6];
    Particle *particles;
    freopen("./inputs.txt","r",stdin);
    scanf("%d %d %d %d %s",&N, &L, &r, &S, mode);
    particles = (Particle *)malloc(N * sizeof(Particle));
    int i =0, idx;
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
            particles[i].x = doubleRand(0.,L);
            particles[i].y = doubleRand(0.,L);
            particles[i].vx = doubleRand(L/(double)8.0/r,L/(double)4.0);
            particles[i].vy = doubleRand(L/(double)8.0/r,L/(double)4.0);
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


    fclose(stdin);
    free(particles);
    return 0;
}