#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <xmmintrin.h>
#include <string.h>
#include <math.h>

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
typedef struct
{
    int pa;
    int pb;
    double time;
} Collision;
int compare (const void * a, const void * b)
{
    Collision *colli_A = a;
    Collision *colli_B = b;
    double cmpf = colli_B->time - colli_A->time;
    if(cmpf!=0)
        return cmpf<0?-1:1;
    else
    {
        int cmpt = colli_B->pa - colli_A->pa;
        if(cmpt!=0)
            return cmpt;
        else
            return colli_B->pb - colli_A->pb;
    }
}
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
    Particle *particles, *P_a, *P_b;
    freopen("./inputs.txt","r",stdin);
    scanf("%d %d %d %d %s",&N, &L, &r, &S, mode);
    particles = (Particle *)malloc(N * sizeof(Particle));
    int i =0,j, t, idx, cnt;
    double x, y, vx, vy, lambda, lambda_1, lambda_2;
    double dx1, dx2, dy1, dy2, Dx, Dy, DDpDD, dDpdD, dDmdD, Delta;
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

    // int **colli_mat;
    // colli_mat = (int**)malloc((N+1)*sizeof(int*));
    // for(i=0; i<N+1; i++)
    // {    
    //     colli_mat[i]=(int*)malloc((N+1)*sizeof(int));
    //     memset(colli_mat[i],0,(N+1)*sizeof(int));
    // }
    int *colli_mat = (int *)malloc(N*sizeof(int));
    memset(colli_mat,0, N*sizeof(int));
    Collision *colli_time = (Collision*)malloc(N*(N+1)/2*sizeof(Collision));
    // Start sim
    double r_sq_4 = 4*r*r;
    for(t=0;t<S;t++)
    {
        // Step 1: speculate no collision happen, get new pos & v.
        for(i=0; i<N; i++)
        {
            particles[i].x_n = particles[i].x + particles[i].vx;
            particles[i].y_n = particles[i].y + particles[i].vy;
        }
        // Step 2: find all possible collision independently. fill colli_mat and colli_time.
        cnt = 0;
        for(i=0; i<N; i++)
        {
            P_a = particles+i;
            //Case 1: collision with wall
            ///////////////
            lambda_1 = lambda_2 = 2;
            if(P_a->x_n<r)
            {
                lambda_1 = (r - P_a->x) / P_a->vx;
            }
            else if(P_a->x_n>L-r)
            {
                lambda_1 = (L-r - P_a->x) / P_a->vx;
            }

            if(P_a->y_n<r)
            {
                lambda_2 = (r - P_a->y) / P_a->vy;
            }
            else if(P_a->y_n>L-r)
            {
                lambda_2 = (L-r - P_a->y) / P_a->vy;
            }
            lambda = lambda_1 < lambda_2? lambda_1:lambda_2;
            if(lambda < 1)
            {
                colli_time[cnt].time = lambda;
                colli_time[cnt].pa = i;
                colli_time[cnt].pb = N;
                cnt++;
            }
            ///////////////
            for(j=i+1; i<N; i++)
            {
                P_b = particles+j;
                // Case 2: overlap at startup
                ////////////////
                dx1 = P_b->x - P_a->x;
                dy1 = P_b->y - P_a->y;
                if(dx1*dx1 + dy1*dy1 - r_sq_4<=0)
                {
                    colli_time[cnt].time = 0;
                    colli_time[cnt].pa = i;
                    colli_time[cnt].pb = j;
                    cnt++;
                    break; // no need to further detect.
                }
                ////////////////
                // Case 3: Normal collision case
                ////////////////
                dx2 = P_b->x_n - P_a->x_n;
                dy2 = P_b->y_n - P_a->y_n;
                Dx = dx2 - dx1;
                Dy = dy2 - dy1;
                DDpDD = Dx*Dx + Dy*Dy;
                dDmdD = dx1*Dy - dy1*Dx;
                Delta = r_sq_4*DDpDD - dDmdD*dDmdD;
                dDpdD = dx1*Dx + dy1*Dy;
                if(Delta<=0||dDpdD>0)
                    continue;
                Delta = sqrt(Delta);
                lambda_1 = (-dDpdD + Delta)/DDpDD;
                lambda_2 = (-dDpdD - Delta)/DDpDD;
                lambda = lambda_1<lambda_2?lambda_1:lambda_2;
                if(lambda<1)
                {
                    colli_time[cnt].time = lambda;
                    colli_time[cnt].pa = i;
                    colli_time[cnt].pb = j;
                    cnt++;
                }
                ////////////////
            }
        }
        // Step 3: sort collision table and process collision
        

    }


    fclose(stdin);
    free(particles);
    free(colli_time);
    free(colli_mat);
    return 0;
}