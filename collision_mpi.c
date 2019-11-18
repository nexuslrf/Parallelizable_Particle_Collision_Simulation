#include <assert.h>
#include <stdio.h>
#include <stddef.h>
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
    int id;
    double x;
    double y;
    double vx;
    double vy;   
    int colli_p;
    int colli_w;
}  Particle;

typedef struct
{
    int pa;
    int pb;
    double time;
} Collision;

int n, l, r, s, bnd_far, r_sq_4, num_cmp, count;
int num_slave, myid, nprocs, slave_id, offset, send_size,
    chunk_size, last_chunk_size, dst_id;
Particle *particles, *P_a, *P_b;
Collision *colli_time, *colli;


MPI_Status Stat;
MPI_Datatype Particle_Type, Colli_Type;
MPI_Comm parti_comm;

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

void randomise_particles()
{
    /* Implement randomisation */
    int i;
    for(i=0; i<n; i++)
    {
        particles[i].id = i;
        particles[i].x = doubleRand(r,bnd_far);
        particles[i].y = doubleRand(r,bnd_far);
        particles[i].vx = (1 - 2*RAND01)*doubleRand(l/(double)8.0/r,l/(double)4.0);
        particles[i].vy = (1 - 2*RAND01)*doubleRand(l/(double)8.0/r,l/(double)4.0);
    }
}

void bound_pos(Particle *p)
{
    double tx=0,ty=0;
    if(p->x>bnd_far)
        tx = (p->x-bnd_far)/p->vx;
    else if(p->x<r)
        tx = (p->x-r)/p->vx;
    if(p->y>bnd_far)
        ty = (p->y-bnd_far)/p->vy;
    else if(p->y<r)
        ty = (p->y-r)/p->vy;

    tx =ty = tx>ty?tx:ty;
    p->x = p->x - tx*p->vx;
    p->y = p->y - ty*p->vy;
}

void check_wall_colli(int chunk_idx, int num_item)
{
    double lambda_1, lambda_2, lambda;
    int i, wall_colli;
    double x_n, y_n;
    int offset = chunk_idx * chunk_size;
    for(i=0; i<num_item; i++)
    {
        P_a = particles + i + offset;
        x_n = P_a->x + P_a->vx;
        y_n = P_a->y + P_a->vy;
        //Case 1: collision with wall
        ///////////////
        lambda_1 = lambda_2 = 2;
        wall_colli = 0;
        if(x_n<r)
        {
            lambda_1 = (r - P_a->x) / P_a->vx;
            wall_colli = 1;
        }
        else if(x_n>bnd_far)
        {
            lambda_1 = (bnd_far - P_a->x) / P_a->vx;
            wall_colli = 1;
        }

        if(y_n<r)
        {
            lambda_2 = (r - P_a->y) / P_a->vy;
            wall_colli = 1;
        }
        else if(y_n>bnd_far)
        {
            lambda_2 = (bnd_far - P_a->y) / P_a->vy;
            wall_colli = 1;
        }
        if(wall_colli)
        {
            colli_time[count].pb = P_a->id;
            lambda = lambda_1-lambda_2;
            if(fabs(lambda)<eps) // Cornor collision!
            {
                colli_time[count].pa = -1; // -1 to present this case.
                colli_time[count].time = lambda_1;
            }
            else if(lambda<0) // x wall collision!
            {
                colli_time[count].pa = -2; // -2 to present this case.
                colli_time[count].time = lambda_1;
            }
            else if(lambda>0) // y wall collision!
            {
                colli_time[count].pa = -3; // -3 to present this case.
                colli_time[count].time = lambda_2;
            }
            count++;
        }
        ///////////////
    }
}

void check_pp_colli(int chunk_idx_A, int chunk_idx_B, int num_item_A, int num_item_B)
{
    double dx1, dy1, Delta, Dx, Dy, dDpdD, dDmdD, DDpDD, lambda;
    int offset_A, offset_B, i, j;
    offset_A = chunk_idx_A * chunk_size;
    offset_B = chunk_idx_B * chunk_size;
    // printf("ID[%d] #Item_A: %d #Item_B: %d\n", myid, num_item_A, num_item_B);
    for(i=0; i<num_item_A; i++)
    {
        if(chunk_idx_A == chunk_idx_B)
            j = i + 1;
        else
            j = 0;
        for(; j<num_item_B; j++)
        {
            P_a = particles + offset_A + i;
            P_b = particles + offset_B + j;
            dx1 = P_b->x - P_a->x;
            dy1 = P_b->y - P_a->y;
            // early stop
            Dx = P_b->vx - P_a->vx;
            Dy = P_b->vy - P_a->vy;
            dDpdD = dx1*Dx + dy1*Dy;
            if(dDpdD<0) // To ensure the right direction
            {
                // Case 2: overlap at startup:
                ////////////////
                Delta = dx1*dx1 + dy1*dy1;
                if(Delta - r_sq_4<=0 && Delta!=0)
                {
                    colli_time[count].time = 0.0;
                    colli_time[count].pa = P_a->id;
                    colli_time[count].pb = P_b->id;
                    count++;
                }
                ////////////////
                else
                {
                    // Case 3: Normal collision case
                    ////////////////
                    DDpDD = Dx*Dx + Dy*Dy;
                    dDmdD = dx1*Dy - dy1*Dx;
                    Delta = r_sq_4*DDpDD - dDmdD*dDmdD;
                    if(Delta>0)
                    {
                        Delta = sqrt(Delta);
                        lambda = (-dDpdD - Delta)/DDpDD;
                        // printf("[Debug:lambda]: %f\n", lambda);
                        if(lambda<1)
                        {
                            colli_time[count].time = lambda;
                            colli_time[count].pa = P_a->id;
                            colli_time[count].pb = P_b->id;
                            count++;
                        }
                    }
                    ////////////////
                }
            }
        }
    }
}

void update_particle(Particle* parti, int num_item, int* colli_mat)
{
    int i;
    for(i=0; i<num_item; i++)
    {
        P_a = parti + i;
        if(!colli_mat[i])
        {
            P_a->x = P_a->x + P_a->vx;
            P_a->y = P_a->y + P_a->vy;
        }
        else
            colli_mat[i] = 0;
    }
}

void proc_collision(Collision* colli)
{
    double Dx ,Dy, Delta, dx1, dy1, dx2, dy2, DDpDD;
    if(colli->pa==-1) // Cornor colli;
    {
        P_a = particles + colli->pb;
        P_a->vx = -1*P_a->vx;
        P_a->vy = -1*P_a->vy;
        P_a->x = P_a->x+(1-2*colli->time)*P_a->vx;
        P_a->y = P_a->y+(1-2*colli->time)*P_a->vy;
        P_a->colli_w++;
        bound_pos(P_a);
    }
    else if(colli->pa==-2)//  X wall colli;
    {
        P_a = particles + colli->pb;
        P_a->vx = -1*P_a->vx;
        P_a->x = P_a->x+(1-2*colli->time)*P_a->vx;
        P_a->y = P_a->y+P_a->vy;
        P_a->colli_w++;
        bound_pos(P_a);
    }
    else if(colli->pa==-3)// Y wall colli;
    {
        P_a = particles + colli->pb;
        P_a->vy = -1*P_a->vy;
        P_a->y = P_a->y+(1-2*colli->time)*P_a->vy;
        P_a->x = P_a->x+P_a->vx;
        P_a->colli_w++;
        bound_pos(P_a);
    }
    else // P-P colli;
    {
        P_a = particles + colli->pa;
        P_b = particles + colli->pb;
        P_a->colli_p++;
        P_b->colli_p++;
        P_a->x = P_a->x + colli->time*P_a->vx;
        P_a->y = P_a->y + colli->time*P_a->vy;
        P_b->x = P_b->x + colli->time*P_b->vx;
        P_b->y = P_b->y + colli->time*P_b->vy;
        Dx = P_b->x - P_a->x;
        Dy = P_b->y - P_a->y;
        Delta = 1 - colli->time;
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
        P_a->x = P_a->x + Delta*P_a->vx;
        P_a->y = P_a->y + Delta*P_a->vy;
        bound_pos(P_a);
        P_b->x = P_b->x + Delta*P_b->vx;
        P_b->y = P_b->y + Delta*P_b->vy;
        bound_pos(P_b);
    }

}

/*
Data distribution policy
vars:
    send_mat: each row i represent the chunk needed by slave i, 
        the last item of each row is the #item needed to be sent;
    job_mat[i][j]: indicate the collision between chunk i & j is 
        processed by slave job_mat[i][j]; 
*/
int gen_comm_mat(int (*send_mat)[nprocs], int (*job_mat)[num_slave], int* pa_idx, int* pb_idx)
{
    int i, j, k, proc, cnt=0, job_cnt=0;
    int flag_mat[num_slave][num_slave];
    memset(flag_mat[0] ,0, num_slave*num_slave*sizeof(int));
    for(k=num_slave; k>0; k--)
    {
        j = num_slave - k;
        for(i=0; i<k; i++)
        {
            proc = cnt % num_slave;
            job_mat[i][j] = proc;
            if(!flag_mat[proc][i])
            {
                flag_mat[proc][i] = 1;
                send_mat[proc][send_mat[proc][num_slave]] = i;
                send_mat[proc][num_slave]++;
            }
            if(!flag_mat[proc][j])
            {
                flag_mat[proc][j] = 1;
                send_mat[proc][send_mat[proc][num_slave]] = j;
                send_mat[proc][num_slave]++;
            }
            if(proc+1==myid)
            {
                pa_idx[job_cnt] = i;
                pb_idx[job_cnt] = j;
                job_cnt ++;
            }
            cnt++;
            j++;
        }
    }
    return job_cnt;
}

void gen_parti_type(MPI_Datatype *Type_ptr)
{
    MPI_Datatype type[7] = { 
        MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
        MPI_INT, MPI_INT };
    int blocklen[7] = { 1, 1, 1, 1, 1, 1, 1 };
    MPI_Aint disp[7];
    disp[0] = offsetof(Particle, id);
    disp[1] = offsetof(Particle, x);
    disp[2] = offsetof(Particle, y);
    disp[3] = offsetof(Particle, vx);
    disp[4] = offsetof(Particle, vy);
    disp[5] = offsetof(Particle, colli_p);
    disp[6] = offsetof(Particle, colli_w);
    MPI_Type_create_struct(7, blocklen, disp, type, Type_ptr);
    MPI_Type_commit(Type_ptr);
}

void gen_colli_type(MPI_Datatype *Type_ptr)
{
    MPI_Datatype type[3] = { MPI_INT, MPI_INT, MPI_DOUBLE };
    int blocklen[3] = { 1, 1, 1 };
    MPI_Aint disp[3];

    disp[0] = offsetof(Collision, pa);
    disp[1] = offsetof(Collision, pb);
    disp[2] = offsetof(Collision, time);
    MPI_Type_create_struct(3, blocklen, disp, type, Type_ptr);
    MPI_Type_commit(Type_ptr);
}

/*
Customized MPI communicator
*/
void gen_mpi_comm(int num_proc, MPI_Comm *new_comm)
{
    int i, group[num_proc];
    for(i=0; i<num_proc; i++)
        group[i] = i;
    // Get the group or processes of the default communicator
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    MPI_Group new_group;
    MPI_Group_incl(world_group, num_proc, group, &new_group);
    // Create the new communicator from that group of processes.
    MPI_Comm new_communicator;
    MPI_Comm_create(MPI_COMM_WORLD, new_group, new_comm);
}

void print_particles(int step)
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%d %d %10.8lf %10.8lf %10.8lf %10.8lf\n", step, i, particles[i].x, particles[i].y,
            particles[i].vx, particles[i].vy);
    }
}

void print_statistics(int num_step)
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%d %d %10.8lf %10.8lf %10.8lf %10.8lf %d %d\n", num_step, i, particles[i].x,
            particles[i].y, particles[i].vx, particles[i].vy,
            particles[i].colli_p, particles[i].colli_w);
    }
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
    char mode_buf[6];
    scanf("%d", &n);
    scanf("%d", &l);
    scanf("%d", &r);
    scanf("%d", &s);
    scanf("%5s", mode_buf);
    
    // Use collective func.
    MPI_Bcast(&n, 1, MPI_INT, MASTER_ID, parti_comm);
    MPI_Bcast(&l, 1, MPI_INT, MASTER_ID, parti_comm);
    MPI_Bcast(&r, 1, MPI_INT, MASTER_ID, parti_comm);
    MPI_Bcast(&s, 1, MPI_INT, MASTER_ID, parti_comm);
    ///////////
    if(n<num_slave)
    {
        num_slave = n;
        nprocs = n+1;
        gen_mpi_comm(nprocs, &parti_comm);
    }
    int i, j, k, m, step, real_colli, extra_cnt;
    double x, y, vx, vy;
    simulation_mode_t mode;
    int send_mat[num_slave][nprocs], job_mat[num_slave][num_slave];
    int countbuf[nprocs];
    Particle parti_buf;
    memset(send_mat[0], 0, num_slave*nprocs*sizeof(int));
    memset(job_mat[0], 0, num_slave*num_slave*sizeof(int));
    gen_comm_mat(send_mat, job_mat, NULL, NULL);

    bnd_far = l - r;
    num_cmp = n * (n+1) / 2;
    chunk_size = (n-1) / num_slave + 1;
    last_chunk_size = n - (num_slave - 1) * chunk_size;
    particles = (Particle *)malloc(n * sizeof(Particle));
    int colli_chunk_queue[num_slave][n], colli_mat[n];
    j = 0;
    while (scanf("%d %lf %lf %lf %lf", &i, &x, &y, &vx, &vy) != EOF) {
        j++;
        particles[i].id = i;
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
    
    colli_time = (Collision *)malloc(num_cmp *sizeof(Collision));
    print_particles(0);
    /*
        start sim
    */
    for(step = 0; step < s; step++)
    {
        // Scatter particles
        for(k=0;k<num_slave;k++)
        {
            for(i=0;i<send_mat[k][num_slave];i++)
            {
                offset = send_mat[k][i] * chunk_size;
                send_size = send_mat[k][i]==num_slave-1?last_chunk_size:chunk_size;
                MPI_Send(particles+offset, send_size, Particle_Type, k+1, i, parti_comm);
            }
        }
        memset(colli_mat, 0, n*sizeof(int));
        // printf("ID[%d] Start gathering!\n", myid);
        MPI_Gather(&count, 1, MPI_INT, countbuf, 1, MPI_INT, MASTER_ID, parti_comm);
        // printf("ID[%d] Get countbuf!\n", myid);
        count = 0;
        extra_cnt = 0;
        for(i=1; i<num_slave+1; i++)
        {
            MPI_Recv(colli_time+count, countbuf[i], Colli_Type, i, i, parti_comm, &Stat);
            count += countbuf[i];
        }
        // printf("ID[%d] Get colli_time!\n", myid);
        qsort(colli_time, count, sizeof(Collision), compare);
        memset(countbuf, 0, nprocs*sizeof(int));
        // printf("ID[%d] Finding Real_colli, #Colli: %d!\n", myid, count);
        for(k=0;k<count;k++)
        {
            colli = colli_time+k;
            ////
            // if(step==0) // && (colli->pa == 0||colli->pb==0))
            // {
            //     printf("[Debug:inconsist] %d %d %10.8f\n",colli->pa, colli->pb, colli->time);
            // }
            ////
            if(colli->pa<0)
            { //wall collision
                if(!colli_mat[colli->pb])
                {
                    colli_mat[colli->pb] = 1;
                    i = colli->pb / chunk_size;
                    colli_chunk_queue[i][countbuf[i + 1]++]=k;
                }
            }
            else if(!colli_mat[colli->pa]) // p-p collision
            {
                if(!colli_mat[colli->pb])
                {
                    colli_mat[colli->pa] = 1;
                    colli_mat[colli->pb] = 1;
                    i = colli->pa / chunk_size;
                    j = colli->pb / chunk_size;
                    m = job_mat[i][j];
                    colli_chunk_queue[m][countbuf[m + 1]++]=k;
                    // if (i!=j)
                    //     colli_chunk_queue[j][countbuf[j + 1]++]=k;
                }
            }
        }
        // printf("ID[%d] Scatter colli_time!\n", myid);
        MPI_Scatter(countbuf, 1, MPI_INT, &count, 1, MPI_INT, MASTER_ID, parti_comm);
        for(i=0; i<num_slave; i++)
        {
            for(j=0; j<countbuf[i+1]; j++)
            {
                k = colli_chunk_queue[i][j];
                MPI_Send(colli_time+k, 1, Colli_Type, i+1, i+1, parti_comm);
            }
        }
        // printf("ID[%d] Gather extra_cnt!\n", myid);
        MPI_Gather(&extra_cnt, 1, MPI_INT, countbuf, 1, MPI_INT, MASTER_ID, parti_comm);
        offset = 0;
        // printf("ID[%d] Gather major chunk!\n", myid);
        for(i=0;i<num_slave;i++)
        {
            send_size = i!=num_slave-1?chunk_size:last_chunk_size;
            MPI_Recv(particles+offset, send_size, Particle_Type, i+1, i+1, parti_comm, &Stat);
            offset += chunk_size;
        }
        // printf("ID[%d] Gather extra parti!\n", myid);
        for(i=0;i<num_slave;i++)
        {
            for(j=0; j<countbuf[i+1]; j++)
            {
                MPI_Recv(&parti_buf, 1, Particle_Type, i+1, i+1, parti_comm, &Stat);
                memcpy(particles+parti_buf.id, &parti_buf, sizeof(Particle));
            }
        }
        if(mode==MODE_PRINT)
            print_particles(step+1);
    }
    print_statistics(s);
    double exec_time=GetTimer();
    printf("Time elapsed: %lf ms",exec_time);
    free(particles);
    free(colli_time);
}

void slave()
{
    MPI_Bcast(&n, 1, MPI_INT, MASTER_ID, parti_comm);
    MPI_Bcast(&l, 1, MPI_INT, MASTER_ID, parti_comm);
    MPI_Bcast(&r, 1, MPI_INT, MASTER_ID, parti_comm);
    MPI_Bcast(&s, 1, MPI_INT, MASTER_ID, parti_comm);
    if(n<num_slave)
    {
        num_slave = n;
        nprocs = n+1;
        gen_mpi_comm(nprocs, &parti_comm);
        if(slave_id>=n)
            return;
    }
    int i, j, k, step, num_chunk, major_chunk, num_chunk_cmp, extra_cnt;
    num_chunk_cmp = num_slave * (num_slave + 1) / 2;
    num_chunk_cmp = (num_chunk_cmp-1) / num_slave +1;
    int chunk_map[num_slave], chunk_size_list[num_slave],
        pa_idx[num_chunk_cmp], pb_idx[num_chunk_cmp];
    int send_mat[num_slave][nprocs], job_mat[num_slave][num_slave];
    int countbuf[nprocs];
    memset(send_mat[0], 0, num_slave*nprocs*sizeof(int));
    memset(job_mat[0], 0, num_slave*num_slave*sizeof(int));
    memset(chunk_map, -1, num_slave*sizeof(int));

    num_chunk_cmp = gen_comm_mat(send_mat, job_mat, pa_idx, pb_idx);
    bnd_far = l - r;
    r_sq_4 = r * r * 4;
    chunk_size = (n-1) / num_slave + 1;
    num_chunk = send_mat[slave_id][num_slave];
    last_chunk_size = n - (num_slave - 1) * chunk_size;
    n = (num_chunk - 1) * chunk_size + last_chunk_size;
    num_cmp = n * (n+1) / 2;
    particles = (Particle *)malloc(num_chunk * chunk_size * sizeof(Particle));
    colli_time = (Collision *)malloc(num_cmp *sizeof(Collision));
    int colli_mat[chunk_size];
    memset(colli_mat, 0, chunk_size * sizeof(int));
    int extra_send[chunk_size * (num_chunk - 1)];
    /*
        start sim
    */
    for(step = 0; step < s; step++)
    {
        // Scatter particles
        for(i=0;i<num_chunk;i++)
        {
            offset = chunk_size * i;
            send_size = send_mat[slave_id][i]==num_slave-1?last_chunk_size:chunk_size;
            MPI_Recv(particles+offset, send_size, Particle_Type, MASTER_ID, i, parti_comm, &Stat);
            if(send_mat[slave_id][i]==slave_id)
            {
                major_chunk = i;
            }
            chunk_map[send_mat[slave_id][i]] = i;
            chunk_size_list[i] = send_size;
        }
        // 
        // printf("ID[%d] Start Detecting!\n", myid);
        count = 0;
        extra_cnt = 0;
        check_wall_colli(major_chunk, chunk_size_list[major_chunk]);
        // printf("ID[%d] Start P-P Detecting!\n", myid);
        for(k=0; k<num_chunk_cmp; k++)
        {
            i = chunk_map[pa_idx[k]];
            j = chunk_map[pb_idx[k]];
            check_pp_colli(i, j, chunk_size_list[i], chunk_size_list[j]);
        }
        // for(i = 0; i<count; i++)
        // {
        //     printf("ID[%d] %d %d %10.8f\n", myid, colli_time[i].pa, colli_time[i].pb, colli_time[i].time);
        // }
        // printf("ID[%d] Start gathering! #Colli: %d\n", myid, count);
        MPI_Gather(&count, 1, MPI_INT, countbuf, 1, MPI_INT, MASTER_ID, parti_comm);
        MPI_Send(colli_time, count, Colli_Type, MASTER_ID, myid, parti_comm);
        MPI_Scatter(countbuf, 1, MPI_INT, &count, 1, MPI_INT, MASTER_ID, parti_comm);
        for(i=0; i<count; i++)
        {
            MPI_Recv(colli_time+i, 1, Colli_Type, MASTER_ID, myid, parti_comm, &Stat);
            // if(step==0)
            //     printf("ID[%d] Colli[%d] %d %d %10.8f\n", myid, i, colli_time[i].pa, colli_time[i].pb, colli_time[i].time);
        }
        for(k=0; k<count; k++)
        {
            colli = colli_time + k;
            if(colli->pa<0) // wall
            {
                // printf("ID[%d] Colli Wall[%d]: %d\n", myid, k, colli->pb);
                colli->pb = chunk_map[colli->pb / chunk_size] * chunk_size + colli->pb % chunk_size;
                colli_mat[colli->pb % chunk_size]=1;
                // printf("ID[%d] Colli Wall_X[%d]: %d \n", myid, k, (particles + colli->pb)->id);
                proc_collision(colli);
            }
            else
            {
                // printf("ID[%d] Colli P-P[%d]: %d %d\n", myid, k, colli->pa, colli->pb);
                i = colli->pa / chunk_size;
                j = colli->pb / chunk_size;
                if(chunk_map[i]>=0 && chunk_map[j]>=0)
                {
                    colli->pa = chunk_map[i] * chunk_size + colli->pa % chunk_size;
                    colli->pb = chunk_map[j] * chunk_size + colli->pb % chunk_size;
                    // printf("ID[%d] Colli P-P_X[%d]: %d %d\n", myid, k, (particles + colli->pa)->id, (particles + colli->pb)->id);
                    if(i==slave_id)
                        colli_mat[colli->pa % chunk_size]=1;
                    else
                        extra_send[extra_cnt++] = colli->pa;
                    if(j==slave_id)
                        colli_mat[colli->pb % chunk_size]=1;
                    else
                        extra_send[extra_cnt++] = colli->pb;

                    proc_collision(colli);
                }
            }
        }
        update_particle(particles+major_chunk*chunk_size, chunk_size_list[major_chunk], colli_mat);
        // for (i = major_chunk*chunk_size; i < major_chunk*chunk_size + chunk_size_list[major_chunk]; i++) 
        // {
        //     printf("ID[%d] Major: %d %d %10.8lf %10.8lf %10.8lf %10.8lf\n", myid, step, particles[i].id, particles[i].x, particles[i].y,
        //         particles[i].vx, particles[i].vy);
        // }
        MPI_Gather(&extra_cnt, 1, MPI_INT, countbuf, 1, MPI_INT, MASTER_ID, parti_comm);
        offset = major_chunk * chunk_size;
        MPI_Send(particles+offset, chunk_size_list[major_chunk], Particle_Type, MASTER_ID, myid, parti_comm);
        for(i=0; i<extra_cnt; i++)
        {
            offset = extra_send[i];
            MPI_Send(particles+offset, 1, Particle_Type, MASTER_ID, myid, parti_comm);
            // printf("ID[%d] Extra: %d %d %10.8lf %10.8lf %10.8lf %10.8lf\n", myid, step, 
            //     particles[offset].id, particles[offset].x, particles[offset].y, particles[offset].vx, particles[offset].vy);
        }
    }
    free(particles);
    free(colli_time);
}


int main(int argc, char ** argv)
{
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (nprocs < 2)
    {
        fprintf(stderr, "#Proc >= 2 !\n");
        return 1;
    }
    parti_comm = MPI_COMM_WORLD;
    num_slave = nprocs - 1;
    slave_id = myid - 1;
    gen_parti_type(&Particle_Type);
    gen_colli_type(&Colli_Type);
    /*
    Master-Slave Pattern
    */
    if (myid == MASTER_ID)
    {
        master();
    }
    else
    {
        slave();
    }
    // printf("ID[%d] Done!\n", myid); 
    MPI_Finalize();
    return 0;
}