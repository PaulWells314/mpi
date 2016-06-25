#include <mpi.h>
#include <stdio.h>
#define SLICE_WIDTH 64
#define OVERLAP 8
#define NUM_SLICES 2

float un[SLICE_WIDTH * NUM_SLICES];
float un_1[SLICE_WIDTH * NUM_SLICES];

float local_data[SLICE_WIDTH];
float local_data_1[SLICE_WIDTH];



void init_signal(float* pu)
{
  int i;
  for (i=0; i< SLICE_WIDTH  * NUM_SLICES - 2 * OVERLAP * (NUM_SLICES-1) ; i++)
  {
     *(pu+i) = 0.0f;
     if (i > 80 && i< 100)
     {
         *(pu+i) = 1.0f;
     }
  }
}

void display_signal(float* pu)
{
  int i;
  for (i=0; i< SLICE_WIDTH * NUM_SLICES - 2 * OVERLAP * (NUM_SLICES-1); i++)
  {
     printf("%d %f \n", i, *(pu + i));
  }
}

void update_signal(float*p1, float*p2,  int start, int finish)
{
   int i;
   for (i=start; i< finish; i++)
   {
      *(p2+i) = *(p1+i) + 0.01f * ( *(p1+i+1) -2.0f * *(p1+i) + *(p1+i-1));
   }
}


int main(int argc, char **argv) {
    int my_rank;
    float * p1;
    float * p2;
    float * ptemp;
    int i, j;
    MPI_Status  stats[NUM_SLICES];
    MPI_Request requests[NUM_SLICES];
    int ierr;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); /*DETERMINE RANK OF THIS PROCESSOR*/
    printf("Hello Paul %d\n", my_rank);

    if (my_rank == 0)
    {
        init_signal(un);
        p1 = un;    
        for (j = 0; j < SLICE_WIDTH; j++)
        {
            local_data[j] = *(un+j);
        }

      
        for (i= 1; i < NUM_SLICES; i++)
        {
            MPI_Isend(un + (SLICE_WIDTH-2* OVERLAP) * i, SLICE_WIDTH, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &requests[i]);
        }
        MPI_Waitall(NUM_SLICES-1, &requests[1], &stats[1]);

        for (j=0; j< 8192; j++)
{

        for (i=0; i < OVERLAP/2; i++)
        {
            update_signal(local_data, local_data_1, 0, SLICE_WIDTH);
        }

        for (i=0; i < OVERLAP/2; i++)
        {
            update_signal(local_data_1, local_data, 0, SLICE_WIDTH);
        }
        MPI_Isend(local_data + (SLICE_WIDTH- 2 *OVERLAP), OVERLAP, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(local_data + (SLICE_WIDTH - OVERLAP), OVERLAP, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &requests[1]);
        MPI_Waitall(2, &requests[0], &stats[0]);
}

        for (i= 1; i < NUM_SLICES; i++)
        {
            MPI_Irecv(un_1 + (SLICE_WIDTH-2* OVERLAP) * i, SLICE_WIDTH, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &requests[i]);
        }
        MPI_Waitall(NUM_SLICES-1, &requests[1], &stats[1]);

        for (j = 0; j < SLICE_WIDTH; j++)
        {
             *(un_1 +j) = local_data[j];
        }
        display_signal(un_1);
    }
    else
    {
        MPI_Irecv(local_data, SLICE_WIDTH, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Wait(&requests[0], &stats[0]);

       for (j=0; j< 8192; j++)
{
        for (i=0; i < OVERLAP/2; i++)
        {
            update_signal(local_data, local_data_1, 0, SLICE_WIDTH);
        }

        for (i=0; i < OVERLAP/2; i++)
        {
            update_signal(local_data_1, local_data, 0, SLICE_WIDTH);
        }
        MPI_Isend(local_data + OVERLAP, OVERLAP, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(local_data, OVERLAP, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Waitall(2, &requests[0], &stats[0]);
}
        MPI_Isend(local_data, SLICE_WIDTH, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Wait(&requests[0], &stats[0]);
    }

    char processor_name[128];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s\n",
           processor_name);

    MPI_Finalize();
}