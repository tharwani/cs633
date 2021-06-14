
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>


// calculates yearwise minimum on local data array.
void find_yearwise_min(float* minarr, float* data, int num_rows, int num_cols){
	int idx = 0;
	for(int x = 0; x < num_rows; x++){
		for(int y = 0; y < num_cols; y++){
			if(minarr[y] > data[idx + y]){
				minarr[y] = data[idx + y];
			}
		}
		idx += num_cols;
	}
}

// given an array of yearwise minimums, calculates the global minimum
float find_minimum(float* input, int num){
    float ans = 1000000000000;
    for(int y = 0; y < num; y++){
			if(ans > input[y]){
				ans = input[y];
			}
	}
    return ans;
}

// reads the file using sequential io and sets num_col 
//with the actual number of columns and mat array with actual data 
// excluding lattitude and longitude and returns total number of rows
int readfile(float* mat, int* num_col, char* filename){
	char line[1000];

	FILE *fptr;
	fptr = fopen(filename,"r");

	fscanf(fptr,"%s",line);
	int c = 0;
	for(int count = 0; count < strlen(line); count++){
		if(line[count] == ','){
			c++;
		}
	}
	c = c-1;
	*num_col = c;
	char comma;
	fscanf(fptr,"%c",&comma);
	int idx = 0;
	while(1){
		if(fscanf(fptr,"%f",&mat[idx]) == EOF) break; // scans the value of latitude
		fscanf(fptr,"%c",&comma);
		fscanf(fptr,"%f",&mat[idx+1]);
		for(int x = 0; x < c; x++){
			fscanf(fptr,"%c",&comma);
			fscanf(fptr,"%f",&mat[idx+x]);	
		}
        idx += c;
		fscanf(fptr,"%c",&comma);
	}
	return idx/c;
}

int main(int argc, char *argv[])
{
	FILE *fptr;
	fptr = fopen("output.txt","w");
	char* filename = argv[1];
	// printf("%s\n",argv[0]);
	printf("%s\n",filename);
	int myrank, size;
	MPI_Status status, sstatus;
	MPI_Request request;
	double time,time2;
	float * mat = (float *) malloc(sizeof(float)*1024*1024*150) ;
	int num_rows;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// read the csv file from one process and a barrier before starting the timer.
	int num_cols;
	if(myrank==0) num_rows  = readfile(mat,&num_cols,filename);
	MPI_Barrier(MPI_COMM_WORLD);
	time = MPI_Wtime();
	
	// sharing number of rows and columns across processes
	MPI_Bcast( &num_rows , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
	MPI_Bcast( &num_cols , 1 , MPI_INT , 0 , MPI_COMM_WORLD);

	float* new_mat;
	int *sendcounts = (int*)malloc(sizeof(int)*size);
	int *displ = (int*)malloc(sizeof(int)*(size+1));
	int base = num_rows/size;
	int l = num_rows%size;
	displ[0] = 0;
	for(int i=0;i<size;i++){
		if(i < l){
			sendcounts[i] = (base + 1)*num_cols;
		}
		else{
			sendcounts[i] = base*num_cols;
		}
		displ[i+1] = displ[i] + sendcounts[i];
	}
	new_mat = (float *) malloc(sizeof(float)*sendcounts[myrank]) ;

	// We scatter data such that the data is equally distributed across all proccesses. If not possible, number of rows to be processed by each process differs by one
	MPI_Scatterv( mat , sendcounts , displ , MPI_FLOAT , new_mat , sendcounts[myrank] , MPI_FLOAT , 0 , MPI_COMM_WORLD);
	float minimum_yearwise[num_cols];
	for(int x = 0; x < num_cols; x++){
		minimum_yearwise[x] = 1000000000.11; // max val
	}
	find_yearwise_min(minimum_yearwise,new_mat,sendcounts[myrank]/num_cols,num_cols);

	// given local values of yearwise minimum, reduce is used to get yerawise local minimum across all processes.
	if(myrank==0) MPI_Reduce( MPI_IN_PLACE , minimum_yearwise , num_cols , MPI_FLOAT , MPI_MIN , 0 , MPI_COMM_WORLD);
	else MPI_Reduce( minimum_yearwise , minimum_yearwise , num_cols , MPI_FLOAT , MPI_MIN , 0 , MPI_COMM_WORLD);

	// global minimum calculated
	if(myrank==0){
		float global_min = find_minimum(minimum_yearwise,num_cols);
		for(int i=0;i<num_cols;i++){
			fprintf(fptr,"%f",minimum_yearwise[i]);
			if(i<num_cols-1) fprintf(fptr,","); 
		}
		fprintf(fptr,"\n");
		fprintf(fptr,"%f\n",global_min);
	}
	time = MPI_Wtime() - time;
	double finaltime = 0;
	MPI_Reduce( &time , &finaltime , 1 , MPI_DOUBLE , MPI_MAX , 0 , MPI_COMM_WORLD);
	if(myrank==0) fprintf(fptr,"%lf\n",finaltime);
	MPI_Finalize();
	return 0;
}
