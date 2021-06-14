#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
// tags
#define LEFT 0
#define RIGHT 1
#define TOP 2
#define BOT 3



// function responsible for stencil computation in data array
void run_stencil_computation(int flag, int datapoints_per_row,double ** receivebuf,int left,int right,int top,int bot,double *** data){
    int flag_new = flag^1;
    // the comments below show the region of array the stencil computation is performed. "X" marks all such areas.
    /*
        oooooo
        oXXXXo
        oXXXXo
        oXXXXo
        oXXXXo
        oooooo
    */
    for(int i = 1; i < datapoints_per_row-1; i++){
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][i][j] = (data[flag][i-1][j] + data[flag][i+1][j] + data[flag][i][j+1] + data[flag][i][j-1])/4;
        }
    }

    /*
        oXXXXo
        oooooo
        oooooo
        oooooo
        oooooo
        oooooo
    */
    if(top != -1){
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][0][j] = (receivebuf[TOP][j] + data[flag][1][j] + data[flag][0][j+1] + data[flag][0][j-1])/4;
        }
    }
    else{
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][0][j] = ( data[flag][1][j] + data[flag][0][j+1] + data[flag][0][j-1])/3;
        }
    }

    /*
        oooooo
        oooooo
        oooooo
        oooooo
        oooooo
        oXXXXo
    */
    if(bot != -1){
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][datapoints_per_row-1][j] = (receivebuf[BOT][j] + data[flag][datapoints_per_row-2][j] + data[flag][datapoints_per_row-1][j+1] + data[flag][datapoints_per_row-1][j-1])/4;
        }
    }
    else{
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][datapoints_per_row-1][j] = (data[flag][datapoints_per_row-2][j] + data[flag][datapoints_per_row-1][j+1] + data[flag][datapoints_per_row-1][j-1])/3;
        }
    }
    
    /*
        oooooo
        Xooooo
        Xooooo
        Xooooo
        Xooooo
        oooooo
    */
    if(left != -1){
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][j][0] = (receivebuf[LEFT][j] + data[flag][j][1] + data[flag][j+1][0] + data[flag][j-1][0])/4;
        }
    }
    else{
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][j][0] = (data[flag][j][1] + data[flag][j+1][0] + data[flag][j-1][0])/3;
        }
    }

    /*
        oooooo
        oooooX
        oooooX
        oooooX
        oooooX
        oooooo
    */
    if(right != -1){
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][j][datapoints_per_row-1] = (receivebuf[RIGHT][j] + data[flag][j][datapoints_per_row-2] + data[flag][j+1][datapoints_per_row-1] + data[flag][j-1][datapoints_per_row-1])/4;
        }
    }
    else{
        for(int j = 1; j < datapoints_per_row-1; j++){
            data[flag_new][j][datapoints_per_row-1] = (data[flag][j][datapoints_per_row-2] + data[flag][j+1][datapoints_per_row-1] + data[flag][j-1][datapoints_per_row-1])/3;
        }
    }

    /*
        Xooooo
        oooooo
        oooooo
        oooooo
        oooooo
        oooooo
    */
    if(left != -1 && top != -1){
        data[flag_new][0][0] = (receivebuf[TOP][0] + data[flag][1][0] + data[flag][0][1] + receivebuf[LEFT][0])/4;
    }
    else if(left != -1){
        data[flag_new][0][0] = ( data[flag][1][0] + data[flag][0][1] + receivebuf[LEFT][0])/3;
    }
    else if(top != -1){
        data[flag_new][0][0] = (receivebuf[TOP][0] + data[flag][1][0] + data[flag][0][1])/3;
    }
    else{
        data[flag_new][0][0] = (data[flag][1][0] + data[flag][0][1])/2;
    }

    /*
        oooooX
        oooooo
        oooooo
        oooooo
        oooooo
        oooooo
    */
    if(right != -1 && top != -1){
        data[flag_new][0][datapoints_per_row-1] = (receivebuf[TOP][datapoints_per_row-1] + data[flag][1][datapoints_per_row-1] + data[flag][0][datapoints_per_row-2] + receivebuf[RIGHT][0])/4;
    }
    else if(right != -1){
        data[flag_new][0][datapoints_per_row-1] = ( data[flag][1][datapoints_per_row-1] + data[flag][0][datapoints_per_row-2] + receivebuf[RIGHT][0])/3;
    }
    else if(top != -1){
        data[flag_new][0][datapoints_per_row-1] = (receivebuf[TOP][datapoints_per_row-1] + data[flag][1][datapoints_per_row-1] + data[flag][0][datapoints_per_row-2])/3;
    }
    else{
        data[flag_new][0][datapoints_per_row-1] = (data[flag][1][datapoints_per_row-1] + data[flag][0][datapoints_per_row-2])/2;
    }

     /*
        oooooo
        oooooo
        oooooo
        oooooo
        oooooo
        Xooooo
    */
    if(left != -1 && bot != -1){
        data[flag_new][datapoints_per_row-1][0] = (receivebuf[BOT][0] + data[flag][datapoints_per_row-2][0] + data[flag][datapoints_per_row-1][1] + receivebuf[LEFT][datapoints_per_row-1])/4;
    }
    else if(left != -1){
        data[flag_new][datapoints_per_row-1][0] = (data[flag][datapoints_per_row-2][0] + data[flag][datapoints_per_row-1][1] + receivebuf[LEFT][datapoints_per_row-1])/3;
    }
    else if(bot != -1){
        data[flag_new][datapoints_per_row-1][0] = (receivebuf[BOT][0] + data[flag][datapoints_per_row-2][0] + data[flag][datapoints_per_row-1][1])/3;
    }
    else{
        data[flag_new][datapoints_per_row-1][0] = (data[flag][datapoints_per_row-2][0] + data[flag][datapoints_per_row-1][1])/2;
    }

    /*
        oooooo
        oooooo
        oooooo
        oooooo
        oooooo
        oooooX
    */
    if(right != -1 && bot != -1){
        data[flag_new][datapoints_per_row-1][datapoints_per_row-1] = (receivebuf[BOT][datapoints_per_row-1] + data[flag][datapoints_per_row-2][datapoints_per_row-1] + data[flag][datapoints_per_row-1][datapoints_per_row-2] + receivebuf[RIGHT][datapoints_per_row-1])/4;
    }
    else if(right != -1){
        data[flag_new][datapoints_per_row-1][datapoints_per_row-1] = (data[flag][datapoints_per_row-2][datapoints_per_row-1] + data[flag][datapoints_per_row-1][datapoints_per_row-2] + receivebuf[RIGHT][datapoints_per_row-1])/3;
    }
    else if(bot != -1){
        data[flag_new][datapoints_per_row-1][datapoints_per_row-1] = (receivebuf[BOT][datapoints_per_row-1] + data[flag][datapoints_per_row-2][datapoints_per_row-1] + data[flag][datapoints_per_row-1][datapoints_per_row-2])/3;
    }
    else{
        data[flag_new][datapoints_per_row-1][datapoints_per_row-1] = (data[flag][datapoints_per_row-2][datapoints_per_row-1] + data[flag][datapoints_per_row-1][datapoints_per_row-2])/2;
    }

}

int main( int argc, char *argv[])
{
    //seed
    srand(time(NULL));
    int myrank;
    
    MPI_Init(&argc,&argv);

    // input
    int processes_per_row = atoi(argv[1]);
    int datapoints_per_row = atoi(argv[2]);
    int given_time = atoi(argv[3]);

    // MPI variables
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Request sendrequest[4*datapoints_per_row];
    MPI_Request recieverequest[4*datapoints_per_row];
    MPI_Status status[4*datapoints_per_row];
    MPI_Status status2[4*datapoints_per_row];

    // buffer arrays to be used in 3 methods
    double sendbuf[4][datapoints_per_row];
    double ** receivebuf = (double **) malloc(sizeof(double*) * 4);
    for (int i = 0; i < 4; i++){
        receivebuf[i] = (double *) malloc(sizeof(double) * datapoints_per_row);
    }

    //reporting times 
    double stime, etime, time1, time2, time3;

    // finding out left, right , top , bottom neighbours. if not there, then the placeholder contains -1.
    // ranks for processes
    int left,right,top,bot;
    left = -1;
    right = -1;
    top = -1;
    bot = -1;
    if (myrank < processes_per_row){
        if(myrank == 0){
            right = myrank + 1;
            bot = myrank + processes_per_row;
        }
        else if(myrank == processes_per_row-1){
            left = myrank - 1;
            bot = myrank + processes_per_row;
        }
        else{
            left = myrank - 1;
            right = myrank + 1;
            bot = myrank + processes_per_row;
        }
    }
    else if(myrank >= processes_per_row*(processes_per_row-1)){
        if(myrank == processes_per_row*(processes_per_row-1)){
            right = myrank + 1;
            top = myrank - processes_per_row;
        }
        else if(myrank == processes_per_row*processes_per_row-1){
            left = myrank - 1;
            top = myrank - processes_per_row;
        }
        else{
            left = myrank - 1;
            right = myrank + 1;
            top = myrank - processes_per_row;
        }
    }
    else if(myrank % processes_per_row == 0){
        if(myrank != 0 && myrank != processes_per_row*(processes_per_row-1)){
            top = myrank - processes_per_row;
            bot = myrank + processes_per_row;
            right = myrank + 1;
        }
    }
    else if((myrank-processes_per_row + 1) % processes_per_row == 0){
        if(myrank != processes_per_row-1 && myrank != processes_per_row*processes_per_row-1){
            top = myrank - processes_per_row;
            bot = myrank + processes_per_row;
            left = myrank - 1;
        }
    }
    else{
        top = myrank - processes_per_row;
        bot = myrank + processes_per_row;
        left = myrank - 1;
        right = myrank + 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    

    /////////////////////////   Data initialisation    /////////////////////////////////////////////////////
    
    // make sure that memory allocated in malloc is contiguous for MPI_PACK buffer
    void * addr = malloc(sizeof(double) * datapoints_per_row * datapoints_per_row * 2);
    int start;
    double *** data = (double ***) malloc(sizeof(double**) * 2);
    for (int i = 0; i < 2; i++){
        data[i] = (double **) malloc(sizeof(double*) * datapoints_per_row);
    }
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < datapoints_per_row; j++){
            data[i][j] = (double *) (addr + ((i*datapoints_per_row*datapoints_per_row) + j*datapoints_per_row)*sizeof(double));
        }
    } 
    int flag = 0;
    for(int i = 0; i < datapoints_per_row; i++){
        for(int j = 0; j < datapoints_per_row; j++){
            // random initialisation
            data[0][i][j] = rand();
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////// method 1 : send each double seperately /////////////////////////////////////////////////////
    start = 0;
    stime = MPI_Wtime();
    // runs for given_time time steps
    for(int time = 0; time < given_time; time++){
        start = 0;
        if(left != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                MPI_Issend(&data[flag][j][0],1,MPI_DOUBLE,left,LEFT*datapoints_per_row+j,MPI_COMM_WORLD,&sendrequest[start+j]);
                MPI_Irecv(&receivebuf[LEFT][j],1,MPI_DOUBLE,left,RIGHT*datapoints_per_row+j,MPI_COMM_WORLD,&recieverequest[start+j]); 
            }
            start += datapoints_per_row;
        }
        if(right != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                MPI_Issend(&data[flag][j][datapoints_per_row-1],1,MPI_DOUBLE,right,RIGHT*datapoints_per_row+j,MPI_COMM_WORLD,&sendrequest[start+j]);
                MPI_Irecv(&receivebuf[RIGHT][j],1,MPI_DOUBLE,right,LEFT*datapoints_per_row+j,MPI_COMM_WORLD,&recieverequest[start+j]);
            }
            start += datapoints_per_row;
        }
        if(top != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                MPI_Issend(&data[flag][0][j],1,MPI_DOUBLE,top,TOP*datapoints_per_row+j,MPI_COMM_WORLD,&sendrequest[start+j]);
                MPI_Irecv(&receivebuf[TOP][j],1,MPI_DOUBLE,top,BOT*datapoints_per_row+j,MPI_COMM_WORLD,&recieverequest[start+j]);      
            }
            start += datapoints_per_row;
        }
        if(bot != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                MPI_Issend(&data[flag][datapoints_per_row-1][j],1,MPI_DOUBLE,bot,BOT*datapoints_per_row+j,MPI_COMM_WORLD,&sendrequest[start+j]);
                MPI_Irecv(&receivebuf[BOT][j],1,MPI_DOUBLE,bot,TOP*datapoints_per_row+j,MPI_COMM_WORLD,&recieverequest[start+j]);
            }
            start += datapoints_per_row;
        }
        
        MPI_Waitall(start, recieverequest,status);
        run_stencil_computation(flag,datapoints_per_row,receivebuf,left,right,top,bot,data);
        MPI_Waitall(start,sendrequest,status2);
        flag = flag ^ 1;
    }
    etime = MPI_Wtime();
    time1 = etime - stime;
 
    ///////////////////////////////////////////////////////////////////////////////////////
    

    /////////////////////////   Data initialisation for method 2   /////////////////////////////////////////////////////
    
    int position;
    int type_size;
    MPI_Type_size(MPI_DOUBLE, &type_size);
    flag = 0;
    for(int i = 0; i < datapoints_per_row; i++){
        for(int j = 0; j < datapoints_per_row; j++){
            // reinitialize
            data[0][i][j] = rand();
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////// method 2 : packing data in user buffer /////////////////////////////////////////////////////
   
    start = 0;
    stime = MPI_Wtime();
    for(int time = 0; time < given_time; time++){
        start = 0;
        if(left != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                // pack call
		        MPI_Pack (&data[flag][j][0], 1, MPI_DOUBLE, &sendbuf[LEFT][0], datapoints_per_row*type_size, &position, MPI_COMM_WORLD);
            }
            MPI_Issend(&sendbuf[LEFT][0],position,MPI_PACKED,left,LEFT,MPI_COMM_WORLD,&sendrequest[start]);
            position = 0;
            MPI_Irecv(&receivebuf[LEFT][0],datapoints_per_row,MPI_DOUBLE,left,RIGHT,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }
        if(right != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                // pack call
                MPI_Pack (&data[flag][j][datapoints_per_row-1], 1, MPI_DOUBLE, &sendbuf[RIGHT][0], datapoints_per_row*type_size, &position, MPI_COMM_WORLD);
            }
            MPI_Issend(&sendbuf[RIGHT][0],position,MPI_PACKED,right,RIGHT,MPI_COMM_WORLD,&sendrequest[start]);
            position = 0;
            MPI_Irecv(&receivebuf[RIGHT][0],datapoints_per_row,MPI_DOUBLE,right,LEFT,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }
        if(top != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                // pack call
                MPI_Pack (&data[flag][0][j], 1, MPI_DOUBLE, &sendbuf[TOP][0], datapoints_per_row*type_size, &position, MPI_COMM_WORLD);
            }
            MPI_Issend(&sendbuf[TOP][0],position,MPI_PACKED,top,TOP,MPI_COMM_WORLD,&sendrequest[start]);
            position = 0;
            MPI_Irecv(&receivebuf[TOP][0],datapoints_per_row,MPI_DOUBLE,top,BOT,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }
        if(bot != -1){
            for(int j = 0; j < datapoints_per_row; j++){
                // pack call
                MPI_Pack (&data[flag][datapoints_per_row-1][j], 1, MPI_DOUBLE, &sendbuf[BOT][0], datapoints_per_row*type_size, &position, MPI_COMM_WORLD);
            }
            MPI_Issend(&sendbuf[BOT][0],position,MPI_PACKED,bot,BOT,MPI_COMM_WORLD,&sendrequest[start]);
            position = 0;
            MPI_Irecv(&receivebuf[BOT][0],datapoints_per_row,MPI_DOUBLE,bot,TOP,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }

        MPI_Waitall(start,recieverequest,status);
        run_stencil_computation(flag,datapoints_per_row,receivebuf,left,right,top,bot,data);
        MPI_Waitall(start,sendrequest,status2);
        flag = flag ^ 1;
    }
    etime = MPI_Wtime();
    time2 = etime - stime;




    ///////////////////////////////////////////////////////////////////////////////////////
    

    /////////////////////////   Data initialisation for method 3   /////////////////////////////////////////////////////
    
    MPI_Datatype type1, type2;
    int count;
    int blocklen;
    int stride;
    
    // creation of column vector datatype
    count = datapoints_per_row;
    blocklen=1;
    stride=datapoints_per_row;
    MPI_Type_vector( count , blocklen , stride , MPI_DOUBLE , &type1);
    MPI_Type_commit(&type1);
    
    
    // creation of row vector datatype
    count = datapoints_per_row;
    blocklen=1;
    stride=1;
    MPI_Type_vector( count , blocklen , stride , MPI_DOUBLE , &type2);
    MPI_Type_commit(&type2);

    flag = 0;
    for(int i = 0; i < datapoints_per_row; i++){
        for(int j = 0; j < datapoints_per_row; j++){
            // reinitialize
            data[0][i][j] = rand();
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////// method 3 : using derived data types /////////////////////////////////////////////////////
    
    //send
    start = 0;
    stime = MPI_Wtime();
    for(int time = 0; time < given_time; time++){
        start = 0;
        if(left != -1){
            // sending left column
            MPI_Isend(&data[flag][0][0],1,type1,left,LEFT,MPI_COMM_WORLD,&sendrequest[start]);
            MPI_Irecv(&receivebuf[LEFT][0],datapoints_per_row,MPI_DOUBLE,left,RIGHT,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }
        if(right != -1){
            // sending right column
            MPI_Isend(&data[flag][0][datapoints_per_row-1],1,type1,right,RIGHT,MPI_COMM_WORLD,&sendrequest[start]);
            MPI_Irecv(&receivebuf[RIGHT][0],datapoints_per_row,MPI_DOUBLE,right,LEFT,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }
        if(top != -1){
            // sending top row
            MPI_Isend(&data[flag][0][0],1,type2,top,TOP,MPI_COMM_WORLD,&sendrequest[start]);
            MPI_Irecv(&receivebuf[TOP][0],datapoints_per_row,MPI_DOUBLE,top,BOT,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }
        if(bot != -1){
            // sending bottom row
            MPI_Isend(&data[flag][datapoints_per_row-1][0],1,type2,bot,BOT,MPI_COMM_WORLD,&sendrequest[start]);
            MPI_Irecv(&receivebuf[BOT][0],datapoints_per_row,MPI_DOUBLE,bot,TOP,MPI_COMM_WORLD,&recieverequest[start]);
            start+=1;
        }

        MPI_Waitall(start,recieverequest,status);
        run_stencil_computation(flag,datapoints_per_row,receivebuf,left,right,top,bot,data);
        MPI_Waitall(start,sendrequest,status2);
        flag = flag ^ 1;
    }
    etime = MPI_Wtime();
    time3 = etime - stime;
    
    // report final times
    if(myrank == 0){
        printf("%lf,%lf,%lf\n",time1,time2,time3);
    }
    MPI_Finalize();
}
