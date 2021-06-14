#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <math.h>


struct host_info{
    int host_no; // initial host assigned by mpich
    int group_no; // host group number
    int global_id; //modified parameter to decide new ranks in type 1 modification
    int old_rank; // rank in MPI_Comm_world
    int new_global_rank; 
    int new_intra_comm_rank;
    int new_inter_comm_rank;
    int is_root; // 0 meaning no, 1 meaning yes
    int total_groups;
    int root_group_no;
} typedef host_info;

int find_my_topology(host_info* hinfo, int roothost_no){

    if(hinfo->host_no < 13){
        hinfo->group_no = 0;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no == 13){
        hinfo->group_no = 1;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no < 17){
        hinfo->group_no = 0;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no < 31){
        hinfo->group_no = 1;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no == 31){
        hinfo->group_no = 0;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no == 32){
        hinfo->group_no = 1;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no < 45){
        hinfo->group_no = 2;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no == 45){
        hinfo->group_no = 3;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no == 46){
        hinfo->group_no = 2;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no < 62){
        hinfo->group_no = 3;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no < 79){
        hinfo->group_no = 4;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else if(hinfo->host_no < 93){
        hinfo->group_no = 5;
        hinfo->global_id = 100*(hinfo->group_no) + hinfo->host_no;
    }
    else{
        return -1;
    }

    if(roothost_no < 0) return 0;
    // root group info
    host_info*rhinfo = (host_info*)malloc(sizeof(host_info)); 
    rhinfo->host_no = roothost_no;
    find_my_topology(rhinfo,-1);
    hinfo->root_group_no=rhinfo->group_no;
    if(hinfo->group_no == rhinfo->group_no){
        hinfo->global_id -= 100*(hinfo->group_no);
    }
    else{
        hinfo->global_id += 100;
    }

    return 0;
}

void new_bcast(void* send_data,int send_count,MPI_Datatype send_datatype,int root,MPI_Comm communicator){
    // initialisation
    int len,roothost_no;
    host_info* hinfo = (host_info*)malloc(sizeof(host_info));
    hinfo->is_root=-1;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(communicator, &(hinfo->old_rank));
    MPI_Get_processor_name (hostname, &len);
    hinfo->host_no = atoi(&hostname[5]);
    if(hinfo->old_rank == root){
        roothost_no = hinfo->host_no;
    }
    MPI_Bcast(&roothost_no,1,MPI_INT, root, communicator);
    
    if(find_my_topology(hinfo,roothost_no) < 0)
        if(hinfo->old_rank == 0)
            printf("failed");
    if(hinfo->old_rank == root){
        hinfo->global_id = 0;
    }

    // new comm
    MPI_Comm globalcomm;
    int color = 0;    
    MPI_Comm_split (communicator, color, hinfo->global_id, &globalcomm);
    MPI_Comm_rank (globalcomm, &(hinfo->new_global_rank));  
    MPI_Bcast(send_data,send_count,send_datatype,0,globalcomm);
    MPI_Comm_free( &globalcomm);
}
void new_bcast_alpha(void* send_data,int send_count,MPI_Datatype send_datatype,int root,MPI_Comm communicator){
    // initialisation
    int len,roothost_no;
    host_info* hinfo = (host_info*)malloc(sizeof(host_info));
    hinfo->is_root=-1;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(communicator, &(hinfo->old_rank));
    MPI_Get_processor_name (hostname, &len);
    hinfo->host_no = atoi(&hostname[5]);
    if(hinfo->old_rank == root){
        roothost_no = hinfo->host_no;
    }
    MPI_Bcast(&roothost_no,1,MPI_INT, root, communicator);
    
    if(find_my_topology(hinfo,roothost_no) < 0)
        if(hinfo->old_rank == 0)
            printf("failed");
    if(hinfo->old_rank == root){
        hinfo->global_id = 0;
    }

    // creating intra group comm
    int color;
    color = hinfo->group_no;
    MPI_Comm intercomm,intracomm;
    MPI_Comm_split (communicator, color, hinfo->global_id, &intracomm);
    MPI_Comm_rank (intracomm, &(hinfo->new_intra_comm_rank)); 

    if(hinfo->new_intra_comm_rank == 0){
        color = 0;
    }
    else{
        color = 1; 
    }
    
    
    MPI_Comm_split (communicator, color, hinfo->global_id, &intercomm);
    MPI_Comm_rank (intercomm, &(hinfo->new_inter_comm_rank)); 
    if(hinfo->new_intra_comm_rank == 0) MPI_Bcast(send_data,send_count,send_datatype,0,intercomm);
    MPI_Bcast(send_data,send_count,send_datatype,0,intracomm);
    MPI_Comm_free( &intracomm);
    MPI_Comm_free( &intercomm);
}
void new_bcast_beta(void* send_data,int send_count,MPI_Datatype send_datatype,int root,MPI_Comm communicator){
    // initialisation
    int len,roothost_no;
    host_info* hinfo = (host_info*)malloc(sizeof(host_info));
    hinfo->is_root=-1;hinfo->total_groups=0;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(communicator, &(hinfo->old_rank));
    MPI_Get_processor_name (hostname, &len);
    hinfo->host_no = atoi(&hostname[5]);
    if(hinfo->old_rank == root){
        roothost_no = hinfo->host_no;
    }
    MPI_Bcast(&roothost_no,1,MPI_INT, root, communicator);
    if(find_my_topology(hinfo,roothost_no) < 0)
        if(hinfo->old_rank == 0)
            printf("failed");
    if(hinfo->old_rank == root){
        hinfo->global_id = 0;
    }

    // new comm
    MPI_Comm intracomm;
    MPI_Comm_split (communicator, hinfo->group_no, hinfo->global_id, &intracomm);
    MPI_Comm_rank(intracomm,&(hinfo->new_intra_comm_rank));
    int size_global;
    MPI_Comm_size(communicator,&size_global);
    
    int send_meta_buf[2];
    send_meta_buf[0]=-1;
    send_meta_buf[1]=-1;
    int *recv_meta_buf=(int*)malloc(sizeof(int)*2*size_global);
    int group_leader_old_rank[6]; // -1 meaning that group no doesn't exist for this case
    for(int i=0;i<6;i++) group_leader_old_rank[i]=-1;
    if(hinfo->new_intra_comm_rank == 0){
        send_meta_buf[0]=hinfo->group_no;
        send_meta_buf[1]=hinfo->old_rank;
    }
    MPI_Allgather( send_meta_buf , 2 , MPI_INT , recv_meta_buf , 2 , MPI_INT , communicator);
    for(int i=0;i<2*size_global;i+=2){
        if(recv_meta_buf[i]!=-1){
            group_leader_old_rank[recv_meta_buf[i]]=recv_meta_buf[i+1];
            hinfo->total_groups+=1;
        }
    }

    MPI_Comm intercomms[5];
    if(hinfo->group_no == hinfo->root_group_no){
        for(int i=0;i<6;i++){
            if(group_leader_old_rank[i]!=-1 && i!=hinfo->group_no){
                MPI_Intercomm_create( intracomm , 0 , communicator , group_leader_old_rank[i] , i+2 , &intercomms[i]);        
            }
        }
    }
    else MPI_Intercomm_create( intracomm , 0 , communicator , group_leader_old_rank[hinfo->root_group_no] , hinfo->group_no+2 , &intercomms[hinfo->group_no]);
    if(hinfo->group_no==hinfo->root_group_no){
        MPI_Bcast(send_data,send_count,send_datatype,0,intracomm);
        if(hinfo->old_rank==root){
            for(int i=0;i<6;i++){
                if(group_leader_old_rank[i]!=-1 && i!=hinfo->root_group_no) MPI_Bcast(send_data,send_count,send_datatype,MPI_ROOT,intercomms[i]); 
            }
        }
        else{
            for(int i=0;i<6 ;i++){
                if(group_leader_old_rank[i]!=-1 && i!=hinfo->root_group_no) MPI_Bcast(send_data,send_count,send_datatype,MPI_PROC_NULL,intercomms[i]); 
            }
        }
            
    }
    else{
        MPI_Bcast(send_data,send_count,send_datatype,0,intercomms[hinfo->group_no]); 
    }
    MPI_Comm_free( &intracomm);
}


void new_gather( void* send_data,int send_count,MPI_Datatype send_datatype,void* recv_data,int recv_count,MPI_Datatype recv_datatype,int root,MPI_Comm communicator) {
    int len,roothost_no;
    host_info* hinfo = (host_info*)malloc(sizeof(host_info));
    hinfo->is_root=-1;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(communicator, &(hinfo->old_rank));
    MPI_Get_processor_name (hostname, &len);
    hinfo->host_no = atoi(&hostname[5]);
    if(hinfo->old_rank == root){
        roothost_no = hinfo->host_no;
    }        
    MPI_Bcast(&roothost_no,1,MPI_INT, root, communicator);

    if(find_my_topology(hinfo,roothost_no) < 0)
        if(hinfo->old_rank == 0)
            printf("failed");
    if(hinfo->old_rank == root){
        hinfo->global_id = 0;
    }

    // creating intra group comm
    int color;
    color = hinfo->group_no;
    MPI_Comm intercomm,intracomm;
    MPI_Comm_split (communicator, color, hinfo->global_id, &intracomm);
    MPI_Comm_rank (intracomm, &(hinfo->new_intra_comm_rank)); 

    if(hinfo->new_intra_comm_rank == 0){
        color = 0;
    }
    else{
        color = 1; 
    }
    
    
    MPI_Comm_split (communicator, color, hinfo->global_id, &intercomm);
    MPI_Comm_rank (intercomm, &(hinfo->new_inter_comm_rank)); 
           
    // int newrank;
    MPI_Comm globalcomm;
    color = 0;    
    MPI_Comm_split (communicator, color, hinfo->global_id, &globalcomm);
    MPI_Comm_rank (globalcomm, &(hinfo->new_global_rank));  

    int intracommsize;
    MPI_Comm_size(intracomm, &intracommsize);

    int intercommsize;
    MPI_Comm_size(intercomm, &intercommsize);

    int globalsize;
    MPI_Comm_size(globalcomm, &globalsize);

    double* recv_intra_data=(double*)malloc(sizeof(double)*recv_count*intracommsize);
    double* recv_inter_data=(double*)malloc(sizeof(double)*recv_count*globalsize);    
    // intra communication
    MPI_Gather(send_data,send_count,send_datatype,recv_intra_data,recv_count,recv_datatype, 0 ,intracomm);
    
    // metadata for intercomm
    int *countarray=(int*)malloc(sizeof(int)*intercommsize);
    int *displarray=(int*)malloc(sizeof(int)*intercommsize);

    int sendval = recv_count*intracommsize;

    if(hinfo->new_intra_comm_rank==0) MPI_Gather(&sendval,1,MPI_INT, countarray,1,MPI_INT,0,intercomm);

    int sum = 0;
    for(int x = 0; x < intercommsize; x++){
        displarray[x] = sum;
        sum += countarray[x];
    }
    // inter commincation
    if(hinfo->new_intra_comm_rank==0) MPI_Gatherv (recv_intra_data, recv_count*intracommsize, send_datatype, recv_inter_data, countarray, displarray, recv_datatype, 0, intercomm);
    
    int rank_pair[2]; // 0 position for old rank, 1 for new global rank
    rank_pair[0]=hinfo->old_rank;
    rank_pair[1]=hinfo->new_global_rank;

    int* rank_pairs=(int*)malloc(sizeof(int)*2*globalsize);
    // for jumbling the array
    MPI_Gather(rank_pair,2,MPI_INT,rank_pairs,2,MPI_INT,0,globalcomm);
    
    // final jumbling
    if(hinfo->new_global_rank==0){
        for(int i=0;i<2*globalsize;i+=2){
            memcpy(recv_data+rank_pairs[i]*recv_count*sizeof(double),recv_inter_data+rank_pairs[i+1]*recv_count,sizeof(double)*recv_count);
        }
    }
    MPI_Comm_free( &globalcomm);
    MPI_Comm_free( &intercomm);
    MPI_Comm_free( &intracomm);
    

}

void new_reduce(void* send_buf,void* recv_buf,int count,MPI_Datatype datatype,MPI_Op op,int root,MPI_Comm communicator){
    // special commutative ops - MPI_MAXLOC, MPI_MINLOC, returns rank number also
    int len,roothost_no;
    host_info* hinfo = (host_info*)malloc(sizeof(host_info));
    hinfo->is_root=-1;
    hinfo->new_inter_comm_rank=-1;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(communicator, &(hinfo->old_rank));
    MPI_Get_processor_name (hostname, &len);
    hinfo->host_no = atoi(&hostname[5]);
    if(hinfo->old_rank == root){
        roothost_no = hinfo->host_no;
    }        
    MPI_Bcast(&roothost_no,1,MPI_INT, root, communicator);

    if(find_my_topology(hinfo,roothost_no) < 0)
        if(hinfo->old_rank == 0)
            printf("failed");
    if(hinfo->old_rank == root){
        hinfo->global_id = 0;
    }

    // creating intra group comm
    int color;
    color = hinfo->group_no;
    MPI_Comm intercomm,intracomm;
    MPI_Comm_split (communicator, color, hinfo->global_id, &intracomm);
    MPI_Comm_rank (intracomm, &(hinfo->new_intra_comm_rank)); 

    if(hinfo->new_intra_comm_rank == 0){
        color = 0;
    }
    else{
        color = 1; 
    }
    
    
    MPI_Comm_split (communicator, color, hinfo->global_id, &intercomm);
    MPI_Comm_rank (intercomm, &(hinfo->new_inter_comm_rank)); 
    
    double *recv_intra_data=(double*)malloc(sizeof(double)*count);
    double *recv_inter_data=(double*)malloc(sizeof(double)*count);    
    // intra communication
    MPI_Reduce(send_buf,recv_intra_data,count,datatype, op,0 ,intracomm);
    
    if(hinfo->new_intra_comm_rank==0) MPI_Reduce(recv_intra_data,recv_buf,count,datatype, op,0 ,intercomm);
    
    MPI_Comm_free( &intercomm);
    MPI_Comm_free( &intracomm);
}
void new_alltoallv(const void *sendbuf, const int *sendcounts,
                  const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
                  const int *recvcounts, const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm){
    int len,roothost_no;
    host_info* hinfo = (host_info*)malloc(sizeof(host_info));
    hinfo->is_root=-1;
    hinfo->new_inter_comm_rank=-1;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(comm, &(hinfo->old_rank));
    MPI_Get_processor_name (hostname, &len);
    hinfo->host_no = atoi(&hostname[5]);
    
    if(find_my_topology(hinfo,-1) < 0)
        if(hinfo->old_rank == 0)
            printf("failed");
    
    // creating intra group comm
    int color;
    color = hinfo->group_no;
    MPI_Comm intercomm,intracomm;


    MPI_Comm_split (comm, color, hinfo->global_id, &intracomm);
    MPI_Comm_rank (intracomm, &(hinfo->new_intra_comm_rank)); 

    if(hinfo->new_intra_comm_rank == 0){
        color = 0;
    }
    else{
        color = 1; 
    }
    
    
    MPI_Comm_split (comm, color, hinfo->global_id, &intercomm);
    MPI_Comm_rank (intercomm, &(hinfo->new_inter_comm_rank)); 
           
    // int newrank;
    MPI_Comm globalcomm;
    color = 0;    
    MPI_Comm_split (comm, color, hinfo->global_id, &globalcomm);
    MPI_Comm_rank (globalcomm, &(hinfo->new_global_rank));  
    
    int intracommsize;
    MPI_Comm_size(intracomm, &intracommsize);

    int intercommsize;
    MPI_Comm_size(intercomm, &intercommsize);

    int globalsize;
    MPI_Comm_size(globalcomm, &globalsize);
    
    int r = hinfo->new_inter_comm_rank;
    MPI_Bcast(&r,1,MPI_INT,0,intracomm);

    int rank_pairs2[3*globalsize];
    int rank_pairs[3*globalsize];
    
    rank_pairs2[0]=hinfo->old_rank;
    rank_pairs2[1]=hinfo->new_global_rank;
    rank_pairs2[2]=r;
    

    MPI_Allgather(rank_pairs2,3,MPI_INT,rank_pairs,3,MPI_INT,globalcomm);
    
    
    int old_to_new[globalsize];
    for(int a =0; a < globalsize; a++){
        old_to_new[rank_pairs[3*a]] = a;
    }

    int* sendcountsarray = (int *) malloc(sizeof(int)*intracommsize*globalsize);
    int* sendcountsarraydisp = (int *) malloc(sizeof(int)*intracommsize*globalsize);
    int* recvcountsarray = (int *) malloc(sizeof(int)*intracommsize*globalsize);
    int* recvcountsarraydisp = (int *) malloc(sizeof(int)*intracommsize*intercommsize);
    int* recvcountsarraynew = (int *) malloc(sizeof(int)*intracommsize*intercommsize);
    
    
    
    MPI_Gather(sendcounts,globalsize,MPI_INT,sendcountsarray,globalsize,MPI_INT,0,intracomm);
    MPI_Gather(recvcounts,globalsize,MPI_INT,recvcountsarray,globalsize,MPI_INT,0,intracomm);
    
    int* leaderrecvcounts = (int *) malloc(sizeof(int)*intracommsize);
    int* leaderrdisps = (int *) malloc(sizeof(int)*intracommsize);

    int zs = 0;
    int temp1 = 0;
    if(hinfo->new_intra_comm_rank == 0) {
        for(int a = 0;a < intracommsize; a++){
            leaderrdisps[a] = zs;
            temp1 = 0;
            for(int b = 0; b < globalsize; b++){
                temp1 += sendcountsarray[a*globalsize+b];
            }
            leaderrecvcounts[a] = temp1;
            zs += temp1;
        }
    }

    

    int intrasize = 0;
    for(int x = 0; x < globalsize; x++){
        intrasize += sendcounts[x];
    }

    int totalsize = 0;
    double* intrabuf;
    if(hinfo->new_intra_comm_rank == 0){
        for(int x = 0; x < globalsize*intracommsize; x++){
            sendcountsarraydisp[x] = totalsize;
            totalsize += sendcountsarray[x];
        }
        intrabuf = (double *) malloc(sizeof(double)*totalsize);
    }
    else{
        intrabuf = (double *) malloc(sizeof(double)*intrasize);
    }

    double* temp = intrabuf;
    double* temp3 = (double *)sendbuf;
    for(int x=0; x < globalsize; x++){
        memcpy(temp,temp3 + sdispls[x],sizeof(double)*sendcounts[x]);
        temp += sendcounts[x];
    }
    
    if(hinfo->new_intra_comm_rank == 0) {
        MPI_Gatherv(MPI_IN_PLACE,0, MPI_DATATYPE_NULL,intrabuf,leaderrecvcounts,leaderrdisps,MPI_DOUBLE,0,intracomm);
    }
    else{
        MPI_Gatherv(intrabuf,intrasize, MPI_DOUBLE,intrabuf,leaderrecvcounts,leaderrdisps,MPI_DOUBLE,0,intracomm); 
    }
    
    
    double* interbuf = (double *) malloc(sizeof(double)*totalsize);
    
    temp = interbuf;

    int* intersendcounts = (int *) malloc(sizeof(int)*intercommsize);
    int* intersdisps = (int *) malloc(sizeof(int)*intercommsize);
    int* interrecvcounts = (int *) malloc(sizeof(int)*intercommsize);
    int* interrdisps = (int *) malloc(sizeof(int)*intercommsize);

    if(hinfo->new_intra_comm_rank == 0){
        intersdisps[0] = 0;    
        for(int x = 0; x < intercommsize; x++){
            intersendcounts[x] = 0;
        }
        int newa = 0;
        for (int a = 0; a < globalsize; a++){
            newa = rank_pairs[3*a];
            for(int b = 0; b < intracommsize; b++){
                memcpy(temp,intrabuf+sendcountsarraydisp[newa+b*globalsize],sizeof(double)*sendcountsarray[newa+b*globalsize]);
                temp += sendcountsarray[newa+b*globalsize];
                //final
                intersendcounts[rank_pairs[3*a+2]] += sendcountsarray[newa+b*globalsize];
            }     
        }
    }

    
    int totalrecvsize = 0;
    if(hinfo->new_intra_comm_rank == 0){


        for(int x = 1; x < intercommsize; x++){
            intersdisps[x] = intersdisps[x-1] + intersendcounts[x-1]  ;
        }


        for(int a = 0; a < intercommsize; a++){
            interrecvcounts[a] = 0;
        }

        for (int a = 0; a < globalsize; a++){
            int newa = rank_pairs[3*a];
            for(int b = 0; b < intracommsize; b++){
                interrecvcounts[rank_pairs[3*a+2]] += recvcountsarray[newa+b*globalsize];
                totalrecvsize += recvcountsarray[newa+b*globalsize];
            }     
        }

        interrdisps[0] = 0;
        for(int x = 1; x < intercommsize; x++){
            interrdisps[x] = interrdisps[x-1] + interrecvcounts[x-1]  ;
        }
    
    }
    double * interrecvbuf = (double*) malloc(sizeof(double)*totalrecvsize);
    if(hinfo->new_intra_comm_rank == 0) MPI_Alltoallv(interbuf,intersendcounts,intersdisps,MPI_DOUBLE,interrecvbuf,interrecvcounts,interrdisps,MPI_DOUBLE,intercomm);
    
    // scatter
//** totalrecvsize not defined for non leader processes

    intrasize = 0;
        for(int x = 0; x < globalsize; x++){
            intrasize += recvcounts[x];
        }

    double * intrarecvbuf;
    if(hinfo->new_intra_comm_rank == 0){
        intrarecvbuf = (double*) malloc(sizeof(double)*totalrecvsize);
    }
    else{
        intrarecvbuf = (double*) malloc(sizeof(double)*intrasize);
    }
     
    temp = intrarecvbuf;
    int ts = 0;
    int y = 0;

    

    if(hinfo->new_intra_comm_rank == 0) {
        for (int a = 0; a < intercommsize*intracommsize; a++){
            recvcountsarraydisp[a] = 0;
            recvcountsarraynew[a] = 0;
        }
        for (int a = 0; a < globalsize; a++){
            for(int b = 0; b < intracommsize; b++){
                int c = rank_pairs[3*a+2];
                // final
                recvcountsarraydisp[b+c*intracommsize] += recvcountsarray[rank_pairs[3*a]+b*globalsize];
                recvcountsarraynew[b+c*intracommsize] += recvcountsarray[rank_pairs[3*a]+b*globalsize];
            }     
        }
        int tp = 0;
        int tp2;
        for (int a = 0; a < intercommsize*intracommsize; a++){
            // if(a%intracommsize == 0) tp = 0;
            tp2 = recvcountsarraydisp[a];
            recvcountsarraydisp[a] = tp;
            tp += tp2;
        }
    }

    int* intrarecvcounts = (int *) malloc(sizeof(int)*intracommsize);
    int* intrardisps = (int *) malloc(sizeof(int)*intracommsize);
    int ys = 0;
    int totals = 0;
    
    if(hinfo->new_intra_comm_rank == 0) {
        for (int a = 0; a < intracommsize; a++){
            intrarecvcounts[a] = 0;
        }
        for (int a = 0; a < intracommsize; a++){
            
            intrardisps[a] = totals;
            ys = 0;
            for(int b = 0; b < intercommsize; b++){
                memcpy(temp,interrecvbuf+recvcountsarraydisp[a+b*intracommsize],sizeof(double)*recvcountsarraynew[a+b*intracommsize]);
                temp += recvcountsarraynew[a+b*intracommsize];
                ys += recvcountsarraynew[a+b*intracommsize];
            }  
            totals += ys;
            intrarecvcounts[a] = ys;    
        }
    }


    if(hinfo->new_intra_comm_rank == 0){
        MPI_Scatterv(intrarecvbuf,intrarecvcounts,intrardisps,MPI_DOUBLE,MPI_IN_PLACE,intrasize, MPI_DOUBLE,0,intracomm);
    }
    else{
        MPI_Scatterv(intrarecvbuf,intrarecvcounts,intrardisps,MPI_DOUBLE,intrarecvbuf,intrasize, MPI_DOUBLE,0,intracomm);
    }

    double* temp2 = (double*) recvbuf;
    temp = intrarecvbuf;  
    for(int b = 0; b < globalsize; b++){
            memcpy(temp2+rdispls[old_to_new[b]],temp,sizeof(double)*recvcounts[old_to_new[b]]);
            temp += recvcounts[old_to_new[b]];
    }
    
}

int main( int argc, char *argv[])
{
    int myrank, len;
    MPI_Init(&argc,&argv);
    int num_process;
    int D = atoi(argv[1])*128; // number of doubles 
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    double otime,mtime,ntime,n0time,n2time;
    double bcast_time,bcast_opt_time,gather_time,gather_opt_time,reduce_time,reduce_opt_time,alltoallv_time,alltoallv_opt_time;

    double* data = (double *) malloc(sizeof(double)*D*num_process); 
    double* data2 =  (double *) malloc(sizeof(double)*D*num_process);
    double* data3 = (double *) malloc(sizeof(double)*D*num_process);
    double* data4 = (double *) malloc(sizeof(double)*D*num_process);
    double* recv = (double *) malloc(sizeof(double)*D*num_process); 
    double* recv2 =  (double *) malloc(sizeof(double)*D*num_process);
    double* recv3 = (double *) malloc(sizeof(double)*D*num_process);

    int root = rand()%num_process;

    srand ( time ( NULL));

    for (int x =0; x < D*num_process; x++){
        data[x] = (double)rand()/RAND_MAX*20000.0-10000.0;
        data2[x] = data[x];
        data3[x] = data[x];
    }

    for (int x =0; x < D; x++){
        data[x] = (double)rand()/RAND_MAX*20000.0-10000.0;
        data2[x] = data[x];
        data3[x] = data[x];
    }

    // bcast
    bcast_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            MPI_Bcast(data, D, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }
    bcast_time = MPI_Wtime() - bcast_time;

    ntime = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            new_bcast(data3, D, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }
    ntime = MPI_Wtime() - ntime;

    bcast_opt_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            new_bcast_alpha(data4, D, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }
    bcast_opt_time = MPI_Wtime() - bcast_opt_time;

    n2time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            new_bcast_beta(data2, D, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }
    n2time = MPI_Wtime() - n2time;

    if(myrank == root)
        printf("%lf,%lf,%lf,%lf\n",bcast_time/5,bcast_opt_time/5,ntime/5,n2time/5);

    
    // gather 
    gather_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            MPI_Gather(data, D, MPI_DOUBLE, recv,D,MPI_DOUBLE,root, MPI_COMM_WORLD);
    }
    gather_time = MPI_Wtime() - gather_time;

    gather_opt_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            new_gather(data, D, MPI_DOUBLE,recv3,D,MPI_DOUBLE, root, MPI_COMM_WORLD);
    }
    gather_opt_time = MPI_Wtime() - gather_opt_time;

    if(myrank == root)
        printf("%lf,%lf\n",gather_time/5,gather_opt_time/5);

    // reduce
    // optimized for all commutative operations, no guarantees of correctness in case of user defined operations
    reduce_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            MPI_Reduce(data, recv,D, MPI_DOUBLE,MPI_SUM,root, MPI_COMM_WORLD);
    }
    reduce_time = MPI_Wtime() - reduce_time;

    reduce_opt_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            new_reduce(data, recv3,D, MPI_DOUBLE,MPI_SUM,root, MPI_COMM_WORLD);
    }
    reduce_opt_time = MPI_Wtime() - reduce_opt_time;

    if(myrank == root)
        printf("%lf,%lf\n",reduce_time/5,reduce_opt_time/5);


    //alltoallv
    int scounts[num_process];
    int sdisps[num_process];
    int rcounts[num_process];
    int rdisps[num_process];

    for(int a = 0; a < num_process; a++){
        scounts[a] = (int)(((double)rand()/RAND_MAX)*D);
    }
    MPI_Alltoall(scounts,1,MPI_INT,rcounts,1,MPI_INT,MPI_COMM_WORLD);
    sdisps[0] = 0;
    rdisps[0] = 0;
    for(int a = 1; a < num_process; a++){
        sdisps[a] = sdisps[a-1] + scounts[a-1];
        rdisps[a] = rdisps[a-1] + rcounts[a-1];
    }
    alltoallv_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            MPI_Alltoallv(data, scounts, sdisps, MPI_DOUBLE, recv,rcounts,rdisps,MPI_DOUBLE, MPI_COMM_WORLD);
    }
    alltoallv_time = MPI_Wtime() - alltoallv_time;

    alltoallv_opt_time = MPI_Wtime();
    for(int step = 0; step < 5; step++ ){
            new_alltoallv(data3, scounts, sdisps, MPI_DOUBLE, recv3,rcounts,rdisps,MPI_DOUBLE, MPI_COMM_WORLD);
    }
    alltoallv_opt_time = MPI_Wtime() - alltoallv_opt_time;

    if(myrank == 0)
        printf("%lf,%lf\n",alltoallv_time/5,alltoallv_opt_time/5);
ret:    
    MPI_Finalize();
}
