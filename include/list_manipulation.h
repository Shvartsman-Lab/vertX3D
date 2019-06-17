#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************LIST MANIPULATION********************************
//****************************************************************************
//****************************************************************************
void remove_edge_from_basal_edges(int edg_id, int cell_id){
    
    for(int i=3; i<=2+basal_edges[cell_id][2]; i++){
        if(abs(basal_edges[cell_id][i])==edg_id){
            for(int j=i; j<2+basal_edges[cell_id][2]; j++){
                basal_edges[cell_id][j]=basal_edges[cell_id][j+1];
            }
        }
    }
    basal_edges[cell_id][2]-=1;
}
//****************************************************************************
int insert_edge_to_basal_edges(int edg_id, int cell_id, int edg_prev, int edg_next){
    
    int vstart;
    
    for(size_t i=3; i<2+basal_edges[cell_id][2]; i++){
        if(abs(basal_edges[cell_id][i])==edg_prev && abs(basal_edges[cell_id][i+1])==edg_next){
            for(size_t j=2+basal_edges[cell_id][2]+1; j>i+2; j--){
                basal_edges[cell_id][j]=basal_edges[cell_id][j-1];
            }
            if(basal_edges[cell_id][i]>0) vstart=e[abs(basal_edges[cell_id][i])][2];
            else if(basal_edges[cell_id][i]<0) vstart=e[abs(basal_edges[cell_id][i])][1];
            if(vstart==e[edg_id][1]){
                basal_edges[cell_id][i+2]=basal_edges[cell_id][i+1];
                basal_edges[cell_id][i+1]=edg_id;
            }
            else if (vstart==e[edg_id][2]){
                basal_edges[cell_id][i+2]=basal_edges[cell_id][i+1];
                basal_edges[cell_id][i+1]=-edg_id;
            }
            basal_edges[cell_id][2]+=1;
            return 0;
        }
        else if(abs(basal_edges[cell_id][i])==edg_next && abs(basal_edges[cell_id][i+1])==edg_prev){
            for(size_t j=2+basal_edges[cell_id][2]+1; j>i+2; j--){
                basal_edges[cell_id][j]=basal_edges[cell_id][j-1];
            }
            if(basal_edges[cell_id][i]>0) vstart=e[abs(basal_edges[cell_id][i])][2];
            else if(basal_edges[cell_id][i]<0) vstart=e[abs(basal_edges[cell_id][i])][1];
            if(vstart==e[edg_id][1]){
                basal_edges[cell_id][i+2]=basal_edges[cell_id][i+1];
                basal_edges[cell_id][i+1]=edg_id;
            }
            else if (vstart==e[edg_id][2]){
                basal_edges[cell_id][i+2]=basal_edges[cell_id][i+1];
                basal_edges[cell_id][i+1]=-edg_id;
            }
            basal_edges[cell_id][2]+=1;
            return 0;
        }
        
    }
    size_t i=2+basal_edges[cell_id][2];
    if(abs(basal_edges[cell_id][i])==edg_prev && abs(basal_edges[cell_id][3])==edg_next){
        if(basal_edges[cell_id][i]>0) vstart=e[abs(basal_edges[cell_id][i])][2];
        else if(basal_edges[cell_id][i]<0) vstart=e[abs(basal_edges[cell_id][i])][1];
        if(vstart==e[edg_id][1]) basal_edges[cell_id][i+1]=edg_id;
        else if (vstart==e[edg_id][2]) basal_edges[cell_id][i+1]=-edg_id;
    }
    else if(abs(basal_edges[cell_id][i])==edg_next && abs(basal_edges[cell_id][3])==edg_prev){
        if(basal_edges[cell_id][i]>0) vstart=e[abs(basal_edges[cell_id][i])][2];
        else if(basal_edges[cell_id][i]<0) vstart=e[abs(basal_edges[cell_id][i])][1];
        if(vstart==e[edg_id][1]) basal_edges[cell_id][i+1]=edg_id;
        else if (vstart==e[edg_id][2]) basal_edges[cell_id][i+1]=-edg_id;
    }
    basal_edges[cell_id][2]+=1;
    return 0;
}
//****************************************************************************
int divide_cell_in_edge_network(int i, int e1, int e2, int ee1, int ee2, int newe){
    
    int basal_edges_copy[15];
    int basal_edges_newc[15];
    
    //MAKES COPY OF CELL i
    int j,k;
    for(j=1; j<=2+basal_edges[i][2]; j++){
        basal_edges_copy[j]=basal_edges[i][j];
    }
    
    //MAKES NEW CELL i
    basal_edges[i][1]=i;
    basal_edges[i][2]=0;
    j=3;
    while(abs(basal_edges_copy[j])!=e1){
        basal_edges[i][j]=basal_edges_copy[j];
        j++;
    }
    //
    k=3;
    if(basal_edges_copy[j]>0){
        basal_edges[i][j]=e1;
        basal_edges_newc[k]=ee1;
    }
    else{
        basal_edges[i][j]=-ee1;
        basal_edges_newc[k]=-e1;
    }
    k++;
    j++;
    //
    basal_edges[i][j]=newe;
    j++;
    //
    int jj=j-1;
    while(abs(basal_edges_copy[jj])!=e2){
        basal_edges_newc[k]=basal_edges_copy[jj];
        k++;
        jj++;
    }
    if(basal_edges_copy[jj]>0){
        basal_edges[i][j]=ee2;
        basal_edges_newc[k]=e2;
    }
    else{
        basal_edges[i][j]=-e2;
        basal_edges_newc[k]=-ee2;
    }
    k++;
    jj++;
    j++;
    //
    while(jj<=2+basal_edges_copy[2]){
        basal_edges[i][j]=basal_edges_copy[jj];
        jj++;
        j++;
    }
    basal_edges[i][2]=j-3;
    //
    basal_edges_newc[k]=-newe;
    k++;
    basal_edges_newc[2]=k-3;
    for(j=k; j<=14; j++){
        basal_edges_newc[j]=0;
    }
    
    //MAKE NEW CELL
    int newc=make_cell(basal_edges_newc[2],basal_edges_newc[3],basal_edges_newc[4],basal_edges_newc[5],basal_edges_newc[6],basal_edges_newc[7],basal_edges_newc[8],basal_edges_newc[9],basal_edges_newc[10],basal_edges_newc[11],basal_edges_newc[12],basal_edges_newc[13],basal_edges_newc[14]);
    
    return newc;
}
//****************************************************************************
void replace_one_edge_with_two_edges(int i, int e1, int ee1){
    
    int j;
    for(j=2+basal_edges[i][2]+1;j>=4;j--){
        if(abs(basal_edges[i][j-1])==e1){
            if(basal_edges[i][j-1]>0) {
                basal_edges[i][j]=ee1;
                basal_edges[i][j-1]=e1;
            }
            else {
                basal_edges[i][j]=-e1;
                basal_edges[i][j-1]=-ee1;
            }
            break;
        }
        else basal_edges[i][j]=basal_edges[i][j-1];
    }
    basal_edges[i][2]+=1;
    
    /*for(j=1; j<=2+basal_edges[i][2]; j++){
     printf("%d ", basal_edges[i][j]);
     }
     printf("\n");*/
    //62	5	936	-937	-711	157	-29
    
    
    
}
//****************************************************************************