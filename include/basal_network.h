#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************BASAL NETWORK************************************
//****************************************************************************
//****************************************************************************
int make_cell(int poly_class, int e1, int e2, int e3, int e4, int e5, int e6, int e7, int e8, int e9, int e10, int e11, int e12){
    
    int c_id=c_freeId.get(Nc+1);
    basal_edges[c_id][1]=c_id;
    basal_edges[c_id][2]=poly_class;
    basal_edges[c_id][3]=e1;
    basal_edges[c_id][4]=e2;
    basal_edges[c_id][5]=e3;
    basal_edges[c_id][6]=e4;
    basal_edges[c_id][7]=e5;
    basal_edges[c_id][8]=e6;
    basal_edges[c_id][9]=e7;
    basal_edges[c_id][10]=e8;
    basal_edges[c_id][11]=e9;
    basal_edges[c_id][12]=e10;
    basal_edges[c_id][13]=e11;
    basal_edges[c_id][14]=e12;
    
    //ATTRIBUTES
    c_V0[c_id]=V0;
    c_kV[c_id]=kV;
    if(c_alpha[c_id]<0.0001) c_alpha[c_id]=alph;
    if(c_beta[c_id]<0.0001) c_beta[c_id]=bet;
    
    if(c_id>Nc) Nc++;
    
    return c_id;
}
//****************************************************************************
void set_v_cell(int vert_id, int i){
    
    if      (v_cell1[vert_id]==0 && v_cell2[vert_id]!=i && v_cell3[vert_id]!=i && v_cell4[vert_id]!=i) v_cell1[vert_id]=i;
    else if (v_cell2[vert_id]==0 && v_cell1[vert_id]!=i && v_cell3[vert_id]!=i && v_cell4[vert_id]!=i) v_cell2[vert_id]=i;
    else if (v_cell3[vert_id]==0 && v_cell1[vert_id]!=i && v_cell2[vert_id]!=i && v_cell4[vert_id]!=i) v_cell3[vert_id]=i;
    else if (v_cell4[vert_id]==0 && v_cell1[vert_id]!=i && v_cell2[vert_id]!=i && v_cell3[vert_id]!=i) v_cell4[vert_id]=i;
    
}
//****************************************************************************
int vertex_valence(int i){
    int nrcell=0;
    if(v_cell1[i]!=0) nrcell++;
    if(v_cell2[i]!=0) nrcell++;
    if(v_cell3[i]!=0) nrcell++;
    if(v_cell4[i]!=0) nrcell++;
    
    if      (nrcell==1) return 2;
    else if (nrcell==2) return 3;
    else if (nrcell==3) return 3;
    else if (nrcell==4) return 4;
    else return 0;
}
//****************************************************************************
void make_basal_vertices(int i){
    
    int vert_id;
    
    basal_vertices[i][1]=basal_edges[i][1];
    basal_vertices[i][2]=basal_edges[i][2];
    
    for(size_t j=3; j <= 2+basal_edges[i][2]; j++){
        
        if(basal_edges[i][j]>0) basal_vertices[i][j]=e[basal_edges[i][j]][1];
        if(basal_edges[i][j]<0) basal_vertices[i][j]=e[-basal_edges[i][j]][2];
        
        //v_cell
        set_v_cell(basal_vertices[i][j],i);
    }
}
//****************************************************************************
