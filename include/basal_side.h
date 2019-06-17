#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************BASAL SIDE***************************************
//****************************************************************************
//****************************************************************************
int basal_center(int i, int make_or_not){
    
    int vert_id,vert_ref_id, pass_id;
    double x_cent=0, y_cent=0, z_cent=0;
    double *dxdydz = new double[3];
    dxdydz[0]=0; dxdydz[1]=0; dxdydz[2]=0;
    
    //vert_ref_id
    vert_ref_id=basal_vertices[i][3];
    x_cent+=v[vert_ref_id][1];
    y_cent+=v[vert_ref_id][2];
    z_cent+=v[vert_ref_id][3];
    
    //vert_id
    for(size_t j=4; j <= 2+basal_edges[i][2]; j++){
        vert_id=basal_vertices[i][j];
        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        x_cent += v[vert_id][1] + dxdydz[0];
        y_cent += v[vert_id][2] + dxdydz[1];
        z_cent += v[vert_id][3] + dxdydz[2];
    }

    //divided by number of edges
    x_cent/=(1.*basal_edges[i][2]);
    y_cent/=(1.*basal_edges[i][2]);
    z_cent/=(1.*basal_edges[i][2]);
    
    //MAKE PASSIVE CENTER VERTEX
    if(make_or_not==1){
        pass_id=make_vertex_pass(x_cent,y_cent,z_cent);
        torus_vertex_pass(pass_id);
        v_pass_cell[pass_id]=i;
        c_cent_basal[i]=pass_id;
    }
    //PLACE PASSIVE CENTER VERTEX
    else{
        pass_id=c_cent_basal[i];
        v_pass[pass_id][1]=x_cent;
        v_pass[pass_id][2]=y_cent;
        v_pass[pass_id][3]=z_cent;
        torus_vertex_pass(pass_id);
    }
    
    delete []dxdydz;
    
    return pass_id;
    
}
//****************************************************************************
void set_e_cell(int edge_id, int cell_id){
    
    if      (e_cell1[edge_id]==0 && e_cell2[edge_id]!=cell_id) e_cell1[edge_id]=cell_id;
    else if (e_cell2[edge_id]==0 && e_cell1[edge_id]!=cell_id) e_cell2[edge_id]=cell_id;
    
}
//****************************************************************************
void make_basal_side(int i){
    
    make_basal_vertices(i);
    
    //CALCULATE CELL CENTER
    int pass_id=basal_center(i,1);
    
    //MAKE TRIANGULAR FACETS
    basal_facets[i][1]=basal_edges[i][1];
    basal_facets[i][2]=basal_edges[i][2];
    
    for(size_t j=3; j < 2+basal_edges[i][2]; j++){
        basal_facets[i][j]=make_facet(
                                      basal_vertices[i][j],
                                      basal_vertices[i][j+1],
                                      pass_id
                                      );
        f_T[basal_facets[i][j]]=c_beta[i];//surface tension
        f_cell[basal_facets[i][j]]=i;//f_cell
        set_e_cell(abs(basal_edges[i][j]),i);//e_cell
    }
    int j = 2+basal_vertices[i][2];
    basal_facets[i][j]=make_facet(
                                  basal_vertices[i][j],
                                  basal_vertices[i][3],
                                  pass_id
                                  );
    f_T[basal_facets[i][j]]=c_beta[i];//surface tension
    f_cell[basal_facets[i][j]]=i;//f_cell
    set_e_cell(abs(basal_edges[i][j]),i);//e_cell
    
}
//****************************************************************************
void make_basal_side_ALL(){
    for(int i=1; i<=Nc; i++) make_basal_side(i);
}
//****************************************************************************
void set_basal_tension(int i, double tens){
    c_beta[i]=tens;
    for(size_t j=3; j <= 2+basal_edges[i][2]; j++) f_T[basal_facets[i][j]]=tens;
}
//****************************************************************************










