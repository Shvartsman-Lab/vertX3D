#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************APICAL SIDE**************************************
//****************************************************************************
//****************************************************************************
int apical_center(int i, int make_or_not){
    
    int vert_id,vert_ref_id, pass_id;
    double x_cent=0, y_cent=0, z_cent=0;
    double *dxdydz = new double[3];
    dxdydz[0]=0; dxdydz[1]=0; dxdydz[2]=0;
    
    //vert_ref_id
    vert_ref_id=v_partner[basal_vertices[i][3]];
    x_cent+=v[vert_ref_id][1];
    y_cent+=v[vert_ref_id][2];
    z_cent+=v[vert_ref_id][3];
    
    //vert_id
    for(size_t j=4; j <= 2+basal_edges[i][2]; j++){
        vert_id=v_partner[basal_vertices[i][j]];
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
        v_pass_type[pass_id]=2;//v_pass_type
        v_pass_cell[pass_id]=i;
        c_cent_apical[i]=pass_id;
    }
    //PLACE PASSIVE CENTER VERTEX
    else{
        pass_id=c_cent_apical[i];
        v_pass[pass_id][1]=x_cent;
        v_pass[pass_id][2]=y_cent;
        v_pass[pass_id][3]=z_cent;
        torus_vertex_pass(pass_id);
    }
    
    delete []dxdydz;
    
    return pass_id;
    
}
//****************************************************************************
void make_apical_side(int i){
    
    //CALCULATE CELL CENTER
    int pass_id=apical_center(i,1);
    
    //MAKE TRIANGULAR FACETS
    apical_facets[i][1]=basal_edges[i][1];
    apical_facets[i][2]=basal_edges[i][2];
    for(size_t j=3; j < 2+basal_edges[i][2]; j++){
        apical_facets[i][j]=make_facet(
                         v_partner[basal_vertices[i][j]],
                         v_partner[basal_vertices[i][j+1]],
                         pass_id
                         );
        f_type[apical_facets[i][j]]=2;//f_type
        f_T[apical_facets[i][j]]=c_alpha[i];//surface tension
        f_cell[apical_facets[i][j]]=i;//f_cell
    }
    int j=2+basal_edges[i][2];
    apical_facets[i][j]=make_facet(
                     v_partner[basal_vertices[i][j]],
                     v_partner[basal_vertices[i][3]],
                     pass_id
                     );
    f_type[apical_facets[i][j]]=2;//f_type
    f_T[apical_facets[i][j]]=c_alpha[i];//surface tension
    f_cell[apical_facets[i][j]]=i;//f_cell
    
}
//****************************************************************************
void make_apical_side_ALL(){
    for(int i=1; i<=Nc; i++) make_apical_side(i);
}
//****************************************************************************
void set_apical_tension(int i, double tens){
    c_alpha[i]=tens;
    for(size_t j=3; j <= 2+basal_edges[i][2]; j++) f_T[apical_facets[i][j]]=tens;
}
//****************************************************************************
