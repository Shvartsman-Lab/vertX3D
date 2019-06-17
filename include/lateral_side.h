#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************LATERAL SIDE*************************************
//****************************************************************************
//****************************************************************************
int lateral_center(int i, int make_or_not){

    int vert1=e[i][1], vert2=e[i][2];
    int vert1a=v_partner[vert1], vert2a=v_partner[vert2];
    double *dxdydz = new double[3];
    dxdydz[0]=0; dxdydz[1]=0; dxdydz[1]=0;
    
    //CALCULATE CENTER
    double x_cent=0, y_cent=0, z_cent=0;
    double ddx=0, ddy=0, ddz=0;
    int vert_ref_id=vert1, vert_id;
    //vert2
    vert_id=vert2;
    torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
    ddx+=dxdydz[0];
    ddy+=dxdydz[1];
    ddz+=dxdydz[2];
    //vert1a
    vert_id=vert1a;
    torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
    ddx+=dxdydz[0];
    ddy+=dxdydz[1];
    ddz+=dxdydz[2];
    //vert2a
    vert_id=vert2a;
    torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
    ddx+=dxdydz[0];
    ddy+=dxdydz[1];
    ddz+=dxdydz[2];
    //CALCULATES (x_cent,y_cent,z_cent)
    x_cent=0.25*(v[vert1][1]+v[vert2][1]+v[vert1a][1]+v[vert2a][1]+ddx);
    y_cent=0.25*(v[vert1][2]+v[vert2][2]+v[vert1a][2]+v[vert2a][2]+ddy);
    z_cent=0.25*(v[vert1][3]+v[vert2][3]+v[vert1a][3]+v[vert2a][3]+ddz);
    
    int pass_id;
    //MAKE PASSIVE CENTER VERTEX
    if(make_or_not==1){
        pass_id=make_vertex_pass(x_cent,y_cent,z_cent);
        torus_vertex_pass(pass_id);
        v_pass_type[pass_id]=3;//v_pass_type
        v_pass_edge[pass_id]=i;
        e_v_pass[i]=pass_id;
    }
    //PLACE PASSIVE CENTER VERTEX
    else{
        pass_id=e_v_pass[i];
        v_pass[pass_id][1]=x_cent;
        v_pass[pass_id][2]=y_cent;
        v_pass[pass_id][3]=z_cent;
        torus_vertex_pass(pass_id);
    }
    
    delete []dxdydz;
    
    return pass_id;
}
//****************************************************************************
void make_lateral_side(int i){
    
    //CALCULATE CENTER
    int pass_id=lateral_center(i,1);
    
    //MAKE TRIANGULAR FACETS
    int facid;
    int vert1=e[i][1], vert2=e[i][2];
    int vert1a=v_partner[vert1], vert2a=v_partner[vert2];
    //1
    facid=make_facet(
                     vert1,
                     vert2,
                     pass_id
                     );
    f_type[facid]=3;//f_type
    f_T[facid]=1.;//surface tension
    f_edge[facid]=i;//f_edge
    e_lateral1[i]=facid;
    
    //2
    facid=make_facet(
                     vert2,
                     vert2a,
                     pass_id
                     );
    f_type[facid]=3;//f_type
    f_T[facid]=1.;//surface tension
    f_edge[facid]=i;//f_edge
    e_lateral2[i]=facid;
    
    //3
    facid=make_facet(
                     vert2a,
                     vert1a,
                     pass_id
                     );
    f_type[facid]=3;//f_type
    f_T[facid]=1.;//surface tension
    f_edge[facid]=i;//f_edge
    e_lateral3[i]=facid;
    
    //4
    facid=make_facet(
                     vert1a,
                     vert1,
                     pass_id
                     );
    f_type[facid]=3;//f_type
    f_T[facid]=1.;//surface tension
    f_edge[facid]=i;//f_edge
    e_lateral4[i]=facid;
    
    
}
//****************************************************************************
void make_lateral_side_ALL(){
    for(int i=1; i<=Ne; i++) make_lateral_side(i);
}
//****************************************************************************
