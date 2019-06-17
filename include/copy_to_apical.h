#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//*************************COPY TO APICAL*************************************
//****************************************************************************
//****************************************************************************
void set_normal_vector(int i){
    
    int iA=v_partner[i];
    
    //x
    double dx;
    double distx=dist_x(i,iA);
    if(fabs(distx)>0.5*perioXYZ[0]){
        if(distx>0) dx=perioXYZ[0];
        else dx=-perioXYZ[0];
    }
    else dx=0;
    
    //y
    double dy;
    double disty=dist_y(i,iA);
    if(fabs(disty)>0.5*perioXYZ[1]){
        if(disty>0) dy=perioXYZ[1];
        else dy=-perioXYZ[1];
    }
    else dy=0;
    
    //z
    double dz;
    double distz=dist_z(i,iA);
    if(fabs(distz)>0.5*perioXYZ[2]){
        if(distz>0) dz=perioXYZ[2];
        else dz=-perioXYZ[2];
    }
    else dz=0;
    
    double nx=v[iA][1]-(v[i][1]+dx);
    double ny=v[iA][2]-(v[i][2]+dy);
    double nz=v[iA][3]-(v[i][3]+dz);
    
    v_height[i]=sqrt(nx*nx+ny*ny+nz*nz);
    v_normal_vector[i][1]=nx/v_height[i];
    v_normal_vector[i][2]=ny/v_height[i];
    v_normal_vector[i][3]=nz/v_height[i];
}
//****************************************************************************
int copy_vertex_to_apical(int vert_id){
    
    int vid=make_vertex(
                    v[vert_id][1]+v_height[vert_id]*v_normal_vector[vert_id][1],
                    v[vert_id][2]+v_height[vert_id]*v_normal_vector[vert_id][2],
                    v[vert_id][3]+v_height[vert_id]*v_normal_vector[vert_id][3]
                    );
    
    //v_type
    v_type[vid]=2;
    
    //v_partner
    v_partner[vert_id]=vid;
    v_partner[vid]=vert_id;
    
    return vid;
}
//****************************************************************************
void copy_vertex_to_apical_ALL(){
    int i;
    for(i=1; i<=Nv; i++){
        if(v_type[i]==1) copy_vertex_to_apical(i);
    }
}
//****************************************************************************