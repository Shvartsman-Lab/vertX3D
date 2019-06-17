#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//*****************************DISTANCES**************************************
//****************************************************************************
//****************************************************************************
double edge_length(int i){
    
    double ddx,ddy,ddz;
    double *dxdydz = new double[3];
    dxdydz[0]=0; dxdydz[1]=0; dxdydz[2]=0;
    
    //basal
    int v1=e[i][1], v2=e[i][2];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
    ddz=v[v2][3]-(v[v1][3]+dxdydz[2]);
    double lb=sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
    
    //apical
    v1=v_partner[e[i][1]], v2=v_partner[e[i][2]];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
    ddz=v[v2][3]-(v[v1][3]+dxdydz[2]);
    double la=sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
    
    delete []dxdydz;
    
    return 0.5*(la+lb);
}
//****************************************************************************
void output_edge_lengths(char *fileName){
    
    FILE *file1;
    file1 = fopen(fileName, "wt");
    for(size_t i = 1; i<=Ne; ++i){
        if(e[i][0]!=0){
            fprintf(file1, "%g\n", edge_length(i));
        }
    }
    fclose(file1);
    
}
//****************************************************************************
double dist(int i, int j){
    double dx=v[i][1]-v[j][1];
    double dy=v[i][2]-v[j][2];
    double dz=v[i][3]-v[j][3];
    return sqrt(dx*dx+dy*dy+dz*dz);
}
//****************************************************************************
double dist_x(int i, int j){
    return v[j][1]-v[i][1];
}
//****************************************************************************
double dist_y(int i, int j){
    return v[j][2]-v[i][2];
}
//****************************************************************************
double dist_z(int i, int j){
    return v[j][3]-v[i][3];
}
//****************************************************************************
