#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***********************VERT, EDG, FAC***************************************
//****************************************************************************
//****************************************************************************
int make_vertex(double x, double y, double z){
    
    int v_id=v_freeId.get(Nv+1);
    v[v_id][0]=v_id;
    v[v_id][1]=x;
    v[v_id][2]=y;
    v[v_id][3]=z;
    torus_vertex(v_id);
    
    //attributes
    v_type[v_id]=1;
    
    if(v_id>Nv) Nv++;
    
    return v_id;
}
//****************************************************************************
int make_vertex_pass(double x, double y, double z){
    
    int v_id=v_pass_freeId.get(Nv_pass+1);
    v_pass[v_id][0]=v_id;
    v_pass[v_id][1]=x;
    v_pass[v_id][2]=y;
    v_pass[v_id][3]=z;
    
    //attributes
    v_pass_type[v_id]=1;
    
    if(v_id>Nv_pass) Nv_pass++;
    
    return v_id;
}
//****************************************************************************
int make_edge(int v1, int v2){
    
    int e_id=e_freeId.get(Ne+1);
    e[e_id][0]=e_id;
    e[e_id][1]=v1;
    e[e_id][2]=v2;
    
    //attributes
    e_length[e_id]=edge_length(e_id);
    
    if(e_id>Ne) Ne++;
    
    return e_id;
}
//****************************************************************************
int make_facet(int v1, int v2, int v3){
    
    int f_id=f_freeId.get(Nf+1);
    f[f_id][0]=f_id;
    f[f_id][1]=v1;
    f[f_id][2]=v2;
    f[f_id][3]=v3;
    
    //attributes
    f_type[f_id]=1;
    
    if(f_id>Nf) Nf++;
    
    return f_id;
}
//****************************************************************************
