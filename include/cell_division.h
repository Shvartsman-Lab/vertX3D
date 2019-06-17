#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//******************************CELL DIVISION*********************************
//****************************************************************************
//****************************************************************************
int new_vertex_at_edge_center(int i){
    
    int v1=e[i][1], v2=e[i][2];
   
    //CALCULATES CENTER POINT OF EDGE i
    double x_cent, y_cent, z_cent;
    //x
    if(fabs(dist_x(v1,v2))>perioXYZ[0]/2.) x_cent=0.5*(v[v1][1]+v[v2][1]+perioXYZ[0]);
    else x_cent=0.5*(v[v1][1]+v[v2][1]);
    //y
    if(fabs(dist_y(v1,v2))>perioXYZ[1]/2.) y_cent=0.5*(v[v1][2]+v[v2][2]+perioXYZ[1]);
    else y_cent=0.5*(v[v1][2]+v[v2][2]);
    //z
    if(fabs(dist_z(v1,v2))>perioXYZ[2]/2.) z_cent=0.5*(v[v1][3]+v[v2][3]+perioXYZ[2]);
    else z_cent=0.5*(v[v1][3]+v[v2][3]);
    
    //CREATES VERTEX
    int vert_id=make_vertex(x_cent,y_cent,z_cent);
    return vert_id;
}
//****************************************************************************
void ReadFromLists_division(int i, int eid1, int eid2, int *elements){
    
    //IDENTIFIES EDGES e1 & e2
    int e1=abs(basal_edges[i][2+eid1]);
    int e2=abs(basal_edges[i][2+eid2]);
    
    //IDENTIFIES CELLS c1 & c2
    int c1=0,c2=0;
    //c1
    if(e_cell1[e1]!=i) c1=e_cell1[e1];
    else c1=e_cell2[e1];
    //c2
    if(e_cell1[e2]!=i) c2=e_cell1[e2];
    else c2=e_cell2[e2];
    
    //IDENTIFIES VERTICES vv1, vvv1, vv2 & vvv2
    int vv1=e[e1][1];
    int vvv1=e[e1][2];
    int vv2=e[e2][1];
    int vvv2=e[e2][2];
    
    elements[1]=e1;  elements[2]=e2;
    elements[3]=c1;  elements[4]=c2;
    elements[5]=vv1; elements[6]=vvv1; elements[7]=vv2; elements[8]=vvv2;
    
    //printf("e1=%d  e2=%d\n", e1, e2);
    //printf("c1=%d  c2=%d\n", c1, c2);
    //printf("vv1=%d  vvv1=%d  vv2=%d  vvv2=%d\n", vv1, vvv1, vv2, vvv2);
    
}
//****************************************************************************
int divide(int i, int eid1, int eid2){
    
    //printf("\n //divide// \ni=%d\n", i);
    
    //FINDS ELEMENTS ASSOCIATED WITH TRANSFORMATION
    int *elements; elements = new int[9];
    ReadFromLists_division(i,eid1,eid2,elements);
    int e1=elements[1];  int e2=elements[2];
    int c1=elements[3];  int c2=elements[4];
    int vv1=elements[5]; int vvv1=elements[6]; int vv2=elements[7]; int vvv2=elements[8];
    delete [] elements;
    
    //DISSOLVES CELLS i, c1, c2 & c3
    dissolve_cell(i);
    if (c1!=0) dissolve_cell(c1);
    if (c2!=0) dissolve_cell(c2);
    
    //DISSOLVES LATERAL SIDES e1, e2, e3 & e4
    dissolve_lateral(e1);
    dissolve_lateral(e2);
    
    //MAKES NEW VERTICES v1 & v2
    //v1
    set_normal_vector(vv1);
    set_normal_vector(vvv1);
    int v1=new_vertex_at_edge_center(e1);
    v_normal_vector[v1][1]=0.5*(v_normal_vector[vv1][1]+v_normal_vector[vvv1][1]);
    v_normal_vector[v1][2]=0.5*(v_normal_vector[vv1][2]+v_normal_vector[vvv1][2]);
    v_normal_vector[v1][3]=0.5*(v_normal_vector[vv1][3]+v_normal_vector[vvv1][3]);
    double norm=sqrt(v_normal_vector[v1][1]*v_normal_vector[v1][1]+
              v_normal_vector[v1][2]*v_normal_vector[v1][2]+
              v_normal_vector[v1][3]*v_normal_vector[v1][3]);
    v_normal_vector[v1][1]/=norm;
    v_normal_vector[v1][2]/=norm;
    v_normal_vector[v1][3]/=norm;
    v_height[v1]=0.5*(v_height[vv1]+v_height[vvv1]);
    //v2
    set_normal_vector(vv2);
    set_normal_vector(vvv2);
    int v2=new_vertex_at_edge_center(e2);
    v_normal_vector[v2][1]=0.5*(v_normal_vector[vv2][1]+v_normal_vector[vvv2][1]);
    v_normal_vector[v2][2]=0.5*(v_normal_vector[vv2][2]+v_normal_vector[vvv2][2]);
    v_normal_vector[v2][3]=0.5*(v_normal_vector[vv2][3]+v_normal_vector[vvv2][3]);
    norm=sqrt(v_normal_vector[v2][1]*v_normal_vector[v2][1]+
              v_normal_vector[v2][2]*v_normal_vector[v2][2]+
              v_normal_vector[v2][3]*v_normal_vector[v2][3]);
    v_normal_vector[v2][1]/=norm;
    v_normal_vector[v2][2]/=norm;
    v_normal_vector[v2][3]/=norm;
    v_height[v2]=0.5*(v_height[vv2]+v_height[vvv2]);
    
    //COPIES VERTICES v1 & v2 TO APICAL
    copy_vertex_to_apical(v1);
    copy_vertex_to_apical(v2);
    
    //RESTITCHES EDGES e1 & e2
    e[e1][1]=vv1;
    e[e1][2]=v1;
    e[e2][1]=vv2;
    e[e2][2]=v2;
    
    //MAKES EDGES ee1 & ee2
    int ee1=make_edge(v1,vvv1);
    int ee2=make_edge(v2,vvv2);
    
    //MAKES EDGE newe
    int newe=make_edge(v1,v2);
    
    //REMAKES LATERAL SIDES OVER e1, e2, ee1, ee2 & newe
    make_lateral_side(e1);
    make_lateral_side(e2);
    make_lateral_side(ee1);
    make_lateral_side(ee2);
    make_lateral_side(newe);
    
    //DIVIDE CELL IN EDGE NETWORK
    int newc=divide_cell_in_edge_network(i,e1,e2,ee1,ee2,newe);

    //REMAKES CELLS c1, c2, i & newc
    //i
    make_basal_side(i);
    make_apical_side(i);
    //newc
    make_basal_side(newc);
    make_apical_side(newc);
    
    //RESTITCH BASAL NETWORK FOR c1 & c2
    //c1
    if (c1!=0) replace_one_edge_with_two_edges(c1,e1,ee1);
    //c2
    if (c2!=0) replace_one_edge_with_two_edges(c2,e2,ee2);
    
    //c1
    if (c1!=0) make_basal_side(c1);
    if (c1!=0) make_apical_side(c1);
    //c2
    if (c2!=0) make_basal_side(c2);
    if (c2!=0) make_apical_side(c2);
    
    //printf("v1=%d  v2=%d\n", v1, v2);
    //printf("ee1=%d  ee2=%d\n", ee1, ee2);
    //printf("newe=%d\n", newe);
    //printf("newc=%d\n", newc);
    
    return newc;
}
//****************************************************************************
int cell_division(int i){
        
    int newc;
    int eid1, eid2;
    double cvol;
    
    if(basal_edges[i][1]!=0){
        
        cvol=cell_volume(i,f);
        
        printf("cell %d divides\n", i);
        
        if(basal_edges[i][2]%2==0){
            eid1=rnd_int(basal_edges[i][2]/2-1);
            eid2=eid1+basal_edges[i][2]/2;
            newc=divide(i,eid1,eid2);
        }
        else {
            eid1=rnd_int((basal_edges[i][2]-1)/2);
            eid2=eid1+(basal_edges[i][2]+rnd_H())/2;
            newc=divide(i,eid1,eid2);
        }
        
        //CORRECTS v_T1dir of 4-WAY VERTICES OF THE CELL newc
        for(size_t j=3; j<=2+basal_edges[newc][2]; j++){
            if( v_cell4[basal_vertices[newc][j]]!=0 && v_T1dir[basal_vertices[newc][j]]==i ) v_T1dir[basal_vertices[newc][j]]=newc;
        }
        
        //RESETS CELL ATTRIBUTES RELATED TO GROWTH
        c_V0[i]=V0;
        c_kV[i]=kV;
    }
    
    return 1;
}
//****************************************************************************
