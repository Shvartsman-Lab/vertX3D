#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//******************************CELL EXTRUSION********************************
//****************************************************************************
//****************************************************************************
void ReadFromLists_T2(int i, int *elements){

    //IDENTIFIES EDGES e1, e2 & e3
    int e1=0,e2=0,e3=0;
    e1=abs(basal_edges[i][3]);
    e2=abs(basal_edges[i][4]);
    e3=abs(basal_edges[i][5]);
    
    //IDENTIFIES VERTICES v1, v2 & v3
    //v1, v2
    int v1=0,v2=0,v3=0;
    if(basal_edges[i][3]>0){
        v1=e[e1][1];
        v2=e[e1][2];
    }
    else{
        v1=e[e1][2];
        v2=e[e1][1];
    }
    //v3
    if(basal_edges[i][4]>0) v3=e[e2][2];
    else v3=e[e2][1];
    
    //IDENTIFIES CELLS c1, c2 & c3
    int c1=0,c2=0,c3=0;
    //c1
    if(e_cell1[e1]!=i) c1=e_cell1[e1];
    else c1=e_cell2[e1];
    //c2
    if(e_cell1[e2]!=i) c2=e_cell1[e2];
    else c2=e_cell2[e2];
    //c3
    if(e_cell1[e3]!=i) c3=e_cell1[e3];
    else c3=e_cell2[e3];
    
    //IDENTIFIES EDGES ee1, ee2 & ee3
    int ee1=0,ee2=0,ee3=0; int edg_id;
    //ee1
    if(c1==0 && c2==0) ee1=0;
    else if ( (c1!=0 && c2==0) || (c1!=0 && c2!=0) ){
        for(int j=3; j<=2+basal_edges[c1][2]; j++){
            edg_id=abs(basal_edges[c1][j]);
            if(edg_id!=e1) if( ( e_cell1[edg_id]==c1 && (e[edg_id][1]==v2 || e[edg_id][2]==v2) ) || ( e_cell2[edg_id]==c1 && (e[edg_id][1]==v2 || e[edg_id][2]==v2) ) ) ee1=edg_id;
        }
    }
    else if (c1==0 && c2!=0){
        for(int j=3; j<=2+basal_edges[c2][2]; j++){
            edg_id=abs(basal_edges[c2][j]);
            if(edg_id!=e2) if( ( e_cell1[edg_id]==c2 && (e[edg_id][1]==v2 || e[edg_id][2]==v2) ) || ( e_cell2[edg_id]==c2 && (e[edg_id][1]==v2 || e[edg_id][2]==v2) ) ) ee1=edg_id;
        }
    }
    //ee2
    if(c1==0 && c3==0) ee2=0;
    else if ( (c1!=0 && c3==0) || (c1!=0 && c3!=0) ){
        for(int j=3; j<=2+basal_edges[c1][2]; j++){
            edg_id=abs(basal_edges[c1][j]);
            if(edg_id!=e1) if( ( e_cell1[edg_id]==c1 && (e[edg_id][1]==v1 || e[edg_id][2]==v1) ) || ( e_cell2[edg_id]==c1 && (e[edg_id][1]==v1 || e[edg_id][2]==v1) ) ) ee2=edg_id;
        }
    }
    else if (c1==0 && c3!=0){
        for(int j=3; j<=2+basal_edges[c3][2]; j++){
            edg_id=abs(basal_edges[c3][j]);
            if(edg_id!=e3) if( ( e_cell1[edg_id]==c3 && (e[edg_id][1]==v1 || e[edg_id][2]==v1) ) || ( e_cell2[edg_id]==c3 && (e[edg_id][1]==v1 || e[edg_id][2]==v1) ) ) ee2=edg_id;
        }
    }
    //ee3
    if(c2==0 && c3==0) ee3=0;
    else if ( (c2!=0 && c3==0) || (c2!=0 && c3!=0) ){
        for(int j=3; j<=2+basal_edges[c2][2]; j++){
            edg_id=abs(basal_edges[c2][j]);
            if(edg_id!=e2) if( ( e_cell1[edg_id]==c2 && (e[edg_id][1]==v3 || e[edg_id][2]==v3) ) || ( e_cell2[edg_id]==c2 && (e[edg_id][1]==v3 || e[edg_id][2]==v3) ) ) ee3=edg_id;
        }
    }
    else if (c2==0 && c3!=0){
        for(int j=3; j<=2+basal_edges[c3][2]; j++){
            edg_id=abs(basal_edges[c3][j]);
            if(edg_id!=e3) if( ( e_cell1[edg_id]==c3 && (e[edg_id][1]==v3 || e[edg_id][2]==v3) ) || ( e_cell2[edg_id]==c3 && (e[edg_id][1]==v3 || e[edg_id][2]==v3) ) ) ee3=edg_id;
        }
    }
    
    elements[1]=e1;  elements[2]=e2;  elements[3]=e3;
    elements[4]=c1;  elements[5]=c2;  elements[6]=c3;
    elements[7]=ee1; elements[8]=ee2; elements[9]=ee3;
    elements[10]=v1; elements[11]=v2; elements[12]=v3;
    
    //printf("e1=%d  e2=%d  e3=%d\n", e1,e2,e3);
    //printf("c1=%d  c2=%d  c3=%d\n", c1,c2,c3);
    //printf("v1=%d  v2=%d  v3=%d\n", v1,v2,v3);
    //printf("ee1=%d  ee2=%d  ee3=%d\n", ee1,ee2,ee3);
    
    
}
//****************************************************************************
void T2(int i){
    
    //printf("T2 of cell %d\n", i);
    
    //FINDS ELEMENTS ASSOCIATED WITH TRANSFORMATION
    int *elements; elements = new int[13];
    ReadFromLists_T2(i,elements);
    int e1=elements[1];  int e2=elements[2];  int e3=elements[3];
    int c1=elements[4];  int c2=elements[5];  int c3=elements[6];
    int ee1=elements[7]; int ee2=elements[8]; int ee3=elements[9];
    int v1=elements[10]; int v2=elements[11]; int v3=elements[12];
    delete [] elements;
    
    //READS CENTER OF EXTRUDING CELL
    double x_cent=v_pass[c_cent_basal[i]][1];
    double y_cent=v_pass[c_cent_basal[i]][2];
    double z_cent=v_pass[c_cent_basal[i]][3];
    
    //DISSOLVES CELLS i, c1, c2 & c3
    dissolve_cell(i);
    if(c1!=0) dissolve_cell(c1);
    if(c2!=0) dissolve_cell(c2);
    if(c3!=0) dissolve_cell(c3);
    
    //DISSOLVES LATERAL SIDES e1, e2, e3 & e4
    dissolve_lateral(e1);
    dissolve_lateral(e2);
    dissolve_lateral(e3);
    
    //DISSOLVES LATERAL SIDES ee1, e2e & ee3
    if(ee1!=0) dissolve_lateral(ee1);
    if(ee2!=0) dissolve_lateral(ee2);
    if(ee3!=0) dissolve_lateral(ee3);
    
    //DISSOLVES EDGES e1,e2 & e3
    dissolve_edge(e1);
    dissolve_edge(e2);
    dissolve_edge(e3);
    
    //RECALCULATE NORMAL VECTOR
    set_normal_vector(v1);
    set_normal_vector(v2);
    set_normal_vector(v3);
    v_normal_vector[v1][1]=0.33333*(v_normal_vector[v1][1]+v_normal_vector[v2][1]+v_normal_vector[v3][1]);
    v_normal_vector[v1][2]=0.33333*(v_normal_vector[v1][2]+v_normal_vector[v2][2]+v_normal_vector[v3][2]);
    v_normal_vector[v1][3]=0.33333*(v_normal_vector[v1][3]+v_normal_vector[v2][3]+v_normal_vector[v3][3]);
    double norm=sqrt(v_normal_vector[v1][1]*v_normal_vector[v1][1]+
                     v_normal_vector[v1][2]*v_normal_vector[v1][2]+
                     v_normal_vector[v1][3]*v_normal_vector[v1][3]);
    v_normal_vector[v1][1]/=norm;
    v_normal_vector[v1][2]/=norm;
    v_normal_vector[v1][3]/=norm;
    v_height[v1]=0.33333*(v_height[v1]+v_height[v2]+v_height[v3]);
    
    //DISSOLVES VERTICES v2 & v3 AND THEIR APICAL PARTNERS
    //v2
    dissolve_vertex(v_partner[v2]);
    dissolve_vertex(v2);
    //v3
    dissolve_vertex(v_partner[v3]);
    dissolve_vertex(v3);
    
    //MOVES VERTEX v1 TO THE CENTER
    v[v1][1]=x_cent;
    v[v1][2]=y_cent;
    v[v1][3]=z_cent;
    
    //MOVE APICAL PARTNER
    v[v_partner[v1]][1]=v[v1][1]+v_height[v1]*v_normal_vector[v1][1];
    v[v_partner[v1]][2]=v[v1][2]+v_height[v1]*v_normal_vector[v1][2];
    v[v_partner[v1]][3]=v[v1][3]+v_height[v1]*v_normal_vector[v1][3];
    
    //RESTITCHES EDGES ee1 & ee3
    //ee1
    if(ee1!=0){
        if(e[ee1][1]==v2) e[ee1][1]=v1;
        else e[ee1][2]=v1;
    }
    //ee3
    if(ee3!=0){
        if(e[ee3][1]==v3) e[ee3][1]=v1;
        else e[ee3][2]=v1;
    }
    
    //REMAKES LATERAL SIDES OVER ee1, ee2 & ee3
    if(ee1!=0) make_lateral_side(ee1);
    if(ee2!=0) make_lateral_side(ee2);
    if(ee3!=0) make_lateral_side(ee3);
    
    //REMOVES EDGES e1, e2 & e3 FROM basal_edges
    if(c1!=0) remove_edge_from_basal_edges(e1,c1);
    if(c2!=0) remove_edge_from_basal_edges(e2,c2);
    if(c3!=0) remove_edge_from_basal_edges(e3,c3);
    
    //DELETES CELL i FROM LIST basal_edges
    delete_cell_from_basal_list(i);
    
    //REMAKES CELLS c1, c2, c3 & c4
    //c1
    if(c1!=0){
        make_basal_side(c1);
        make_apical_side(c1);
    }
    //c2
    if(c2!=0){
        make_basal_side(c2);
        make_apical_side(c2);
    }
    //c3
    if(c3!=0){
        make_basal_side(c3);
        make_apical_side(c3);
    }
}
//****************************************************************************
int cell_extrusion(int i){
    
    int edg_id;
    
    if( basal_edges[i][1]!=0 ){
        
        for(size_t j=3; j<=2+basal_edges[i][2]; j++){
            edg_id=abs(basal_edges[i][j]);
            //checks if triangles
            if      (e_cell1[edg_id]==i && e_cell2[edg_id]!=0 ){
                if (basal_edges[e_cell2[edg_id]][2]==3) return 0;
            }
            else if (e_cell2[edg_id]==i && e_cell1[edg_id]!=0 ){
                if (basal_edges[e_cell1[edg_id]][2]==3) return 0;
            }
            //checks if 4-way vertices
            if(v_cell4[basal_vertices[i][j]]!=0) return 0;
        }
        
        if(basal_edges[i][2]==3){
            printf("cell %d extrudes\n", i);
            T2(i);
        }
        
        else{
            while(basal_edges[i][2]>3){
                int j=3;
                while(j<=2+basal_edges[i][2]){
                    T1(abs(basal_edges[i][j]),0.01);
                    j++;
                }
            }
            printf("cell %d extrudes\n", i);
            T2(i);
        }
        
    }
    
    return 1;
}
//****************************************************************************
