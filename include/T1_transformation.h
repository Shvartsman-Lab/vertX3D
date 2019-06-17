#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//*********************TOPOLOGICAL TRANSITIONS********************************
//****************************************************************************
//****************************************************************************
void shrink_edge(int v1, int v2){
    
    //CALCULATES CENTER POINT (NEW POSITION OF VERTEX v1)
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
    
    //PUTS VERTEX v1 TO (x_cent,y_cent,z_cent)
    v[v1][1]=x_cent;
    v[v1][2]=y_cent;
    v[v1][3]=z_cent;
    torus_vertex(v1);
    
    //CALCULATES AND SETS NORMAL VECTOR AND HEIGHT TO VERTEX v1
    v_normal_vector[v1][1]=0.5*(v_normal_vector[v1][1]+v_normal_vector[v2][1]);
    v_normal_vector[v1][2]=0.5*(v_normal_vector[v1][2]+v_normal_vector[v2][2]);
    v_normal_vector[v1][3]=0.5*(v_normal_vector[v1][3]+v_normal_vector[v2][3]);
    double norm=sqrt(v_normal_vector[v1][1]*v_normal_vector[v1][1]+
                     v_normal_vector[v1][2]*v_normal_vector[v1][2]+
                     v_normal_vector[v1][3]*v_normal_vector[v1][3]);
    v_normal_vector[v1][1]/=norm;
    v_normal_vector[v1][2]/=norm;
    v_normal_vector[v1][3]/=norm;
    v_height[v1]=0.5*(v_height[v1]+v_height[v2]);
    
    //REPOSITIONS VERTEX v_partner[v1]
    int vert_id=v_partner[v1];
    v[vert_id][1]=v[v1][1]+v_height[v1]*v_normal_vector[v1][1];
    v[vert_id][2]=v[v1][2]+v_height[v1]*v_normal_vector[v1][2];
    v[vert_id][3]=v[v1][3]+v_height[v1]*v_normal_vector[v1][3];
    torus_vertex(vert_id);
    
}
//****************************************************************************
void ReadFromLists_merge_vertices(int i, int *elements){
    
    //FIGURES OUT VERTICES v1 & v2
    int v1=e[i][1], v2=e[i][2];
    
    //FIGURES OUT CELLS c1 & c2
    int c1=0, c2=0;
    c1=e_cell1[i];
    c2=e_cell2[i];
    
    //FIGURES OUT CELLS c3 & c4
    int c3=0, c4=0;
    //c3
    if     (v_cell1[v1]!=c1 && v_cell1[v1]!=c2) c3=v_cell1[v1];
    else if(v_cell2[v1]!=c1 && v_cell2[v1]!=c2) c3=v_cell2[v1];
    else if(v_cell3[v1]!=c1 && v_cell3[v1]!=c2) c3=v_cell3[v1];
    //c4
    if     (v_cell1[v2]!=c1 && v_cell1[v2]!=c2) c4=v_cell1[v2];
    else if(v_cell2[v2]!=c1 && v_cell2[v2]!=c2) c4=v_cell2[v2];
    else if(v_cell3[v2]!=c1 && v_cell3[v2]!=c2) c4=v_cell3[v2];
    
    //FIGURES OUT EDGES e1, e2, e3 & e4
    int edge_id, cell_id;
    int e1=0, e2=0, e3=0, e4=0;
    //e1
    if (c1!=0){
        cell_id=c1;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && edge_id!=i) e1=edge_id;
        }
    }
    else if (c3!=0){
        cell_id=c3;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && (e_cell1[edge_id]!=c2 && e_cell2[edge_id]!=c2)) e1=edge_id;
        }
        
    }
    //e2
    if (c2!=0){
        cell_id=c2;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && edge_id!=i) e2=edge_id;
        }
    }
    else if (c3!=0){
        cell_id=c3;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && (e_cell1[edge_id]!=c1 && e_cell2[edge_id]!=c1)) e2=edge_id;
        }
        
    }
    //e3
    if (c1!=0){
        cell_id=c1;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v2 || e[edge_id][2]==v2) && edge_id!=i) e3=edge_id;
        }
    }
    else if (c4!=0){
        cell_id=c4;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v2 || e[edge_id][2]==v2) && (e_cell1[edge_id]!=c2 && e_cell2[edge_id]!=c2)) e3=edge_id;
        }
        
    }
    //e4
    if (c2!=0){
        cell_id=c2;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v2 || e[edge_id][2]==v2) && edge_id!=i) e4=edge_id;
        }
    }
    else if (c4!=0){
        cell_id=c4;
        for(int j=3; j<=2+basal_edges[cell_id][2]; j++){
            edge_id=abs(basal_edges[cell_id][j]);
            if((e[edge_id][1]==v2 || e[edge_id][2]==v2) && (e_cell1[edge_id]!=c1 && e_cell2[edge_id]!=c1)) e4=edge_id;
        }
        
    }
    
    elements[1]=v1; elements[2]=v2;
    elements[3]=e1; elements[4]=e2; elements[5]=e3; elements[6]=e4;
    elements[7]=c1; elements[8]=c2; elements[9]=c3; elements[10]=c4;
    
    //printf("v1=%d v2=%d\n", v1, v2);
    //printf("c1=%d c2=%d c3=%d c4=%d\n", c1, c2, c3, c4);
    //printf("e1=%d e2=%d e3=%d e4=%d\n", e1, e2, e3, e4);
}
//****************************************************************************
int merge_vertices(int i){//vertex v1 stays for valence_reduction!
    
    //FINDS ELEMENTS ASSOCIATED WITH TRANSFORMATION
    int *elements; elements = new int[11];
    ReadFromLists_merge_vertices(i,elements);
    int v1=elements[1]; int v2=elements[2];
    int e1=elements[3]; int e2=elements[4]; int e3=elements[5]; int e4=elements[6];
    int c1=elements[7]; int c2=elements[8]; int c3=elements[9]; int c4=elements[10];
    delete [] elements;
    
    //CHECKS IF VALENCE REDUCTION IS NEEDED AFTER merge_vertices
    if (vertex_valence(v1)==3 && vertex_valence(v2)==3){
        
        //CELL DIRECTION
        if      (c1!=0) v_T1dir[v1]=c1;
        else if (c2!=0) v_T1dir[v1]=c2;
    
        //RESERVES FOR VALENCE REDUCTION
        v_vertT1[v1]=v2;
        v_edgeT1[v1]=i;
    }
    
    //DISSOLVES CELLS c1, c2 & c4
    if (c1!=0) dissolve_cell(c1);
    if (c2!=0) dissolve_cell(c2);
    if (c3!=0) dissolve_cell(c3);
    if (c4!=0) dissolve_cell(c4);
    
    //DISSOLVES LATERAL SIDE
    dissolve_lateral(i);
    if (e3!=0) dissolve_lateral(e3);
    if (e4!=0) dissolve_lateral(e4);
    
    //DISSOLVES EDGE i
    dissolve_edge(i);
    
    //SHRINKS EDGE
    set_normal_vector(v1);
    set_normal_vector(v2);
    shrink_edge(v1,v2);
    
    //DISSOLVES VERTEX v2 & HIS APICAL PARTNER
    dissolve_vertex(v_partner[v2]);
    dissolve_vertex(v2);
    
    //RESTITCHES EDGES e3 & e4
    //e3
    if (e3!=0){
        if(e[e3][1]==v2) e[e3][1]=v1;
        else if(e[e3][2]==v2) e[e3][2]=v1;
    }
    //e4
    if (e4!=0){
        if(e[e4][1]==v2) e[e4][1]=v1;
        else if(e[e4][2]==v2) e[e4][2]=v1;
    }
    
    //REMAKES LATERAL SIDES OVER e3 & e4
    if (e3!=0) make_lateral_side(e3);
    if (e4!=0) make_lateral_side(e4);
    
    //REMOVES EDGE i FROM basal_edges
    if (c1!=0) remove_edge_from_basal_edges(i,c1);
    if (c2!=0) remove_edge_from_basal_edges(i,c2);
    
    //REMAKES CELLS c1, c2, c3 & c4
    //c1
    if (c1!=0){
        make_basal_side(c1);
        make_apical_side(c1);
    }
    //c2
    if (c2!=0){
        make_basal_side(c2);
        make_apical_side(c2);
    }
    //c3
    if (c3!=0){
        make_basal_side(c3);
        make_apical_side(c3);
    }
    //c4
    if (c4!=0){
        make_basal_side(c4);
        make_apical_side(c4);
    }
    
    return v1;
    
}
//****************************************************************************
void ReadFromLists_valence_reduction(int v1, int *elements){
    
    //v2
    int v2=v_vertT1[v1];
    
    //c1
    int c1=v_T1dir[v1];
    
    //FIGURES OUT EDGES e1 & e2
    int e1=0, e2=0;
    int edge_id;
    for(int j=3; j<=2+basal_edges[c1][2]; j++){
        edge_id=abs(basal_edges[c1][j]);
        if(e[edge_id][1]==v1 || e[edge_id][2]==v1){
            if(e1==0) e1=edge_id;
            else{
                e2=edge_id;
                break;
            }
        }
    }
    
    //FIGURES OUT CELLS c2, c3 & c4
    int c2=0, c3=0, c4=0;
    //c3
    if(e_cell1[e1]==c1) c3=e_cell2[e1];
    if(e_cell2[e1]==c1) c3=e_cell1[e1];
    //c4
    if(e_cell1[e2]==c1) c4=e_cell2[e2];
    if(e_cell2[e2]==c1) c4=e_cell1[e2];
    //c2
    if      (v_cell1[v1]!=c1 && v_cell1[v1]!=c3 && v_cell1[v1]!=c4) c2=v_cell1[v1];
    else if (v_cell2[v1]!=c1 && v_cell2[v1]!=c3 && v_cell2[v1]!=c4) c2=v_cell2[v1];
    else if (v_cell3[v1]!=c1 && v_cell3[v1]!=c3 && v_cell3[v1]!=c4) c2=v_cell3[v1];
    else if (v_cell4[v1]!=c1 && v_cell4[v1]!=c3 && v_cell4[v1]!=c4) c2=v_cell4[v1];
    
    //FIGURES OUT EDGES e3 & e4
    int e3=0, e4=0;
    if (c2!=0){
        for(int j=3; j<=2+basal_edges[c2][2]; j++){
            edge_id=abs(basal_edges[c2][j]);
            if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && (e_cell1[edge_id]!=c4 && e_cell2[edge_id]!=c4)) e3=edge_id;
            if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && (e_cell1[edge_id]!=c3 && e_cell2[edge_id]!=c3)) e4=edge_id;
            if(e3!=0 && e4!=0) break;
        }
    }
    else{
        if(c3!=0){
            for(int j=3; j<=2+basal_edges[c3][2]; j++){
                edge_id=abs(basal_edges[c3][j]);
                if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && (e_cell1[edge_id]!=c1 && e_cell2[edge_id]!=c1)) {
                    e3=edge_id;
                    break;
                }
            }
        }
        if(c4!=0){
            for(int j=3; j<=2+basal_edges[c4][2]; j++){
                edge_id=abs(basal_edges[c4][j]);
                if((e[edge_id][1]==v1 || e[edge_id][2]==v1) && (e_cell1[edge_id]!=c1 && e_cell2[edge_id]!=c1)) {
                    e4=edge_id;
                    break;
                }
            }
        }
    }

    elements[2]=v2;
    elements[3]=e1; elements[4]=e2; elements[5]=e3; elements[6]=e4;
    elements[7]=c1; elements[8]=c2; elements[9]=c3; elements[10]=c4;
    
    //printf("v1=%d v2=%d\n", v1, v2);
    //printf("c1=%d c2=%d c3=%d c4=%d\n", c1, c2, c3, c4);
    //printf("e1=%d e2=%d e3=%d e4=%d\n\n", e1, e2, e3, e4);
}
//****************************************************************************
int valence_reduction(int v1, double fin_len){
    
    //i
    int i=v_edgeT1[v1];
    
    //FINDS ELEMENTS ASSOCIATED WITH TRANSFORMATION
    int *elements; elements = new int[11];
    ReadFromLists_valence_reduction(v1,elements);
    int v2=elements[2];
    int e1=elements[3]; int e2=elements[4]; int e3=elements[5]; int e4=elements[6];
    int c1=elements[7]; int c2=elements[8]; int c3=elements[9]; int c4=elements[10];
    delete [] elements;
    
    //CALCULATES DIRECITONS OF CELLS c1 & c2
    double *cdir; cdir=new double[4];
    //x
    if(fabs(v_pass[c_cent_basal[c1]][1]-v[v1][1])>perioXYZ[0]/2.){
        if(v_pass[c_cent_basal[c1]][1]>v[v1][1]) cdir[1]=v_pass[c_cent_basal[c1]][1]-(v[v1][1]+perioXYZ[0]);
        if(v_pass[c_cent_basal[c1]][1]<v[v1][1]) cdir[1]=v_pass[c_cent_basal[c1]][1]-(v[v1][1]-perioXYZ[0]);
    }
    else cdir[1]=v_pass[c_cent_basal[c1]][1]-v[v1][1];
    //y
    if(fabs(v_pass[c_cent_basal[c1]][2]-v[v1][2])>perioXYZ[1]/2.){
        if(v_pass[c_cent_basal[c1]][2]>v[v1][2]) cdir[2]=v_pass[c_cent_basal[c1]][2]-(v[v1][2]+perioXYZ[1]);
        if(v_pass[c_cent_basal[c1]][2]<v[v1][2]) cdir[2]=v_pass[c_cent_basal[c1]][2]-(v[v1][2]-perioXYZ[1]);
    }
    else cdir[2]=v_pass[c_cent_basal[c1]][2]-v[v1][2];
    //z
    if(fabs(v_pass[c_cent_basal[c1]][3]-v[v1][3])>perioXYZ[2]/2.){
        if(v_pass[c_cent_basal[c1]][3]>v[v1][3]) cdir[3]=v_pass[c_cent_basal[c1]][3]-(v[v1][3]+perioXYZ[2]);
        if(v_pass[c_cent_basal[c1]][3]<v[v1][3]) cdir[3]=v_pass[c_cent_basal[c1]][3]-(v[v1][3]-perioXYZ[2]);
    }
    else cdir[3]=v_pass[c_cent_basal[c1]][3]-v[v1][3];
    double norm=sqrt(cdir[1]*cdir[1]+cdir[2]*cdir[2]+cdir[3]*cdir[3]);
    
    //DISSOLVES CELLS c1, c2, c3 & c4
    if (c1!=0) dissolve_cell(c1);
    if (c2!=0) dissolve_cell(c2);
    if (c3!=0) dissolve_cell(c3);
    if (c4!=0) dissolve_cell(c4);
    
    //DISSOLVES LATERAL SIDES e1, e2, e3 & e4
    if (e1!=0) dissolve_lateral(e1);
    if (e2!=0) dissolve_lateral(e2);
    if (e3!=0) dissolve_lateral(e3);
    if (e4!=0) dissolve_lateral(e4);
    
    //MAKES VERTEX v2
    v2=make_vertex(v[v1][1]+0.5*fin_len*cdir[1]/norm,
                v[v1][2]+0.5*fin_len*cdir[2]/norm,
                v[v1][3]+0.5*fin_len*cdir[3]/norm
                );
    set_normal_vector(v1);
    v_height[v2]=v_height[v1];
    v_normal_vector[v2][1]=v_normal_vector[v1][1];
    v_normal_vector[v2][2]=v_normal_vector[v1][2];
    v_normal_vector[v2][3]=v_normal_vector[v1][3];
    
    //COPIES VERTEX v2 TO APICAL
    copy_vertex_to_apical(v2);
    
    //RESTITCHES EDGES e1 & e2
    //e1
    if (e1!=0){
        if(e[e1][1]==v1) e[e1][1]=v2;
        else e[e1][2]=v2;
    }
    //e2
    if (e2!=0){
        if(e[e2][1]==v1) e[e2][1]=v2;
        else e[e2][2]=v2;
    }
    
    //MOVES VERTEX v1
    double v1_x=v[v1][1]-0.5*fin_len*cdir[1]/norm;
    double v1_y=v[v1][2]-0.5*fin_len*cdir[2]/norm;
    double v1_z=v[v1][3]-0.5*fin_len*cdir[3]/norm;
    delete [] cdir;
    dissolve_vertex(v_partner[v1]);
    dissolve_vertex(v1);//1-reserved
    v1=make_vertex(v1_x,v1_y,v1_z);
    v_height[v1]=v_height[v2];
    v_normal_vector[v1][1]=v_normal_vector[v2][1];
    v_normal_vector[v1][2]=v_normal_vector[v2][2];
    v_normal_vector[v1][3]=v_normal_vector[v2][3];
    copy_vertex_to_apical(v1);
    
    //MAKES EDGE i
    i=make_edge(v1,v2);
    
    //INSERTS EDGE i IN CELLS c3 & c4
    if (c3!=0 && e1!=0 && e3!=0) insert_edge_to_basal_edges(i,c3,e1,e3);
    if (c4!=0 && e2!=0 && e4!=0) insert_edge_to_basal_edges(i,c4,e2,e4);
    
    //REMAKES LATERAL SIDES OVER e1, e2, e3, e4 & i
    if (e1!=0) make_lateral_side(e1);
    if (e2!=0) make_lateral_side(e2);
    if (e3!=0) make_lateral_side(e3);
    if (e4!=0) make_lateral_side(e4);
    make_lateral_side(i);
    
    //REMAKES CELLS c1, c2, c3 & c4
    //c1
    if (c1!=0){
        make_basal_side(c1);
        make_apical_side(c1);
    }
    //c2
    if (c2!=0){
        make_basal_side(c2);
        make_apical_side(c2);
    }
    //c3
    if (c3!=0){
        make_basal_side(c3);
        make_apical_side(c3);
    }
    //c4
    if (c4!=0){
        make_basal_side(c4);
        make_apical_side(c4);
    }
    
    return 1;
    
}
//****************************************************************************
int T1(int i, double fin_len){
    
    if( v_cell4[e[i][1]]==0 && v_cell4[e[i][2]]==0 ){
        
        if( e_cell1[i]!=0 && e_cell2[i]!=0 ){
            if( basal_edges[e_cell1[i]][2]==3  || basal_edges[e_cell2[i]][2]==3 ) return 0;
            else {
                printf("T1 on edge %d\n", i);
                int v1=merge_vertices(i);
                if(v_edgeT1[v1]!=0) valence_reduction(v1,fin_len);
                return 1;
            }
        }
        
        else if( e_cell1[i]==0 && e_cell2[i]!=0 ){
            if( basal_edges[e_cell2[i]][2]==3 ) return 0;
            else {
                printf("T1 on edge %d\n", i);
                int v1=merge_vertices(i);
                if(v_edgeT1[v1]!=0) valence_reduction(v1,fin_len);
                return 1;
            }
        }
        
        else if( e_cell1[i]!=0 && e_cell2[i]==0 ){
            if( basal_edges[e_cell1[i]][2]==3 ) return 0;
            else {
                printf("T1 on edge %d\n", i);
                int v1=merge_vertices(i);
                if(v_edgeT1[v1]!=0) valence_reduction(v1,fin_len);
                return 1;
            }
        }
        else return 0;
    }
    
    else return 0;
}
//****************************************************************************
