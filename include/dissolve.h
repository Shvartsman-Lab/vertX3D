#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************DISSOLVE*****************************************
//****************************************************************************
//****************************************************************************
void dissolve_vertex(int i){
    
    //PASSIVE VERTEX
    v[i][0]=0;
    v[i][1]=0;
    v[i][2]=0;
    v[i][3]=0;
    
    //ATTRIBUTES OF OTHER ELEMENTS
    if(v_partner[i]!=0) v_partner[v_partner[i]]=0;
    
    //ATTRIBUTES
    v_F[i][1]=0;
    v_F[i][2]=0;
    v_F[i][3]=0;
    v_normal_vector[i][1]=0;
    v_normal_vector[i][2]=0;
    v_normal_vector[i][3]=0;
    v_height[i]=0;
    v_type[i]=0;
    v_partner[i]=0;
    v_cell1[i]=0;
    v_cell2[i]=0;
    v_cell3[i]=0;
    v_cell4[i]=0;
    v_T1dir[i]=0;
    v_vertT1[i]=0;
    v_edgeT1[i]=0;
    v_clock[i]=0;
    
    v_freeId.add(i);
}
//****************************************************************************
void dissolve_vertex_pass(int i){
    
    //PASSIVE VERTEX
    v_pass[i][0]=0;
    v_pass[i][1]=0;
    v_pass[i][2]=0;
    v_pass[i][3]=0;
    
    //ATTRIBUTES OF OTHER ELEMENTS
    if(v_pass_edge[i]!=0) e_v_pass[v_pass_edge[i]]=0;
    if(v_pass_cell[i]!=0 && v_pass_type[i]==1) c_cent_basal[v_pass_cell[i]]=0;
    if(v_pass_cell[i]!=0 && v_pass_type[i]==2) c_cent_apical[v_pass_cell[i]]=0;
    
    //ATTRIBUTES
    v_pass_type[i]=0;
    v_pass_cell[i]=0;
    v_pass_edge[i]=0;
    
    v_pass_freeId.add(i);
}
//****************************************************************************
void dissolve_edge(int i){
    
    //EDGE
    e[i][0]=0;
    e[i][1]=0;
    e[i][2]=0;
    
    //ATTRIBUTES OF OTHER ELEMENTS
    if(e_v_pass[i]!=0) v_pass_edge[e_v_pass[i]]=0;
    if(e_lateral1[i]!=0) f_edge[e_lateral1[i]]=0;
    if(e_lateral2[i]!=0) f_edge[e_lateral2[i]]=0;
    if(e_lateral3[i]!=0) f_edge[e_lateral3[i]]=0;
    if(e_lateral4[i]!=0) f_edge[e_lateral4[i]]=0;
    
    //ATTRIBUTES
    e_v_pass[i]=0;
    e_lateral1[i]=0;
    e_lateral2[i]=0;
    e_lateral3[i]=0;
    e_lateral4[i]=0;
    e_cell1[i]=0;
    e_cell2[i]=0;
    e_length[i]=0;
    
    e_freeId.add(i);
}
//****************************************************************************
void dissolve_facet(int i){
    
    //FACET
    f[i][0]=0;
    f[i][1]=0;
    f[i][2]=0;
    f[i][3]=0;
    
    //ATTRIBUTES OF OTHER ELEMENTS
    if(f_edge[i]!=0){
        if(e_lateral1[f_edge[i]]==i) e_lateral1[f_edge[i]]=0;
        else if(e_lateral2[f_edge[i]]==i) e_lateral2[f_edge[i]]=0;
        else if(e_lateral3[f_edge[i]]==i) e_lateral3[f_edge[i]]=0;
        else if(e_lateral4[f_edge[i]]==i) e_lateral4[f_edge[i]]=0;
    }
    
    //ATTRIBUTES
    f_type[i]=0;
    f_cell[i]=0;
    f_edge[i]=0;
    f_T[i]=0;
    
    f_freeId.add(i);
}
//****************************************************************************
void dissolve_apical_basal_sides(int i){
    
    int pass_vert_basal=f[basal_facets[i][3]][3];
    int pass_vert_apical=f[apical_facets[i][3]][3];
    
    //DISSOLVE FACETS
    //basal
    for(size_t j=3; j<=2+basal_facets[i][2]; j++){
        dissolve_facet(basal_facets[i][j]);
    }
    //apical
    for(size_t j=3; j<=2+apical_facets[i][2]; j++){
        dissolve_facet(apical_facets[i][j]);
    }
    
    //DELETES FACETS FROM LISTS
    for(size_t j=0; j<=14; j++){
        //basal
        basal_facets[i][j]=0;
        //apical
        apical_facets[i][j]=0;
    }
    
    //PASSIVE VERTICES
    dissolve_vertex_pass(pass_vert_basal);
    dissolve_vertex_pass(pass_vert_apical);
    
}
//****************************************************************************
void dissolve_cell(int i){
    
    dissolve_apical_basal_sides(i);
    
    int vert_id, edge_id;
    
    //DELETES LIST basal_vertices ( is recreated by make_basal_side() )
    for(size_t j=3; j<=2+basal_edges[i][2]; j++){
        //v_cell
        vert_id=basal_vertices[i][j];
        if(v_cell1[vert_id]==i) v_cell1[vert_id]=0;
        if(v_cell2[vert_id]==i) v_cell2[vert_id]=0;
        if(v_cell3[vert_id]==i) v_cell3[vert_id]=0;
        if(v_cell4[vert_id]==i) v_cell4[vert_id]=0;
        basal_vertices[i][j]=0;
        //e_cell
        edge_id=abs(basal_edges[i][j]);
        if(e_cell1[edge_id]==i) e_cell1[edge_id]=0;
        if(e_cell2[edge_id]==i) e_cell2[edge_id]=0;
    }
    basal_vertices[i][1]=0;
    basal_vertices[i][2]=0;
    
    //ATTRIBUTES
    c_cent_basal[i]=0;
    c_cent_apical[i]=0;
}
//****************************************************************************
void delete_cell_from_basal_list(int i){
    
    for(size_t j=0; j<=14; j++){
        basal_edges[i][j]=0;
    }
    
    //ATTRIBUTES
    c_V0[i]=0;
    c_kV[i]=0;
    c_alpha[i]=0;
    c_beta[i]=0;
    
    c_freeId.add(i);
}
//****************************************************************************
void dissolve_lateral(int i){

    //IDENTIFIES ALL 4 FACETS
    int fac1, fac2, fac3, fac4;
    fac1=e_lateral1[i];
    fac2=e_lateral2[i];
    fac3=e_lateral3[i];
    fac4=e_lateral4[i];
    
    //IDENTIFIES PASSIVE VERTEX
    int pass_vert=f[fac1][3];
    
    //DISSOLVES ALL 4 FACETS
    dissolve_facet(fac1);
    dissolve_facet(fac2);
    dissolve_facet(fac3);
    dissolve_facet(fac4);
    
    //DISSOLVES PASSIVE VERTEX
    dissolve_vertex_pass(pass_vert);
    
}
//****************************************************************************
