#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//******************************OUTPUT****************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

struct Point {
    Point() {x = 0; y = 0; z = 0;};
    Point(double _x, double _y, double _z) {x = _x; y = _y; z = _z; };
    double x, y, z;
};

void write_vtk(int kcount) {
    //Fix Periodic Boundary
    double dxdydz[3] = {0., 0., 0.};
    std::vector<Point> vert_temp;
    vert_temp.resize(Nv);
    std::vector<std::vector<int>> poly_basal, poly_apical, poly_lateral;
    poly_basal.resize(Nc);
    poly_apical.resize(Nc);
    poly_lateral.resize(Ne);
    for (int i = 0; i < Nv; i++) {
        if (v[i+1][0] > 0.5) {
            vert_temp[i].x = v[i+1][1];
            vert_temp[i].y = v[i+1][2];
            vert_temp[i].z = v[i+1][3];
        }
    }
    for (int m = 0; m < Nc; m++) {
        //copy basal_verticies to poly_basal
        poly_basal[m].resize(basal_vertices[m+1][2]);
        for (int i = 3; i <= 2+basal_vertices[m+1][2]; i++) {
            poly_basal[m][i-3] = basal_vertices[m+1][i]-1;
        }

        //copy v_partner[basal_verticies] to poly_apical
        poly_apical[m].resize(basal_vertices[m+1][2]);
        for (int i = 3; i <= 2+basal_vertices[m+1][2]; i++) {
            poly_apical[m][i-3] = v_partner[basal_vertices[m+1][i]]-1;
        }

        int vert_ref_id = poly_basal[m][0]+1;
        bool new_vert;
        for (int i = 0; i < poly_basal[m].size(); i++) {
            //basal side
            new_vert = false;
            int vert_id = poly_basal[m][i]+1;
            torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
            for (int j = 0; j < 3; j++) {
                if (fabs(dxdydz[j]) > 1e-6) {
                    new_vert = true;
                    break;
                }
            }
            if (new_vert) {
                vert_temp.push_back(Point(v[vert_id][1] + dxdydz[0], v[vert_id][2] + dxdydz[1], v[vert_id][3] + dxdydz[2]));
                poly_basal[m][i] = vert_temp.size()-1;
            }

            //apical side
            new_vert = false;
            vert_id = poly_apical[m][i]+1;
            torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
            for (int j = 0; j < 3; j++) {
                if (fabs(dxdydz[j]) > 1e-6) {
                    new_vert = true;
                    break;
                }
            }
            if (new_vert) {
                vert_temp.push_back(Point(v[vert_id][1] + dxdydz[0], v[vert_id][2] + dxdydz[1], v[vert_id][3] + dxdydz[2]));
                poly_apical[m][i] = vert_temp.size()-1;
            }
        }
    }
    for (int m = 0; m < Ne; m++) {
        //create poly_lateral
        poly_lateral[m].resize(4);
        poly_lateral[m][0] = e[m+1][1]-1;
        poly_lateral[m][1] = e[m+1][2]-1;
        poly_lateral[m][2] = v_partner[e[m+1][2]]-1;
        poly_lateral[m][3] = v_partner[e[m+1][1]]-1;

        int vert_ref_id = poly_lateral[m][0]+1;
        bool new_vert;
        for (int i = 1; i < 4; i++) {
            //lateral side
            new_vert = false;
            int vert_id = poly_lateral[m][i]+1;
            torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
            for (int j = 0; j < 3; j++) {
                if (fabs(dxdydz[j]) > 1e-6) {
                    new_vert = true;
                    break;
                }
            }
            if (new_vert) {
                vert_temp.push_back(Point(v[vert_id][1] + dxdydz[0], v[vert_id][2] + dxdydz[1], v[vert_id][3] + dxdydz[2]));
                poly_lateral[m][i] = vert_temp.size()-1;
            }
        }

    }

    const int filenameSize = 200;
    char filename[filenameSize];

    int cell_size = 0;
    for (int i = 0; i < poly_basal.size(); i++) {
        cell_size += poly_basal[i].size() + 1;
    }
    
    //BASAL SIDES
    snprintf(filename, sizeof(char) * filenameSize, "./output/outBasal_%d.vtk", kcount);
    FILE *file1; 
    file1 = fopen(filename, "w");
    fprintf(file1, "# vtk DataFile Version 3.0\n");
    fprintf(file1, "Cell Mesh\n");
    fprintf(file1, "ASCII\n\n");
    fprintf(file1, "DATASET POLYDATA\n");
    fprintf(file1, "POINTS %d double\n", vert_temp.size());
    //Write out verts
    for(int i = 0; i < vert_temp.size(); i++) {
        fprintf(file1, "%f %f %f\n", vert_temp[i].x, vert_temp[i].y, vert_temp[i].z);
    }
    fprintf(file1, "\nPOLYGONS\t%d\t%d\n", poly_basal.size(), cell_size);
    for (int i = 0; i < poly_basal.size(); i++) {
        fprintf(file1, "%d ", poly_basal[i].size());
        for(size_t j = 0; j < poly_basal[i].size(); j++) {
            fprintf(file1, "%d ", poly_basal[i][j]);
        }
        fprintf(file1, "\n");
    }
    fprintf(file1, "\nCELL_DATA\t%d\n", poly_basal.size());
    fprintf(file1, "FIELD Cell_properties 1\n");
    fprintf(file1, "Cell_sides 1 %d double\n", poly_basal.size());
    for (int i = 0; i < poly_basal.size(); i++) {
        fprintf(file1, "%d\n", poly_basal[i].size());
    }
    fclose(file1);
    
    //APICAL SIDES
    snprintf(filename, sizeof(char) * filenameSize, "./output/outApical_%d.vtk", kcount);
    FILE *file2;
    file2 = fopen(filename, "w");
    fprintf(file2, "# vtk DataFile Version 3.0\n");
    fprintf(file2, "Cell Mesh\n");
    fprintf(file2, "ASCII\n\n");
    fprintf(file2, "DATASET POLYDATA\n");
    fprintf(file2, "POINTS %d double\n", vert_temp.size());
    //Write out verts
    for(int i = 0; i < vert_temp.size(); i++) {
        fprintf(file2, "%f %f %f\n", vert_temp[i].x, vert_temp[i].y, vert_temp[i].z);
    }
    fprintf(file2, "\nPOLYGONS\t%d\t%d\n", poly_apical.size(), cell_size);
    for (int i = 0; i < poly_apical.size(); i++) {
        fprintf(file2, "%d ", poly_apical[i].size());
        for(size_t j = 0; j < poly_apical[i].size(); j++) {
            fprintf(file2, "%d ", poly_apical[i][j]);
        }
        fprintf(file2, "\n");
    }
    fprintf(file2, "\nCELL_DATA\t%d\n", poly_apical.size());
    fprintf(file2, "FIELD Cell_properties 1\n");
    fprintf(file2, "Cell_sides 1 %d double\n", poly_apical.size());
    for (int i = 0; i < poly_apical.size(); i++) {
        fprintf(file2, "%d\n", poly_apical[i].size());
    }
    fclose(file2);

    //LATERAL SIDES
    snprintf(filename, sizeof(char) * filenameSize, "./output/outLateral_%d.vtk", kcount);
    FILE *file3;
    file3 = fopen(filename, "w");
    fprintf(file3, "# vtk DataFile Version 3.0\n");
    fprintf(file3, "Cell Mesh\n");
    fprintf(file3, "ASCII\n\n");
    fprintf(file3, "DATASET POLYDATA\n");
    fprintf(file3, "POINTS %d double\n", vert_temp.size());
    //Write out verts
    for(int i = 0; i < vert_temp.size(); i++) {
        fprintf(file3, "%f %f %f\n", vert_temp[i].x, vert_temp[i].y, vert_temp[i].z);
    }
    fprintf(file3, "\nPOLYGONS\t%d\t%d\n", poly_lateral.size(), 5*poly_lateral.size());
    for(int i = 0; i < poly_lateral.size(); i++) {
        fprintf(file3, "4 %d %d %d %d\n", poly_lateral[i][0], poly_lateral[i][1], poly_lateral[i][2], poly_lateral[i][3]);
    }
    fclose(file3);
}

double out_structureB(int intTime){
    
    double *dxdydz = new double[3];
    dxdydz[0]=0;dxdydz[1]=0;dxdydz[2]=0;
    
    char fileName[50];
    snprintf(fileName, sizeof(char) * 50, "./output/structureB_%d.x", intTime);
    FILE *file1; file1 = fopen(fileName, "wt");
    snprintf(fileName, sizeof(char) * 50, "./output/structureB_%d.y", intTime);
    FILE *file2; file2 = fopen(fileName, "wt");
    snprintf(fileName, sizeof(char) * 50, "./output/structureB_%d.z", intTime);
    FILE *file3; file3 = fopen(fileName, "wt");
    
    int vert_ref_id, vert_id;
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            fprintf(file1, "%d ", basal_vertices[i][2]);
            fprintf(file2, "%d ", basal_vertices[i][2]);
            fprintf(file3, "%d ", basal_vertices[i][2]);
        
            vert_ref_id=basal_vertices[i][3];
            fprintf(file1, "%g ", v[vert_ref_id][1]);
            fprintf(file2, "%g ", v[vert_ref_id][2]);
            fprintf(file3, "%g ", v[vert_ref_id][3]);
        
            for(size_t j = 4; j <= 2+basal_edges[i][2]; j++){
                vert_id=basal_vertices[i][j];
                torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
                fprintf(file1, "%g ", v[vert_id][1] + dxdydz[0]);
                fprintf(file2, "%g ", v[vert_id][2] + dxdydz[1]);
                fprintf(file3, "%g ", v[vert_id][3] + dxdydz[2]);
            }
        
            for(size_t j = basal_edges[i][2]; j <= 15; j++){
                fprintf(file1, "0 ");
                fprintf(file2, "0 ");
                fprintf(file3, "0 ");
            }
        
            fprintf(file1, "\n");
            fprintf(file2, "\n");
            fprintf(file3, "\n");
        }
    }
    fclose(file1);
    fclose(file2);
    fclose(file3);
    
    delete []dxdydz;
    
    return 0;
}
//****************************************************************************
double out_structureA(int intTime){
    
    double *dxdydz = new double[3];
    dxdydz[0]=0;dxdydz[1]=0;dxdydz[2]=0;
    
    char fileName[50];
    snprintf(fileName, sizeof(char) * 50, "./output/structureA_%d.x", intTime);
    FILE *file1; file1 = fopen(fileName, "wt");
    snprintf(fileName, sizeof(char) * 50, "./output/structureA_%d.y", intTime);
    FILE *file2; file2 = fopen(fileName, "wt");
    snprintf(fileName, sizeof(char) * 50, "./output/structureA_%d.z", intTime);
    FILE *file3; file3 = fopen(fileName, "wt");
    
    int vert_ref_id, vert_id;
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            fprintf(file1, "%d ", basal_vertices[i][2]);
            fprintf(file2, "%d ", basal_vertices[i][2]);
            fprintf(file3, "%d ", basal_vertices[i][2]);
        
            vert_ref_id=v_partner[basal_vertices[i][3]];
            fprintf(file1, "%g ", v[vert_ref_id][1]);
            fprintf(file2, "%g ", v[vert_ref_id][2]);
            fprintf(file3, "%g ", v[vert_ref_id][3]);
        
            for(size_t j = 4; j <= 2+basal_edges[i][2]; j++){
                vert_id=v_partner[basal_vertices[i][j]];
                torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
                fprintf(file1, "%g ", v[vert_id][1] + dxdydz[0]);
                fprintf(file2, "%g ", v[vert_id][2] + dxdydz[1]);
                fprintf(file3, "%g ", v[vert_id][3] + dxdydz[2]);
            }
        
            for(size_t j = basal_edges[i][2]; j <= 15; j++){
                fprintf(file1, "0 ");
                fprintf(file2, "0 ");
                fprintf(file3, "0 ");
            }
        
            fprintf(file1, "\n");
            fprintf(file2, "\n");
            fprintf(file3, "\n");
        }
    }
    fclose(file1);
    fclose(file2);
    fclose(file3);
    
    delete []dxdydz;
    
    return 0;
}
//****************************************************************************
double out_structureL(int intTime){
    
    double *dxdydz = new double[3];
    dxdydz[0]=0;dxdydz[1]=0;dxdydz[2]=0;
    
    char fileName[50];
    snprintf(fileName, sizeof(char) * 50, "./output/structureL_%d.x", intTime);
    FILE *file1; file1 = fopen(fileName, "wt");
    snprintf(fileName, sizeof(char) * 50, "./output/structureL_%d.y", intTime);
    FILE *file2; file2 = fopen(fileName, "wt");
    snprintf(fileName, sizeof(char) * 50, "./output/structureL_%d.z", intTime);
    FILE *file3; file3 = fopen(fileName, "wt");
    
    int vert_ref_id, vert_id;
    for(int i=1; i<=Ne; i++){
        if(e[i][0]!=0){
            fprintf(file1, "4 ");
            fprintf(file2, "4 ");
            fprintf(file3, "4 ");
        
            //1
            vert_ref_id=e[i][1];
            fprintf(file1, "%g ", v[vert_ref_id][1]);
            fprintf(file2, "%g ", v[vert_ref_id][2]);
            fprintf(file3, "%g ", v[vert_ref_id][3]);
            
            //2
            vert_id=e[i][2];
            torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
            fprintf(file1, "%g ", v[vert_id][1] + dxdydz[0]);
            fprintf(file2, "%g ", v[vert_id][2] + dxdydz[1]);
            fprintf(file3, "%g ", v[vert_id][3] + dxdydz[2]);
            
            //3
            vert_id=v_partner[e[i][2]];
            torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
            fprintf(file1, "%g ", v[vert_id][1] + dxdydz[0]);
            fprintf(file2, "%g ", v[vert_id][2] + dxdydz[1]);
            fprintf(file3, "%g ", v[vert_id][3] + dxdydz[2]);
            
            //4
            vert_id=v_partner[e[i][1]];
            torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
            fprintf(file1, "%g ", v[vert_id][1] + dxdydz[0]);
            fprintf(file2, "%g ", v[vert_id][2] + dxdydz[1]);
            fprintf(file3, "%g ", v[vert_id][3] + dxdydz[2]);
        
            for(size_t j = 4; j <= 15; j++){
                fprintf(file1, "0 ");
                fprintf(file2, "0 ");
                fprintf(file3, "0 ");
            }
        
            fprintf(file1, "\n");
            fprintf(file2, "\n");
            fprintf(file3, "\n");
        }
    }
    fclose(file1);
    fclose(file2);
    fclose(file3);
    
    delete []dxdydz;
    
    return 0;
}
//****************************************************************************
void outStructure(int intTime){
    out_structureA(intTime);
    out_structureB(intTime);
    out_structureL(intTime);
}
//****************************************************************************
void out_Vertissue3D(char *fileName){
    
    FILE *file1;
    file1 = fopen(fileName, "wt");
    int nrV=0, nrE=0, nrC=0;
    
    //ALLOCATE INTERNAL VARIABLES
    int **v_natural_id;
    int **e_natural_id;
    int **c_natural_id;
    v_natural_id  = new int*[Nv+1];
    e_natural_id  = new int*[Ne+1];
    c_natural_id  = new int*[Nc+1];
    for(size_t i=0; i < Nv+1; i++){
        v_natural_id[i]  = new int[3];
    }
    for(size_t i=0; i < Ne+1; i++){
        e_natural_id[i]  = new int[3];
    }
    for(size_t i=0; i < Nc+1; i++){
        c_natural_id[i]  = new int[3];
    }
    double *dxdydz = new double[3];
    dxdydz[0]=0;dxdydz[1]=0;dxdydz[2]=0;
    
    //RESET NATURAL IDS
    for(size_t i=1; i<=Nv; i++){
        v_natural_id[i][1]=0;
        v_natural_id[i][2]=0;
    }
    for(size_t i=1; i<=Ne; i++){
        e_natural_id[i][1]=0;
        e_natural_id[i][2]=0;
    }
    for(size_t i=1; i<=Nc; i++){
        c_natural_id[i][1]=0;
        c_natural_id[i][2]=0;
    }
    
    //VERTICES
    int cnt=0;
    for(size_t i=1; i<=Nv; i++){
        if(v[i][0]!=0 && v_type[i]==1){
            cnt++;
            v_natural_id[i][1]=cnt;
            v_natural_id[cnt][2]=i;
        }
    }
    nrV=cnt;
    
    //EDGES
    cnt=0;
    for(size_t i=1; i<=Ne; i++){
        if(e[i][0]!=0){
            cnt++;
            e_natural_id[i][1]=cnt;
            e_natural_id[cnt][2]=i;
        }
    }
    nrE=cnt;
    
    //CELLS
    cnt=0;
    for(size_t i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            cnt++;
            c_natural_id[i][1]=cnt;
            c_natural_id[cnt][2]=i;
        }
    }
    nrC=cnt;
    fprintf(file1, "%d  %d  %d\n\n", nrV, nrE, nrC);
    fprintf(file1, "%g  %g  %g\n\n", perioXYZ[0], perioXYZ[1], perioXYZ[2]);
    //*********************************************************
    //*********************************************************
    //VERTICES
    for(size_t i=1; i<=Nv; i++){
        if(v[i][0]!=0 && v_type[i]==1){
            
            torus_dx_dy_dz(dxdydz,v_partner[i],i);
            
            v_normal_vector[i][1]=v[v_partner[i]][1]+dxdydz[0]-v[i][1];
            v_normal_vector[i][2]=v[v_partner[i]][2]+dxdydz[1]-v[i][2];
            v_normal_vector[i][3]=v[v_partner[i]][3]+dxdydz[2]-v[i][3];
            
            v_height[i]=sqrt(v_normal_vector[i][1]*v_normal_vector[i][1]+v_normal_vector[i][2]*v_normal_vector[i][2]+v_normal_vector[i][3]*v_normal_vector[i][3]);
            
            v_normal_vector[i][1]/=v_height[i];
            v_normal_vector[i][2]/=v_height[i];
            v_normal_vector[i][3]/=v_height[i];
            
            fprintf(file1, "%f  %f  %f  %f  %f  %f  %f\n", v[i][1], v[i][2], v[i][3], v_normal_vector[i][1], v_normal_vector[i][2], v_normal_vector[i][3], v_height[i]);
        }
    }
    fprintf(file1, "\n");
    
    //EDGES
    for(size_t i=1; i<=Ne; i++){
        if(e[i][0]!=0){
            fprintf(file1, "%d  %d\n", v_natural_id[e[i][1]][1], v_natural_id[e[i][2]][1]);
        }
    }
    fprintf(file1, "\n");
    
    //CELLS
    for(size_t i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            fprintf(file1, "%d  ", basal_edges[i][2]);
            for(int j=3; j<=2+basal_edges[i][2]; j++){
                if(basal_edges[i][j]<0) fprintf(file1, "%d  ", -e_natural_id[abs(basal_edges[i][j])][1]);
                else fprintf(file1, "%d  ", e_natural_id[abs(basal_edges[i][j])][1]);
            }
            for(int j=12; j>basal_edges[i][2]; j--){
                fprintf(file1, "0  ");
            }
            fprintf(file1, "\n");
        }
    }
    
    //DEALLOCATE INTERNAL VARIABLES
    for(size_t i=0; i < Nv+1; i++){
        delete [] v_natural_id[i] ;
    }
    for(size_t i=0; i < Ne+1; i++){
        delete [] e_natural_id[i] ;
    }
    for(size_t i=0; i < Nc+1; i++){
        delete [] c_natural_id[i] ;
    }
    
    delete [] v_natural_id;
    delete [] e_natural_id;
    delete [] c_natural_id;
    delete []dxdydz;
    
    fclose(file1);
}
//****************************************************************************
